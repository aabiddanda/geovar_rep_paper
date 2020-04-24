"""
Compute the expected fraction of mutations in each geodist category.

Assumes a simple model of 2-population split with constant population size,
using the spectral representation of the Wright-Fisher transition density
from Song and Steinruecken 2012.
"""

import numpy as np
from numpy import ndarray
from scipy.integrate import quad
from scipy.special import betainc, eval_jacobi, loggamma, roots_jacobi


def __jacobi_norm_nonzero(n, theta):
    num = 2 * loggamma(n + theta)
    den = loggamma(n + 2 * theta - 1) + loggamma(n + 1)
    return np.exp(num - den) / (2 * (n + theta) - 1)


def __jacobi_norm_zero(theta):
    num = 2 * loggamma(theta)
    den = loggamma(2 * theta)
    return np.exp(num - den)


def jacobi_norm(n, theta):
    """Compute the value of <R_n, R_n>."""
    return np.where(n == 0, __jacobi_norm_zero(theta),
                    __jacobi_norm_nonzero(n, theta))


def jacobi_eigenval(n, theta):
    """Compute the eigenvalue of R_n in the WF transition density."""
    return n * (n + 2 * theta - 1)


def binom_cdf(k, n, p):
    """Compute the CDF of the binomial distribution."""
    return np.where(k == n, 1.0, betainc(n - k, k + 1, 1 - p))


def __integrand(x1, x2, k1, k2, n, c, t, theta, m, weights):
    x1_eval = eval_jacobi(m, theta - 1, theta - 1,
                          2 * np.expand_dims(x1, -1) - 1)
    x2_eval = eval_jacobi(m, theta - 1, theta - 1,
                          2 * np.expand_dims(x2, -1) - 1)
    jacobi_total = np.sum(x2_eval * x1_eval * weights, axis=-1)
    bc1 = binom_cdf(k1, n, (1 - c) * x1 + c * x2)
    bc2 = binom_cdf(k2, n, c * x1 + (1 - c) * x2)
    return bc1 * bc2 * jacobi_total


def __inner_quad(*args):
    """Deprecated."""
    theta = args[-3]
    val, err = quad(__integrand,
                    0,
                    1,
                    weight='alg',
                    wvar=(theta - 1, theta - 1),
                    args=args)
    return val


def cumulative_sfs(k1: int,
                   k2: int,
                   n: int,
                   c: float,
                   t: float,
                   theta: float,
                   m_max: int,
                   order: int = None) -> float:
    """Compute the CDF of the expected joint site frequency spectrum.

    Parameters
    ----------
    k1 : int
        the allele count in population 1
    k2 : int
        the allele count in population 2
    n : int
        the sample size in each population
    c : float
        the admixture proportion in the present
    t : float
        the population-scaled split time
    theta : float
        the population-scaled mutation rate
    m_max : int
        the largest basis function to compute
    order : int, optional
        the quadrature order (default=2n)

    Returns
    -------
    float

    """
    m = np.arange(m_max + 1)
    lambdas = jacobi_eigenval(m, theta)
    norms = jacobi_norm(m, theta)
    # The jacobi_norm(0, theta) scales the density to integrate to 1
    weights = np.exp(-2 * lambdas * t) / norms / jacobi_norm(0, theta)

    if order is None:
        order = 2 * n
    x, w = roots_jacobi(order, theta - 1, theta - 1)
    x = (x + 1) / 2
    w /= 2**(2 * theta - 1)
    W = w * w[:, np.newaxis]

    evals = __integrand(x, x[:, np.newaxis], k1, k2, n, c, t, theta, m,
                        weights)
    return np.sum(W * evals)


def __cumulative_sfs(k1: int, k2: int, n: int, c: float, t: float,
                     theta: float, m_max: int) -> float:
    """Deprecated."""
    m = np.arange(m_max + 1)
    lambdas = jacobi_eigenval(m, theta)
    norms = jacobi_norm(m, theta)
    # The jacobi_norm(0, theta) scales the density to integrate to 1
    weights = np.exp(-2 * lambdas * t) / norms / jacobi_norm(0, theta)
    # Parameters to pass through to integrand
    # With each call to quad, will prepend a dummy variable to serve as x1, x2
    args = (k1, k2, n, c, t, theta, m, weights)

    val_outer, err_outer = quad(__inner_quad,
                                0,
                                1,
                                weight='alg',
                                wvar=(theta - 1, theta - 1),
                                args=args)
    return val_outer


def binned_joint_sfs(K1, K2, *args) -> ndarray:
    """Compute the joint SFS for allele count bins.

    Parameters
    ----------
    K1 : array_like
        The upper bounds (inclusive) of allele count bins in population 1
    K2 : array_like
        The upper bounds (inclusive) of allele count bins in population 2
    args :
        Arguments to pass to cumulative_sfs

    Returns
    -------
    ndarray

    """
    csfs = np.zeros((len(K1), len(K2)))
    for i, k1 in enumerate(K1):
        for j, k2 in enumerate(K2):
            csfs[i, j] = cumulative_sfs(k1, k2, *args)
    diff1 = np.diff(csfs, axis=0, prepend=0)
    return np.diff(diff1, axis=1, prepend=0)


def __sfs2geodist(sfs):
    gc_counts = np.zeros((3, 3))
    gc_counts[:3, :3] = sfs[:3, :3]
    # UU
    gc_counts[0, 0] += sfs[4, 4]
    # UR
    gc_counts[0, 1] += sfs[4, 3]
    # RU
    gc_counts[1, 0] += sfs[3, 4]
    # UC
    gc_counts[0, 2] += np.sum([sfs[4, 1:3], sfs[0, 3:]])
    # CU
    gc_counts[2, 0] += np.sum([sfs[1:3, 4], sfs[3:, 0]])
    # RR
    gc_counts[1, 1] += sfs[3, 3]
    # NOTE: the last two lines make an approximation re the minor allele
    # that will be canceled by the symmetry.
    # RC.
    gc_counts[1, 2] += np.sum(sfs[3, 1:3])
    # CR
    gc_counts[2, 1] += np.sum(sfs[1:3, 3])
    return gc_counts


def geodist_counts(k_rare: int,
                   n: int,
                   c: float,
                   t: float,
                   theta: float,
                   m_max: int = None,
                   symmetrize: bool = True,
                   normalize: bool = True) -> ndarray:
    """Compute the relative counts of mutations in different geodist categories.

    Parameters
    ----------
    k_rare : int
        the allele count cutoff between rare mutations (min(k, n-k) <= k_rare)
        and common mutations (min(k, n-k) > k_rare)
    n : int
        the sample size
    c : float
        the symmetric admixture coefficient between populations in the present
    t : float
        the population-scaled split time
    theta : float
        the population-scaled mutation rate
    m_max : int
        the number of basis functions to evaluate (default = n)
    symmetrize : bool
        if True, force the populations to be symmetric to make up for numerical
        error (default = True)
    normalize : bool
        if True, set the UU category to zero and normalize the remaining counts
        to sum to one (default = True)

    Returns
    -------
    ndarray
        3x3 matrix of counts of two-population allele frequency categories.
        The layout is:
            UU  UR  UC
            RU  RR  RC
            CU  CR  CC

    """
    if not m_max:
        m_max = n
    K = [0, k_rare, n - k_rare, n - 1, n]
    sfs = binned_joint_sfs(K, K, n, c, t, theta, m_max)
    gc_counts = __sfs2geodist(sfs)

    if symmetrize:
        gc_counts = (gc_counts + gc_counts.T) / 2

    if normalize:
        gc_counts[0, 0] = 0
        gc_counts /= np.sum(gc_counts)

    return gc_counts


def test_run():
    # TODO: write some real tests
    """Run on some test parameters. Print the output."""
    n = 100
    c = 0.01
    t = 0.01
    theta = 1e-3

    k1 = 4
    k2 = 10
    m_max = n
    I1 = cumulative_sfs(k1, k2, n, c, t, theta, m_max)
    I2 = __cumulative_sfs(k1, k2, n, c, t, theta, m_max)
    print(I1)
    print(I2)
    print(np.isclose(I1, I2))

    k_rare = 4
    gdcounts = geodist_counts(k_rare, n, c, t, theta, normalize=True)
    print(gdcounts)


if __name__ == "__main__":
    test_run()
