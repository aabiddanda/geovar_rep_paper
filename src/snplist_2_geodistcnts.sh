#/bin/bash

# Takes a reference file of category assignments and list of snps by
# chromosome & position and returns the counts of each category

# The groupbyfield variable below directs the script as to which column to use
# for category assignments

if [  $# -le 1 ]
then
    echo -e "\n Usage: $0 [geodist.txt.gz] [snplist file] "
    exit 1
fi

GZFILE=$1
SNPLISTFILE=$2
GROUPBYFIELD=5

awk 'NR==FNR{a[$1,$2]++ ; next} ($1,$2) in a' $SNPLISTFILE  <(zcat $GZFILE) | awk -v grpby=$GROUPBYFIELD '{cnt[$grpby]++}END{ for(i in cnt) {print i,cnt[i]}}' | sort -n -k1,1

