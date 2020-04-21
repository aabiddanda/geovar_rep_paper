#!/bin/python3



# Including all of the different snakefiles
include: 'snakefiles/geodist.snake'
include: 'snakefiles/subsets.snake'
include: 'snakefiles/plotting.snake'





rule download_data:
  shell:
    """
      wget <data.tar.gz>
      tar -xvf <data.tar.gz>
    """
  
    
    
rule download_data_w_geodist:
  shell:
    """
      wget <data_geodist.tar.gz>
      tar -xvf <data_geodist.tar.gz>
    """

    
# rule total_clean:
#   shell:
#     """
#     """
    