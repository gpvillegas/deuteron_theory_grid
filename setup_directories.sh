#!/usr/bin/bash

# example to prepare a calculation to run in directory
# /Users/boeglinw/Documents/boeglin.1/user/deut/misak/grid/

python ../edepngo4_v1/python/calc_crosecs.py \
    --base_dir /Users/boeglinw/Documents/boeglin.1/ \
        --edpngo_dir /user/deut/misak/edepngo4_v1/ \
# needs to be adjusted to each case
            --local_dir /user/deut/misak/grid/ \
                --result_dir /results/ -S 



