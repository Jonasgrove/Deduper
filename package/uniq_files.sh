#!/bin/bash

# this script makes a file of umi_chromosome prefixes in order to merge database sub-files

# arguments
# 1 = database_dir. directory which databases are made, specified by -o

ls -1 $1Database | cut -c 1,2,3,4,5,6,7,8,9,10 | sort | uniq > $1uniq_output.txt