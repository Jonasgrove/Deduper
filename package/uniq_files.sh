#!/bin/bash

# arguments
# 1 = database_dir

ls -1 $1Database | cut -c 1,2,3,4,5,6,7,8,9,10 | sort | uniq > $1uniq_output.txt