#!/bin/bash

# arguments
# 1 = database_directory
# 2 = number of lines to split by
# 3 = file to split 

# make build database directory
mkdir $1/database_build

# split file into subfiles to sort in parallell
split -l $2 $3 $1/database_build/pref_ 

# make file of database file names
ls -1 $1/database_build > $1/metadata_build.txt