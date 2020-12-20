# Deduper

## Version 1 

* algorithm for removing PCR duplicates
* files: 
*  deduper.py: main program
** sam_record_class.py: class to store records as objects
** database_build.sh: 
* includes options for parallel processing (time efficient) or non parallel processing (memory efficient). Parallel processing processes database files in parallel. Database is made in sequential process. 


## Version 2

* algorithm for removing PCR duplicates
* files: 
** deduper.py: main program
** sam_record_class.py: class to store records as objects
** database_slurm.sh: bash script which splits input file into multiple smaller files using bash **split** command.
** uniq_files.sh: bash script which runs **uniq** command to find unique file prefixes. used to merge sorted database files. 
* includes options for parallel processing or non parallel processing. Parallel processing creates database in parallel process and removes PCR duplicates in parallel process (analogous to version 1).
