#!/usr/bin/env python 

#import dependencies
import os
import argparse
import sam_record_class
import glob
import shutil
from multiprocessing import Pool
import time
import subprocess # didnt actually use

# assign class to store sam records 
SamRecord = sam_record_class.SamRecord

# make dictionary of valid barcodes
def set_barcodes(index_file):

    # initialize dictionary for indexes {"AGTCG": True}
    index_dic = {}  
    index_fh = open(index_file, "r")
    for index in index_fh:
        index = index.strip()
        index_dic[index] = True
    
    return index_dic

## make database of umi_chromsome sorted files in parallel

'''
idea is to build database in a parallel process:
    1. first the file is broken into sub-files using split
        see: database_build.sh
        use linux split in bash script and output file of file names (ls)
        see: uniq_files.sh

    2. each sub-file is sorted by UMI_chrom
        see: make_database()
        
    3. all files are merged together based on matching UMI_chromosme identifier after being sorted
'''

def parallel_database(data_base_dir, umi_file, file_in, threads, size):

    '''
    splt up files 
    '''
    # determine split size
    split_size = int(size/10)
    
    # split file into files of x percent of reads
    # arguments
    # 1 = data_base_dir | 2 = .10 of read number | 3 = file_in
    build_command = "./database_build.sh " + data_base_dir + " " + str(split_size) + " " + file_in
    os.system(build_command)

    # get all sub-file names and make_database arguments
    build_files = open(data_base_dir + "/metadata_build.txt", "r")
    build_args = [(data_base_dir, umi_file, file_name.strip(), True, i) for i, file_name in enumerate(build_files)]


    '''
    sort files in parallel process 
    '''
    with Pool(threads) as p:
        p.starmap(make_database, build_args)

    '''
    use bash script to get all unique prefixes
    '''
    uniq_arg = "./uniq_files.sh " + data_base_dir
    os.system(uniq_arg)

    # make list of all prefix files
    prefix_file = open(data_base_dir + "uniq_output.txt", "r")
    

    '''
    merge sorted subfiles based on matching UMI_chromosome
    '''
    # iterate through prefixes and merge files based on prefix
    for prefix in prefix_file:
        prefix = prefix.strip()

        # get list of filenames  
        glob_arg = data_base_dir + "Database/" + prefix + "_*"    
        read_files = glob.glob(glob_arg)

        # merge all files with same UMI_chromosome_*
        with open(data_base_dir + "Database/" + prefix, "wb") as outfile:
            for f in read_files:
                with open(f, "rb") as infile:
                    outfile.write(infile.read())
    
    '''
    Delete database_build and sorted subfiles. 
    update metadata.txt file for threading deuplication
    '''
    # delete sub-files
    rm_command = "rm " + data_base_dir + "Database/*_*_*" 
    os.system(rm_command)

    # adjust metadata file
    metafile_adjust = "ls -1 " + data_base_dir + "Database > " + data_base_dir + "metadata.txt"
    os.system(metafile_adjust)

    # delete build directory
    path = os.path.join(data_base_dir, "database_build")
    shutil.rmtree(path)

    return




# make database directory
'''
This function will iterate through the file and place reads in database
structure: 
                Database
                ___|___
               |       |
            UMI_1       UMI_n   
           ___|___      ___|___
          |       |    |       |
        CH_1   CH_n  CH_1    CH_n
        _|_    _|_    _|_     _|_
       |   |  |   |  |   |   |   |
      rc   n rc   n rc   n  rc   n  *leaf level not yet implemented

This will allow the program to operate on each of the groups of
potential duplicates in parallel, or process each of them sequentially
if memory is limited.  
'''
def make_database(data_base_dir, umi_file, file_in, parallel, identifyer):

    # make database directory if not already present
    path = data_base_dir + "Database"
    
    try:
        os.mkdir(path)
    except:
        pass

    # dictionary stores open files for writing
    # form: {index_seq_chr: openfile()}, ex. {AGCT_1: open(AGCT_1,"w")}
    db_writing_dic = {}

    # dictionary for storing valid barcode sequences 
    umi_dic = set_barcodes(umi_file)

    # open input file
    if parallel == True:
        file_in = data_base_dir + "/database_build/" + file_in
        file_in = open(file_in, "r")
    else:    
        file_in = open(file_in, "r")    # open input file

    # open file to store database file names
    meta_data = data_base_dir + "/metadata.txt"
    meta_database = open(meta_data, "w")
    
    # iterate through each line and sort it into a file 
    # based on it's UMI and chromosome
    for line in file_in:
        if line[0] != "@": 
            record = SamRecord(line)

            # if parallel: add unique string suffix to key
            if parallel == True:
                record.dic_key = record.dic_key + "_" + str(identifyer)

            # if key has not been seen yet, and the umi is valid
            if record.dic_key not in db_writing_dic and record.umi in umi_dic:
                file_name = data_base_dir + "Database/" + record.dic_key
                db_writing_dic[record.dic_key] = open(file_name, "w")
                db_writing_dic[record.dic_key].write(record.line)

                # write out file name to metadata
                meta_database.write(record.dic_key + "\n")

            elif record.umi in umi_dic:
                db_writing_dic[record.dic_key].write(record.line)
    
    # close all files held in dictionary
    for db_file in db_writing_dic.keys():
        db_writing_dic[db_file].close()

    # close input file
    file_in.close()  
    meta_database.close()             

    return 

# function which finds duplicates
'''
iterate through a file that contains all the reads of particular UMI 
and a particular chromosome and store the adjusted position into 
a dictionar. If the position is already in the dictionary, continue
else, write the record out to the file.
'''
def find_duplicates(dup_file):

    #strip \n
    dup_file = dup_file.strip()

    # initializations
    dup_fh = open(dup_file, "r") # open file
    dup_dic = {}                 # init empty dictionary to be {"adjusted position" : True}
    out_file_name = dup_file + "_filtered"
    out_fh = open(out_file_name, "w")  # file to write out unique reads

    # iterate thorugh database file (umi(n)_chr(n))
    for line in dup_fh:
        record = SamRecord(line)        # make record obj

        # check to see if same position was alread seen
        if record.position_adj not in dup_dic:
            dup_dic[record.position_adj] = True
            out_fh.write(record.line)

    # close files
    dup_fh.close()
    out_fh.close()

    return 


# import arguments from command line
def get_args():
    parser = argparse.ArgumentParser(description='PCR deduplicate read remover. Requires sam file and UMI file as input')
    parser.add_argument("-db", "--database", type=str, default="./", help="specifies the directory which the database should be created in. default = ./ ")
    parser.add_argument("-i", "--file_in", type=str, help='specifies input file. must be SAM format')
    parser.add_argument("-u", "--UMI", type=str, help='specifies the barcodes (indexes) used in experiment')
    parser.add_argument("-p", "--parallel", default=False, type=bool, help="boolean (True or False ) which specifies if parallel processing shoule be used. default is False")
    parser.add_argument("-t", "--threads", type=int, default=8, help="specify the number of cores/threads to run multiprocessing. Only applicable when parallel=True")
    parser.add_argument("-s", "--size", type=int, default=None, help="specify size of the file, so that multiprocessing can break acordingly. only applicable when parallel=True")
    parser.add_argument("-o", "--output", type=str, default="./deduped.sam", help="specify directory and name of output file. default = ./deduped.sam")
    parser.add_argument("-b", "--break_factor", type=float, default= .10, help="specify the percentage use to break up file size. default is .10 (i.e subfiles are 10 percent the size of original")

    return parser.parse_args()

parseArgs = get_args()

data_base_dir = parseArgs.database
file_in = parseArgs.file_in
umi_file = parseArgs.UMI 
parallel = parseArgs.parallel
threads = parseArgs.threads
size = parseArgs.size
output = parseArgs.output
break_factor = parseArgs.break_factor

# main function
'''
main function will 
    1. first make a database by calling either:
        call prallel_databases which will call make_databases in parallell (-p True)
        or
        call make_databases on whole file (-p False)
    2. eliminate duplicates from each Database file (either sequentially of in parallel)
    3. merge together all output sub-files
'''
def main():
    start = time.time()

    ## make database ...O(N)
    '''
    if parallel is chosen: the function parallel_datbases is called which,
        1. breaks input randomly into sub-files based on the size of the file (|sub-file| = |file|*.10 )
        2. call make_database using all sub-files as arguments, in parallel process
        3. join files back together based on UMI_chromosome identifyer
    else: 
        iterate through file using make_database function
    '''
    if parallel == True:
        parallel_database(data_base_dir, umi_file, file_in, threads, size)
    else:
        make_database(data_base_dir, umi_file, file_in, parallel, 0)

    ## process all files in database (in parallell) O(N / |UMI| / |chr|)

    # open meta_data file and make output directory
    meta_database_file = open(data_base_dir + "metadata.txt", "r") 

    # if parallel option is specified as True
    # process database files in parallell
    if parallel == True:
    
        # get database files
        meta_database_list = [data_base_dir + "Database/" + file_name.strip() for file_name in meta_database_file]


        with Pool(threads) as p:
            p.map(find_duplicates, meta_database_list)
        end = time.time()
        print("Parallel run time: ", end - start)

    # else process files sequencially
    else:

        for file_name in meta_database_file:
            #file_name = file_name.strip()
            find_duplicates(data_base_dir + "Database/" + file_name.strip())
        end = time.time()
        print("Sequential run time: ", end - start)

    # close meta database file
    meta_database_file.close()

    ## cat together all output_filtered files O(pretty fast..)
    read_files = glob.glob("./Database/*_*_filtered")

    with open(output, "wb") as outfile:
        for f in read_files:
            with open(f, "rb") as infile:
                outfile.write(infile.read())

    # delete database
    path = os.path.join(data_base_dir, "Database")
    shutil.rmtree(path)

    # delete metafiles
    command = "rm " + data_base_dir + "metadata.txt " + data_base_dir + "metadata_build.txt " + data_base_dir + "uniq_output.txt"
    os.system(command)
    print("Process complete. Database Deleted.")


    return
# run program 
main()

# TO DO
'''
(done) make option for directing output to specific directory
(cancelled) get rid of database directory option
add umi error correcting
add functionality for paired end sequence data
add functionality for ramdom indexes
(eh) add functionality to write out duplicates
add functionality to take duplicate with highest quality (or something)
add functionality to deal with other cigar string characters (*, D, etc)
(done) make the databse creation multi threaded. break up the big files, sort 
them into their umis and chromosomes and then bring the files together.
'''
