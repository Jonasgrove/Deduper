#!/usr/bin/env python 

#import dependencies
import os
import argparse
import sam_record_class
import glob
import shutil
from multiprocessing import Pool
import time

# assign classes 
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
      rc   n rc   n rc   n  rc   n  *third level not yet implemented

This will allow the program to operate on each of the groups of
potential duplicates in parallel, or process each of them sequentially
if memory is limited.  
'''
def make_database(data_base_dir, umi_file, file_in):

    # make database directory if not already present
    path = data_base_dir + "Database"
    
    #try:
    os.mkdir(path)

    # dictionary stores open files for writing
    # form: {index_seq_chr: openfile()}, ex. {AGCT_1: open(AGCT_1,"w")}
    db_writing_dic = {}

    # dictionary for storing valid barcode sequences 
    umi_dic = set_barcodes(umi_file)

    # open input file   
    file_in = open(file_in, "r")    # open input file

    # open file to store database file names
    meta_data = path + "/meta_data.txt"
    meta_database = open(meta_data, "w")
    
    # iterate through each line and sort it into a file 
    # based on it's UMI and chromosome
    for line in file_in:
        if line[0] != "@": 
            record = SamRecord(line)
            if record.dic_key not in db_writing_dic and record.umi in umi_dic:
                file_name = data_base_dir + "Database/" + record.dic_key
                db_writing_dic[record.dic_key] = open(file_name, "w")
                db_writing_dic[record.dic_key].write(record.line)

                # write out file name to metadata
                meta_database.write(file_name + "\n")

            elif record.umi in umi_dic:
                db_writing_dic[record.dic_key].write(record.line)
    
    # close all files held in dictionary
    for db_file in db_writing_dic.keys():
        db_writing_dic[db_file].close()

    # close input file
    file_in.close()  
    meta_database.close()             

    #except:
    print("Database is already made")

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
    parser.add_argument("-db", "--database", type=str, default="./", help="specifies the directory which the database should be created in")
    parser.add_argument("-i", "--file_in", type=str, help='specifies input file. must be SAM format')
    parser.add_argument("-u", "--UMI", type=str, help='specifies the barcodes (indexes) used in experiment')
    parser.add_argument("-p", "--parallel", default=False, type=bool, help="boolean (True or False ) which specifies if parallel processing shoule be used. default is False")
    parser.add_argument("-t", "--threads", type=int, default=8, help="specify the number of cores/threads to run multiprocessing. Only applicable when parallel=True")

    return parser.parse_args()

parseArgs = get_args()

data_base_dir = parseArgs.database
file_in = parseArgs.file_in
umi_file = parseArgs.UMI 
parallel = parseArgs.parallel
threads = parseArgs.threads

# main function
'''
main function will 
    1. first make a database by calling make_databse
    2. eliminate duplicates from each Databse file (either sequentially of in parallel)
    3. cat together all individual output files
'''
def main():
    start = time.time()
    ## make database ...O(N)
    make_database(data_base_dir, umi_file, file_in)

    ## process all files in database (in parallell) O(N / |UMI| / |chr| / 2)

    # open meta_data file and make output directory
    meta_database_file = open("./Database/meta_data.txt", "r") 

    # if parallel option is specified as True
    # process database files in parallell
    if parallel == True:
        #start = time.time()
        meta_database_list = meta_database_file.readlines()
        with Pool(threads) as p:
            p.map(find_duplicates, meta_database_list)
        end = time.time()
        print("Parallel run time: ", end - start)

    # else process files sequencially
    else:
        #start = time.time()
        for file_name in meta_database_file:
            #file_name = file_name.strip()
            find_duplicates(file_name)
        end = time.time()
        print("Sequential run time: ", end - start)

    # close meta database file
    meta_database_file.close()

    ## cat together all output_filtered files O(pretty fast..)
    read_files = glob.glob("./Database/*_*_filtered")

    with open("./filtered_14000.sam", "wb") as outfile:
        for f in read_files:
            with open(f, "rb") as infile:
                outfile.write(infile.read())

    # delete database
    path = os.path.join(data_base_dir, "Database")
    shutil.rmtree(path)
    print("Process complete. Database Deleted.")

# run program 
main()

# TO DO
'''
make option for directing output to specific directory
get rid of database directory option
add umi error correcting
add functionality for paired end sequence data
add functionality for ramdom indexes
add functionality to write out duplicates
add functionality to take duplicate with highest quality (or something)
add functionality to deal with other cigar string characters (*, D, etc)
make the databse creation multi threaded. break up the big files, sort 
them into their umis and chromosomes and then bring the files together.
'''