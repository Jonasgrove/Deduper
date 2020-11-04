#!/usr/bin/env python 

#import dependencies
import os
import argparse
import sam_record_class

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
potential duplicates in parallel.  
'''
def make_database(data_base_dir, umi_file, file_in):
    # make database directory
    path = data_base_dir + "/Database"
    try:
        os.mkdir(path)
    except:
        print("Database is already present")

    # dictionary stores open files for writing
    # form: {index_seq_chr: openfile()}, ex. {AGCT_1: open(AGCT_1,"w")}
    db_writing_dic = {}

    # dictionary for storing valid barcode sequences
    umi_dic = set_barcodes(umi_file)
                      
    file_in = open(file_in, "r")    # open input file
    
    # iterate through each line and sort it into a file 
    # based on it's UMI and chromosome
    for line in file_in:
        if line[0] != "@": 
            record = SamRecord(line)
            if record.dic_key not in db_writing_dic and record.umi in umi_dic:
                file_name = "./Database/" + record.dic_key
                db_writing_dic[record.dic_key] = open(file_name, "w")
                db_writing_dic[record.dic_key].write(record.line)
            else:
                db_writing_dic[record.dic_key].write(record.line)
    
    # close all files held in dictionary
    for db_file in db_writing_dic.keys():
        db_writing_dic[db_file].close()

    # close input file
    file_in.close()               

    return

# function which finds duplicates
'''
iterate through a file that contains all the reads of particular UMI 
and a particular chromosome and store the adjusted position into 
a dictionar. If the position is already in the dictionary, continue
else, write the record out to the file.
'''
def find_duplicates(dup_file):

    # initializations
    dup_fh = open(dup_file, "r") # open file
    dup_dic = {}                 # init empty dictionary to be {"adjusted position" : True}
    out_file_name = dup_file + "_filtered"
    out_fh = open(out_file_name)  # file to write out unique reads

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
    parser = argparse.ArgumentParser(description='deduplicate read remover')
    parser.add_argument("-db", "--database", type=str, help="specifies the directory which the database should be created in")
    parser.add_argument("-i", "--file_in", type=str, help='specifies input file. must be SAM format')
    parser.add_argument("-u", "--UMI", type=str, help='specifies the barcodes (indexes) used in experiment')


    return parser.parse_args()

parseArgs = get_args()

data_base_dir = parseArgs.database
file_in = parseArgs.file_in
umi_file = parseArgs.UMI 

# main function
'''
main function will 
    1. first make a database by calling make_databse
    2. eliminate duplicates from each Databse file
    3. cat together all individual output files
'''
#def main():

    #make database ...O(N)?



    # process all files in database (in parallell) O(N / |UMI| / |chr| / 2)?

    # cat together all output_filtered files O(?)



make_database(data_base_dir, umi_file, file_in)


