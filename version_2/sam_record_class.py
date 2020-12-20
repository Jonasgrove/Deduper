#!/usr/bin/env python 

import statistics
import re

# sam file record class
'''
This is a class which serves the purpose of storing 
all the desired information from a sam file alignment
line. 
'''

class SamRecord:

    def __init__(self, sam_line):
        
        self.line = sam_line
        self.line_list = sam_line.split()
        self.umi = self.get_umi()
        self.chr = self.line_list[2]
        self.bit_flag = int(self.line_list[1])
        self.cigar = self.line_list[5]
        self.left_most = int(self.line_list[3])
        self.dic_key = self.umi + "_" + self.chr
        self.qual = self.get_qual()
        self.rev_comp = (self.bit_flag & 16) == 16
        self.position_adj = str(self.get_position()) + "_" + str(self.rev_comp)
        

    # pulls UMI off of header
    def get_umi(self):     

        header = self.line_list[0]
        umi = header.split(":")[-1]

        return umi
    
    # take average quality of whole read
    def get_qual(self):
        qual = self.line_list[10]
        qual_list = [ord(asci) - 33 for asci in qual]

        return statistics.mean(qual_list)

    # method to get cigar sting adjusted starting position
    def get_position(self):
        
        # calculate soft clipping if reverse complemented
        
        if self.rev_comp == True:
            # remove S on left side if present
            # assuming no soft clipping beyond 2 digits
            for i in range(3):  

                # if there is soft clipping on left side, trim it
                if self.cigar[i] == "S":
                    cigar_trimmed = self.cigar[i+1:]
                    break

                # else use whole sting
                else:
                    cigar_trimmed = self.cigar
            #print(cigar_trimmed)
            # include   M, N, S, I, D

            # now iterate through trimmed cigar string
            # turn numbers into ints and add them together
            number_string = ""
            align_len = 0
            for i in range(len(cigar_trimmed)):    
                
                # if i is a digit          
                if cigar_trimmed[i] != "M" and cigar_trimmed[i] != "N" and cigar_trimmed[i] != "S" and cigar_trimmed[i] != "I" and cigar_trimmed[i] != "D":
                    #print(number_string)
                    number_string += cigar_trimmed[i]
               
                elif cigar_trimmed[i] == "I":
                    #print("elif")
                    #print(number_string)
                    number_string = ""
               

                # else i is a character, add number_string to total
                else:
                    #print(number_string)
                    align_len += int(number_string)
                    number_string = ""

            cigar_adjusted_position = self.left_most + align_len    


        # calculate soft clipping if NOT rev comp
        else:
            if "S" in self.cigar[0:3]:
                clipping = int(self.cigar.split("S")[0])
            else:
                clipping = 0

            cigar_adjusted_position = self.left_most - clipping

        return cigar_adjusted_position


# testing 
'''
sam_record = "_:_:_:_:_:_:70:CTGTTCAC	16	5	1000	36	100M	*	0	0	*	*"

sam_record_obj = SamRecord(sam_record)

print(sam_record_obj.line)
print(sam_record_obj.line_list) 
print(sam_record_obj.umi)
print(sam_record_obj.chr)
print(sam_record_obj.bit_flag) 
print(sam_record_obj.cigar)
print(sam_record_obj.left_most) 
print(sam_record_obj.dic_key)
print(sam_record_obj.qual)
print(sam_record_obj.rev_comp)
print(sam_record_obj.position_adj)

sam_record = "_:_:_:_:_:_:72:CTGTTCAC	16	5	1001	36	4S99M	*	0	0	*	*"
sam_record_obj = SamRecord(sam_record)

print(sam_record_obj.line)
print(sam_record_obj.line_list) 
print(sam_record_obj.umi)
print(sam_record_obj.chr)
print(sam_record_obj.bit_flag) 
print(sam_record_obj.cigar)
print(sam_record_obj.left_most) 
print(sam_record_obj.dic_key)
print(sam_record_obj.qual)
print(sam_record_obj.rev_comp)
print(sam_record_obj.position_adj)
'''


