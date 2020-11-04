#!/usr/bin/env python 

import statistics

# sam file record class

class SamRecord:

    def __init__(self, sam_line):
        
        self.line = sam_line
        self.line_list = sam_line.split()
        self.umi = self.get_umi()
        self.chr = self.line_list[2]
        self.bit_flag = self.line_list[4]
        self.cigar = self.line_list[5]
        self.left_most = self.line_list[3]
        self.dic_key = self.umi + "_" + self.chr
        self.qual = self.get_qual()
        self.rev_comp = (self.bit_flag & 16) == 16
        self.position_adj = self.get_position()
        
    # need to add adjusted starting position for actual alorithm
    #self.position = self.get_position() 

    # pulls UMI off of header
    def get_umi(self):     

        header = self.line_list[0]
        umi = header.split(":")[-1]

        return umi
    
    def get_qual(self):

        qual = self.line_list[10]
        qual_list = [ord(asci) - 33 for asci in qual]

        return statistics.mean(qual_list)

    def get_position(self):

        return None
    # will be function to get position
    # def get_position(self):


# testing 
'''
sam_record = "NS500451:154:HWKTMBGXX:1:11101:21621:1145:AGGTTGCT	0	2	93022350	36	71M	*	0	0	TTCCACTGTTGCTTCATAACTGCAGTCCTAACATAAATGTCTGACATGTAGGATGATCTTAAGCAACCCCT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEAAAEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU"

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
'''