#!/usr/bin/env python 

# linux command functions

'''
this is a customized package consisting of python wrappers of common 
linux/unix shell commands. enjoy. 
'''

import subprocess 
import os 
  
# a function to list the files in 
# the current directory and  
# parse the output. 
def list_command(args = '-l'): 
  
    # the ls command 
    cmd = 'ls'
  
    # using the Popen function to execute the 
    # command and store the result in temp. 
    # it returns a tuple that contains the  
    # data and the error if any. 
    temp = subprocess.Popen([cmd, args], stdout = subprocess.PIPE) 
      
    # we use the communicate function 
    # to fetch the output 
    output = str(temp.communicate()) 
      
    # splitting the output so that 
    # we can parse them line by line 
    output = output.split("\n") 
      
    output = output[0].split('\\') 
  
    # a variable to store the output 
    res = [] 
  
    # iterate through the output 
    # line by line 
    for line in output: 
        res.append(line) 
  
    # print the output 
    #for i in range(1, len(res) - 1): 
    #    print(res[i]) 
  
    return res 

def line_count(file_name, args = "-l"): 
    
    # the ls command 
    cmd = 'wc'
  
    # using the Popen function to execute the 
    # command and store the result in temp. 
    # it returns a tuple that contains the  
    # data and the error if any. 
    command = cmd + " " + args + " " + file_name 
    temp = subprocess.Popen(command, stdout = subprocess.PIPE) 
    print("temp")
    print(temp)
    # we use the communicate function 
    # to fetch the output 
    output = str(temp.communicate()) 
    print("output")
    print(output)
    # splitting the output so that 
    # we can parse them line by line 
    output = output.split("\n") 
    print("splitout")
    print(output)
      
    #output = output[0].split('\\') 
  
    # a variable to store the output 
    #res = [] 
  
    # iterate through the output 
    # line by line 
    #for line in output: 
    #    res.append(line) 
  
    # print the output 
    #for i in range(1, len(res) - 1): 
    #    print(res[i]) 
  
    return res 


x = line_count("filtered_big.sam", args = '-l')
print(x)
