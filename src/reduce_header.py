#!/usr/bin/env python

#-- Copyright (c) 2017 Diego Montiel
#-- Bioinformatics Group @ Wageningen University  --
#--
#-- Script for reducing header in a fasta file 


import sys

import re

def parseFastaFile(filename):
    
    #stores temporally and sum up each line of one sequence
    nucleotides = ''
    #validate the parse to read the label or the sequence
    header = []
    sequence = []
    flag = False
    for line in open(filename):
        if re.match(r'>',line):                      
            key_header = line.rstrip()
            header.append(line.rstrip())
            if flag:
                sequence.append(nucleotides)
                nucleotides = ''              
                flag = False
        else:
           flag = True
           nucleotides = nucleotides + line
    sequence.append(nucleotides)    
    return header, sequence

def change_header(header,sequence, output_file):
    
    output_handel =  open(output_file, "w")
    c = 0
    
    for index in range(len(header)):
            
        protein_id = header[index].split("|")          
	
	##protein_header = c
        ##c+=1
        ##protein_header += protein_id[2]
        ##protein_header += protein_id[2]+"_"+protein_id[5]

	protein_header = protein_id[0]+"|"+protein_id[1]+"|"+protein_id[5]+"|"+protein_id[6]
        #print protein_header
	
	#protein_header = protein_id[0]+"|"+protein_id[1]+"|"+protein_id[5]+"|"+protein_id[6]+"|"+protein_id[7]+"|"+protein_id[8]
        
	protein_seq = sequence[index]
        ##output_handel.write(">%s\n%s\n" % (protein_header,protein_seq))
        output_handel.write("%s\n%s" % (protein_header,protein_seq))                            
                                
if __name__ == "__main__":
    "Run this file from here as a script"
    #Check if parameters are provided; if not, exit with explanation
    if len(sys.argv) < 2:
        print """Please provide as parameters FASTA file and an output file"""
        sys.exit()
    
    fastafile = sys.argv[1]
    output_file = sys.argv[2]
    header, sequence = parseFastaFile(fastafile)
    #print a_domain
    change_header(header, sequence, output_file)
  
    
    
    
    
    
    
    
    
    
    