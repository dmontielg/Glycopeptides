#!/usr/bin/env python
"""
Author: Diego Fernando Montiel Gonzalez
File number: 880505580110
Bioinformatics Group, Wageningen University
BLAST the absence presence from gene precursors
"""

from sys import argv
import get_gene_precursors as my_functions
import subprocess
import os.path
import re

def runningNeedle(related_fasta, ref_fasta, gapopen=8):

    #Check if the output file already exists, if not, execute needle    
    if not os.path.isfile('out.needle'):
        #print(subprocess.check_call('needle -asequence ' + ref_fasta + ' -bsequence ' + related_fasta + ' -gapopen ' , int(gapopen), ' -gapextend 0.5 -outfile out.needle',shell=True) == 0)
       
        if subprocess.check_call('needle -asequence ' +ref_fasta + ' -bsequence ' +related_fasta+ ' -gapopen 8 -gapextend 0.5 -outfile out.needle',shell=True) == 0:
                        
            return 'Success! align sequence of protein with Needle..'
        else: 
            return 'Something went wrong, file missing..'

def getOriginGenbank(gbk_file_dir,output_file):
    
    gbk_filelist = [f for f in os.listdir(gbk_file_dir) \
        if os.path.isfile('/'.join((gbk_file_dir, f))) and not f.startswith('.')]

    for gbk in gbk_filelist:
        gbk_filename = os.path.join(gbk_file_dir, gbk)
        gbk_file = open(gbk_filename)
        accession, organism, sequence = my_functions.parseGenBankFile(gbk_file)
        my_functions.generateFastaFile(accession, organism, sequence,output_file)

if __name__ == '__main__':
    """
    Inputs:(4 arguments, script, related fasta file, reference fasta file and gapopen number
    ex. python p5.py related.fasta ref.fasta 8
    """     
    
    folder_genbank   = argv[1]
    output_file      = argv[2]
    getOriginGenbank(folder_genbank,output_file)
    
    
    """
    related_fasta   = argv[1]    
    ref_fasta       = argv[2]    
    gapopen         = argv[3]    
    runningNeedle(related_fasta, ref_fasta, gapopen)    
    ref_parsed_fasta     = parseFastaFile(ref_fasta)
    related_parsed_fasta = parseFastaFile(related_fasta)
    if os.path.isfile('out.needle'):
        parsed_needle_dictionary = parseNeedle('out.needle')
        showOutput(parsed_needle_dictionary)       
        
    """