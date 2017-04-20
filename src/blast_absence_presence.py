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