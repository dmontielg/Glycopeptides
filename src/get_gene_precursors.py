#!/usr/bin/env python
"""
Author: Diego Fernando Montiel Gonzalez
File number: 880505580110
Bioinformatics Group, Wageningen University
Get gene precursors. Retrieve Genbank NCBI files 
from  from an external ACCESSION number list
and generate a multi-fasta.
"""

from Bio import Entrez
from sys import argv
import os.path
import re

def processOrigin(origin_line):
    """
    -Function: processOrigin
    Description: 
        Function called at the parseGenbankFile for parsing the origin part of the
        file, getting only the nucleotides
    Input: 
        fragment of the genbankfile where the origin is located
    Output: 
           Returns a complete string of all the nucleotides
    """    
    #get the nucleotides of each line except numbers
    line = re.findall(r"[^0-9]+", origin_line, re.I)
    output = ''
    for x in line:
        output = output + x.replace(" ","")
    return output

def parseGenBankFile(filename):
    """
    -Function: parseGenBankFile
    Description: 
        Function that parses a genbankfile separating the alignments and 
        sequences intor records
        
    Input: 
        A complete genbankfile
    Output: 
        accession,organism & origin
    """    
    nucleotides = ''    
    accession = ''
    organism = ''
    origin = ''
    #Use as a condition to stio the searching and if it finds the // string
    #means that has get to the end of that sequence and now all the string is
    #append into the list origin
    flag        = 0
    for line in filename:
        if re.match(r'ACCESSION',line):
            tmp = line.replace("ACCESSION","")
            accession = tmp.strip()
            
        elif re.match(r'  ORGANISM',line):
            tmp = line.replace("ORGANISM","")
            organism = tmp.replace("\n","")
                        
        elif re.match(r'ORIGIN',line):                                                    
            flag = 1
            
        elif re.match(r'//',line):
            origin = processOrigin(nucleotides)
            nucleotides = ''            
            flag = 0             
                                        
        elif flag == 1:         	
            nucleotides = nucleotides + line.strip()
            
    return accession, organism, origin

def generateFastaFile(accession, organism, sequence,output_file):

    with open(output_file,"a+") as file_fasta:
        file_fasta.write('>'+accession)      
        file_fasta.write('\n')
        file_fasta.write(sequence.upper())
        file_fasta.write('\n')
        print 'Registry added to the multi-fasta!'
        return True
    
def readFileInput(list_genes):
    """
    Function: readFileInput
    Description: Get all the sequences of each identifier in the inputfile
    as a fasta format
    Input: Text file with the identifiers to look for: Ex. "NM_1799883"
    Output: List with all the sequences (similar to a fasta file)
    """
    Entrez.email = 'diego.montielgonzalez@wur.nl'
    #Read a file line by line    
    genbank_files = []
    for iterator in list_genes:
        #Remove all spaces
        iterator = iterator.strip()
        handle = Entrez.efetch(db="protein",id=[iterator], rettype="gb")
        gb_records = handle.read()
        iterator = iterator+".gb"
        generateGenBankFile(iterator,gb_records)
        genbank_files.append(iterator)
    
    return genbank_files

def generateGenBankFile(genbank_file_name,gb_records):
    """
    Function: generateGenBankFile
    to generate the genbankfile and then use the function of the p2 script
    to parse the genbank file
    Input: String with the Genbank records
    Output: Genbank file or an False boolean statement
    """
    if not os.path.isfile(genbank_file_name):
        with open(genbank_file_name,"a+") as gb_file:
            gb_file.write(gb_records)
            if os.path.isfile(genbank_file_name):
                print "Succes!, Genbank File succesfully generated!"
                return genbank_file_name
            else:
                print "Failed to generated Genbank File" 
                return False
        
if __name__ == '__main__':
    """
    Inputs: Script and text file with genbank identifiers
        Ex. python p3.py p3input.txt
    """
    file_genes      = argv[1]
    output_file     = argv[2]    
    list_genes      = [line.rstrip('\n') for line in open(file_genes)]
    genbank_files   = readFileInput(list_genes)
    
    for gb in genbank_files:
        file_fasta = open(gb)
        accession, organism, sequence = parseGenBankFile(file_fasta)      
        generateFastaFile(accession, organism, sequence,output_file)
       
       
    
    
    
    
    