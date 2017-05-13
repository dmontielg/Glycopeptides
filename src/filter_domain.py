#!/usr/bin/env python
#Libraries to import

import sys

#Get genes found in two different fasta files

# Functions here
def parse_fasta(fastafile):
    "Get sequences from fasta file and store in dictionary"

    infile = open(fastafile).read()
    entries = infile.split(">")[1:]
    #print entries
    fastadict = {} #Accession as key, sequence as value
    for entry in entries:
        accession = entry.partition("\n")[0].replace("\r","")
        #print accession
        sequence = entry.partition("\n")[2].partition("\n")[0]
        #print sequence
        fastadict[accession] = sequence  

    return fastadict

def get_gene_intersection(a_domain,t_domain, output_file):
    
    output_handel =  open(output_file, "w")

    for a_gene in a_domain.items():  
        a_protein_id = a_gene[0].split("|")  
        a_protein_name = a_gene[0]
        a_protein_seq = a_gene[1]
        for t_gene in t_domain.items():
            t_protein_id = t_gene[0].split("|")            
            if a_protein_id[2] in t_protein_id[2] \
            and a_protein_id[3][-1] in t_protein_id[3][-1]:
                #print a_protein_name
                #print a_protein_seq
                
                output_handel.write(">%s\n%s\n" % (
                            a_protein_name,
                            a_protein_seq))
                
                            
def get_gene_diff(a_domain,t_domain, output_file):
    
    output_handel =  open(output_file, "w")
    for a_gene in a_domain.items():  
        a_protein_id = a_gene[0].split("|")  
        a_protein_name = a_gene[0]
        a_protein_seq = a_gene[1]        
        counter = False
        for t_gene in t_domain.items():
            t_protein_id = t_gene[0].split("|")                 
            if a_protein_id[2] in t_protein_id[2] \
            and a_protein_id[3][-1] in t_protein_id[3][-1]:                        
                counter = True                
        if counter is False:
            output_handel.write(">%s\n%s\n" % (
                            a_protein_name,
                            a_protein_seq))        
                                
if __name__ == "__main__":
    "Run this file from here as a script"
    #Check if parameters are provided; if not, exit with explanation
    if len(sys.argv) <= 3:
        print """Please provide as parameters a two domains sequences and an output file"""
        sys.exit()
    #Read command-line parameters
    a_domain_fasta = sys.argv[1]
    t_domain_fasta = sys.argv[2]
    output_file1 = sys.argv[3]
    #output_file2 = sys.argv[4]

    a_domain = parse_fasta(a_domain_fasta)
    t_domain = parse_fasta(t_domain_fasta)
    
    get_gene_intersection(a_domain,t_domain, output_file1)
    #get_gene_diff(a_domain,t_domain, output_file2)
    
    
    
    
    
    
    
    
    
    