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
                            
def add_substrate(a_domain,substrate, output_file):
        
    output_handel =  open(output_file, "w")
    print "Start checking headers..."
    
    for header,sequence in a_domain.items():
        h_a_domain = header.split("|")
        for sub in substrate.keys():
            s_substrate = sub.split("|")
            #print h_a_domain[2]+h_a_domain[3]
            #print s_substrate[2]+h_a_domain[3]
            #print header+"|"+s_substrate[4]+"|"+s_substrate[5]
            if h_a_domain[2]+h_a_domain[3] == s_substrate[2]+s_substrate[3]:
                new_header = header+"|"+s_substrate[4]+"|"+s_substrate[5]
                #print h_a_domain[2]+h_a_domain[3]
                #print s_substrate[2]+s_substrate[3]
                
                output_handel.write(">%s\n%s\n" % (
                            new_header,
                            sequence))
      
    print "Finish see output generated."
    return True
                                
if __name__ == "__main__":
    "Run this file from here as a script"
    #Check if parameters are provided; if not, exit with explanation
    
    if len(sys.argv) < 3:
        print "Please provide as parameters a HMM file, a FASTA file and an output file"
        sys.exit()
    
    #Read command-line parameters
    a_domain_fasta = sys.argv[1]
    a_domain_substrate = sys.argv[2]
    #output = "list_BGC_.txt"
    output = sys.argv[3]
    
    a_domain = parse_fasta(a_domain_fasta)
    substrate = parse_fasta(a_domain_substrate)
    
    add_substrate(a_domain, substrate, output)
    
    
    
    
    
    
    
    
    
    
    
