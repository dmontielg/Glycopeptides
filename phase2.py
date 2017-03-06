
#!/usr/bin/env python

"""
Bioinformatics Group Wageningen UR
Author: Diego Montiel 
"""
import time
import sys

def parse_fasta(fastafile):
    "Get sequences from fasta file and store in dictionary"
    infile = open(fastafile).read()
    entries = infile.split(">")[1:]
    fastadict = {} #Accession as key, sequence as value
    key_translation_dict = {}
    for entry in entries:
        header = entry.partition("\n")[0].replace("\r","")
        header_pieces = header.split("|")
        key_translation_dict[">"+header_pieces[0]+"-"+header_pieces[-1]] = header
        sequence = entry.partition("\n")[2].partition("\n")[0]
        fastadict[header] = sequence  

    return fastadict, key_translation_dict

def write_fasta(fasta_dict, key_translation_dict, index_list,output):
    #Writes a Fasta with the new key identified in the main dictionary
    output =  open(output, "w")

    for new_key in index_list:
        #print new_key
        if new_key in key_translation_dict:
            old_key = key_translation_dict[new_key]
            output.write(">"+old_key+"\n")
            output.write(fasta_dict[old_key]+"\n")
            
def parse_index(headers):
    #Parse the headers file in an index list
    index_list = []
    output = open("list_nrps_index.txt","w")
    for entry in headers:
        line = entry.split("\t")
        proteins = line[-1].split(";")
        if "nrps" == str(line[3]):                                                                                
            output.write(entry+"\n")
	    for prot in filter(None, proteins):         
                index_list.append(">" + line[0]+"-" +prot) 
             
    output.close()
    return index_list

if __name__ == "__main__":
    "Run this file from here as a script"
    #Check if parameters are provided; if not, exit with explanation
    
    #Example to run the script:
    
    #python phase2.py geneclusters.txt geneclusterprots.fasta
    
    start_time = time.time() #Global variable
    if len(sys.argv) < 1:
        print "Please provide three parameters Blast output file output \
        file name and cutoff"
        sys.exit()
    else:        
        gene_cluster_index = sys.argv[1]
        gene_cluster_fasta = sys.argv[2]
        output = sys.argv[3]
	lines = [line.rstrip('\n') for line in open(gene_cluster_index)]        
        index_list = parse_index(lines)
        fasta_dict, key_translation_dict = parse_fasta(gene_cluster_fasta)
        write_fasta(fasta_dict, key_translation_dict, index_list,output)
        
        print("--- %s seconds in parse index file ---" % (time.time() - start_time))
        
       
