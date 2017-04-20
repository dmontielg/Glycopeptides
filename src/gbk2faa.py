#!/usr/bin/env python
#-- Copyright (c) 2015 Xiaowen Lu
#-- Bioinformatics Group @ Wageningen University  --
#--
#--Script that extract CDS from gbk file
#Libraries to import
import sys
import os
import string
from Bio import SeqIO, SeqFeature

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import warnings
with warnings.catch_warnings(): # Don't display the SearchIO experimental warning.
    warnings.simplefilter("ignore")
    from Bio import SearchIO

def get_CDS_fasta(gbk_file_dir, faa_filename):

    output_handel =  open(faa_filename, "w")
    gbk_filelist = [f for f in os.listdir(gbk_file_dir) \
    if os.path.isfile('/'.join((gbk_file_dir, f))) and not f.startswith('.')]

    for gbk in gbk_filelist:
        gbk_filename = os.path.join(gbk_file_dir, gbk)
        record = SeqIO.parse(gbk_filename, "genbank")

        for seq_record in record:
            print "Dealing with gbk record %s" % seq_record.id
            
            for seq_feature in seq_record.features:
                #screen the feather for each record, and get the CDS
                antibiotic = seq_record.description.split()[0]
                #print seq_feature
                
                if seq_feature.type == "CDS":
                    
                    if "translation" in seq_feature.qualifiers:
                    
                        assert len(seq_feature.qualifiers['translation']) == 1
                        
                        seq_qualifiers = seq_feature.qualifiers.keys()  
                        
                        if "protein_id" in seq_qualifiers:
                            protein_id = seq_feature.qualifiers['protein_id'][0]
                        else:
                            protein_id = seq_feature.qualifiers['locus_tag'][0]

                        if "gene" in seq_qualifiers:
                            gene = "|"+seq_feature.qualifiers['gene'][0]+"_"
                        else:
                            #gene = "?"
                            gene = "|_"
                            #gene = seq_feature.qualifiers['product'][0].replace(" ","_")      
                            
                            #gbk.split(".")[0],
                        output_handel.write(">%s|%s|%s%s\n%s\n" % (                                                        
                                seq_record.name,                                                        
                                antibiotic,                            
                                #seq_record.id,                            
                                protein_id,
                                
                                gene,
                                seq_feature.qualifiers['translation'][0]))
                    

    output_handel.close()
    
# gbk_file_dir = "/glycopeptides/gene_clusters/"
# fasta_filename = "/glycopeptides/glycopeptides.fasta"
# get_CDS_fasta(gbk_file_dir, fasta_filename)


if __name__ == "__main__":
    "Run this file from here as a script"
    #Check if parameters are provided; if not, exit with explanation
    if len(sys.argv) <= 2:
        print "Please provide as parameters a folder with GENBANK files and\
        a FASTA file name an output "
        sys.exit()

    #Read command-line parameters
    gbk_file_dir = sys.argv[1]
    fasta_filename = sys.argv[2]
    #Run get_CDS_fasta
    get_CDS_fasta(gbk_file_dir, fasta_filename)    
    
    
    
    
    
    
    
    