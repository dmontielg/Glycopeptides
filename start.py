#!/usr/bin/env python

## Copyright (c) 2017 Diego Montiel
## Wageningen University
## Bioinformatics Group

###Script for prediction of Glycopeptide Gene Clusters 

from sys import path
path.append("../src")

import os
import os.path
import sys
import subprocess
import time
import search_glycopeptides as sg
import gbk2faa
import extract_domain
import filter_domain

def menu():

    ans=True
    while ans:
        print ("""
        #### Wageningen UR ####
        #### Bioinformatics - Marnix Medema group ####
        Prediction of putative Glycopeptide Biosynthetic Gene Cluster
        1. Retrieve CDS from genbank
        2. Predict putative GBGC
        3. Extract AT domain from sequences
        4. Perform MSA with MUSCLE
        5. Maximum-Likelihood with RaxML
        6. Exit/Quit
        """)
        ans=raw_input("What would you like to do? ") 
        if ans=="1":             
            print("\n Retrieve CDS from a genbank folder")           
            gbk_file_dir = raw_input("(Ex. data/novel_gene_clusters)\n Please provide Genbank folder path: ")            
            output_fastafile = raw_input("\n Please provide output name for fasta file from Genbanks: ")
            #Run get_CDS_fasta
            gbk2faa.get_CDS_fasta(gbk_file_dir, output_fastafile)                  
            if os.path.exists(output_fastafile):
                print "Success! \nFile "+output_fastafile+" generated"
            else: 
                print "Error! \nSomething went wrong, no file was generated!"
                
        elif ans=="2":
            print("\n Predict putative GBGC")                         
                        
            hmmversion = "v3"
            hmmfile = "../data/hmm/nrpspksdomains.hmm"
            cut_off = "0.01"
            
            fastafile = raw_input("\n Please provide sequence query for domain profile: ")
            hmm_output = raw_input("\n Please provide hmm output name: ")            
            fastadict = extract_domain.read_fasta_file(fastafile)            
            hmmer_results = extract_domain.run_HMMer(hmmversion, hmmfile, fastafile, cut_off,hmm_output)
            if os.path.exists(hmm_output):
                print "Success! \nFile "+hmm_output+" generated"
            else: 
                print "Error! \nSomething went wrong, no file was generated!"        
                    
            print "#################################"
            print "### Looking for putative GBGC..."
            print "#################################"
            start_time = time.time() #Global variable
                        
            #hmm_output = "hmmsearch3_cutoff.out"
            #fastafile = "reduce_gene_clusters_nrps.fasta"
            
            lines = [line.rstrip('\n') for line in open(hmm_output)]
            dictionary = sg.parse_hmmersearch(lines)        
            putative_glycopeptides = sg.find_products(dictionary)        
            print("--- %s seconds in parse hmmer file ---" % (time.time() - start_time))
            if len(putative_glycopeptides) > 0:
                #output for fasta file
                print "\nSuccess! We found putative Glycopeptide Gene clusters!"                
                out_fasta_file = raw_input("\n Please provide a name for this fasta file: ")
                sg.extract_sequence(putative_glycopeptides,fastafile,out_fasta_file)      
            else:
                print "\n Failed! to found putative Glycopeptide Gene Clusters!"
        
        elif ans=="3":
            print("\n Filter sequences with AT domains")                         

            hmmversion = "v2"
            cut_off = "0.01"
            
            hmmfile = "../data/hmm/A_domain_232.hmm"
            domain_ab = "A"
            flag = False
            
            fastafile = raw_input("\n Please provide sequence query for domain profile: ")
            #at_output = raw_input("\n Please provide the output name for AT domains: ")    
	    at_output = "AT_domain_hmm2_"+fastafile 
            a_output = "output_A_domain_hmm2_"+fastafile
            a_output_fastafile = "A_domain_hmm2_"+fastafile            
            
            fastadict = extract_domain.read_fasta_file(fastafile)            
            hmmer_results = extract_domain.run_HMMer(hmmversion, hmmfile, fastafile, cut_off,a_output)
            extract_domain.write_fasta_hmm2(fastadict, a_output_fastafile, hmmer_results, domain_ab)
            if os.path.exists(a_output_fastafile):
                print "Success! \nFile "+a_output_fastafile+" generated"
                flag = True
            else: 
                print "Error! \nSomething went wrong, no file was generated!"    
                flag = False
            
            hmmfile = "../data/hmm/T_domain_232.hmm"
            domain_ab = "T"
            flag = False
            
            t_output = "output_T_domain_hmm2_"+fastafile            
            t_output_fastafile = "T_domain_hmm2_"+fastafile            
            
            fastadict = extract_domain.read_fasta_file(fastafile)            
            hmmer_results = extract_domain.run_HMMer(hmmversion, hmmfile, fastafile, cut_off,t_output)
            extract_domain.write_fasta_hmm2(fastadict, t_output_fastafile, hmmer_results, domain_ab)
            if os.path.exists(t_output_fastafile):
                print "Success! \nFile "+t_output_fastafile+" generated"
                flag = True
            else: 
                print "Error! \nSomething went wrong, no file was generated!"    
                flag = False            
            if flag:
                a_domain = filter_domain.parse_fasta(a_output_fastafile)
                t_domain = filter_domain.parse_fasta(t_output_fastafile)                
                filter_domain.get_gene_intersection(a_domain,t_domain, at_output)
                if os.path.exists(at_output):
                    print "Success! \nFile "+at_output+" generated"


                    filenames = [at_output, '../data/glycopeptides/Prediction_Glycopeptide_AT_domain.fasta']
		    #filenames = [at_output, '../data/glycopeptides/Glycopeptide_AT_domain.fasta']
                    #filenames = [at_output, '../data/glycopeptides/Glycopeptide_Novel_AT_domain.fasta']

                    with open('Glycopeptide_Novel_AT_domain.fasta', 'w') as outfile:
                        for fname in filenames:
                            with open(fname) as infile:
                                outfile.write(infile.read())
                
                print "Success! \nFile Glycopeptide_Novel_AT_domain.fasta generated"

            else:
                print "Failed! Something went wrong, please try again..." 
        
        elif ans=="4":
            
            print("\nPerform MSA with MUSCLE") 
            fastafile = raw_input("\n Please provide multi-fasta sequence for MSA: ")
            output_fastafile = raw_input("\n Please provide output for MSA: ")

            if subprocess.check_call('muscle -in ' +fastafile + ' -out ' +output_fastafile, shell=True) == 0:                        
                print 'Success! MSA performed'
                print "\nFile "+output_fastafile+" generated"
                if subprocess.check_call('perl ../src/Fasta2Phylip.pl '+output_fastafile + ' ' +output_fastafile+".phy", shell=True) == 0:       
                    print "\nFile "+output_fastafile+".phy generated"
            else: 
                print 'Something went wrong, check your multi-fasta file..'
        
        elif ans=="5":  
            print "Phylogeny analysis with Maximum-Likelihood using RaxML"
            
            align_fasta = raw_input("\n Please provide alignment multi-fasta: ")
            cpu = '4'
            bb = '100'
            if subprocess.check_call('../external/standard-RAxML-master/raxmlHPC-PTHREADS-AVX -T '+cpu+' -f a -m PROTGAMMAIWAGF -s '+align_fasta+' -p 15000 -x 16000 -#'+bb+' -n RaxML_'+align_fasta, shell=True) == 0:      
                 print 'Success! RaxML performed'
               
        elif ans=="6":
            print("\nGoodbye") 
            ans = None            
        elif ans !="":
            print("\nWarning!: Not Valid Choice Try again")
        
if __name__ == "__main__":

    menu()
    
    
    
    
    