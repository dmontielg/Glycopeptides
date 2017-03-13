
#!/usr/bin/env python
"""
Bioinformatics Group Wageningen UR
Author: Diego Montiel 
"""
import sys
import time
import numpy as np
import collections
import operator
import itertools
import logging
import os
import shutil
from os import path
from operator import itemgetter
import subprocess
import string

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
from argparse import Namespace

import warnings
# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO

from Bio.Alphabet import generic_protein
from zipfile import ZipFile, ZIP_DEFLATED, LargeZipFile
from Bio.Seq import Seq
import re

def parse_hmmersearch(hmm_output):
    dictionary = {} #Accession as key, sequence as value
    for entry in hmm_output:
        if "#" not in entry:
            item = []            
            line = entry.split()            
            tmp_target = line[0].split("|")            
            target_name = line[0]            
            domain = line[3]
            ievalue = float(line[12])
            score = float(line[13])      
            cord_start = int(line[19])      
            cord_end = int(line[20])#evelope
            locus = tmp_target[-1]            
            item.append(locus)            
            item.append(domain)
            item.append(cord_start)
            item.append(cord_end)
            item.append(ievalue)
            item.append(score)

            #Use locus_Tag+cluster+protein_id
            key_header = (tmp_target[0],tmp_target[1])                                    
            if key_header in dictionary:
                dictionary[key_header].append(item)           
            else:
                dictionary[key_header] = [item]     
    return dictionary

def find_products(dictionary):
    # filter by at least 4 A domains
    #AEK43341
    putative_glycopeptides = []
    #Order the dictionary by organism and cluster
    dictionary = collections.OrderedDict(sorted(dictionary.items()))
    for key in dictionary:                 
        #Each gene cluster
        contador = 0
        found = 0
        c_domain = []
        counts = {'AMP-binding':0,'Epimerization':0,'Cglyc':0,'X':0}
        ##loop all the items in the cluster/key
        for item in dictionary[key]:   
            if item[1] in counts:
                counts[item[1]] += 1      
        #Condition for each gene cluster
        if counts["AMP-binding"] > 4 and counts["Epimerization"]>=1\
            and counts["X"] >= 1 and counts["Cglyc"] >= 1:
 
            domainlist = []
            dictionary[key].sort(key=lambda x:(x[0],x[1],x[3],x[4]))
            for array_items in dictionary[key]:
                domainlist.append(array_items)
                if "Condensation" in array_items[1]:                    
                    c_domain.append(array_items)
            domain = check_overlap(domainlist)
            #print key
            #for i in domain:
                ##See how the gene cluster assembly line purify
                #print i
            #check overlaps between condensation domains
            domain1 = check_domain_nrps(domain,c_domain)
            e = "Epimerization" 
            found += found_in_list(domain1,e)
            x = "X"
            found += found_in_list(domain1,x)
            g = "Cglyc"
            found += found_in_list(domain1,g)
            if found ==3:
                #print key
                putative_glycopeptides.append(key)
          
    # Check for the 4 A domains
    return putative_glycopeptides

def found_in_list(domain1,search):
    contador = 0
    for sublist in domain1:
        if sublist[1] == search:
            #print the position 1 in the list = domain
            #print sublist[1]
            contador+=1
            break
    return contador

def check_overlap(pfd_matrix):
    """Check if domains overlap for a certain overlap_cutoff.
     If so, remove the domain(s) with the lower score."""
    
    overlap_cutoff = 0.1
    delete_list = []
    for i in range(len(pfd_matrix)-1):
        for j in range(i+1, len(pfd_matrix)):
            row1 = pfd_matrix[i]
            row2 = pfd_matrix[j]
            
            #check if we are the same CDS
            if row1[0] == row2[0]:
                #check if there is overlap between the domains
                if no_overlap(int(row1[2]), int(row1[3]), int(row2[2]), int(row2[3])) == False:
                    overlapping_nucleotides = overlap(int(row1[2]), int(row1[3]), int(row2[2]), int(row2[3]))
                    overlap_perc_loc1 = overlap_perc(overlapping_nucleotides, int(row1[3])-int(row1[2]))
                    overlap_perc_loc2 = overlap_perc(overlapping_nucleotides, int(row2[3])-int(row2[2]))
                    #check if the amount of overlap is significant
                    if overlap_perc_loc1 > overlap_cutoff or overlap_perc_loc2 > overlap_cutoff:
                        if float(row1[5]) > float(row2[5]): #see which has a better score
                            delete_list.append(row2)
                        elif float(row1[5]) < float(row2[5]):
                            delete_list.append(row1)
    for lst in delete_list:
        try:
            pfd_matrix.remove(lst)
        except ValueError:
            pass
        
 
    return pfd_matrix  

def no_overlap(locA1, locA2, locB1, locB2):    
    """Return True if there is no overlap between two regions"""
    if locA1 < locB1 and locA2 < locB1:
        return True
    elif locA1 > locB2 and locA2 > locB2:
        return True
    else:
        return False
    
def overlap_perc(overlap, len_seq):
    return float(overlap) / len_seq
    
def overlap(locA1, locA2, locB1, locB2):
    """Returns the amount of overlapping nucleotides"""
    if locA1 < locB1:
        cor1 = locA1
    else:
        cor1 = locB1
    if locA2 > locB2:
        cor2 = locA2
    else:
        cor2 = locB2
    total_region = cor2 - cor1
    sum_len = (locA2 - locA1) + (locB2 - locB1)
    return sum_len - total_region

def check_domain_nrps(domain_a, domain_b):
    domain_list = []
    #lista = [(3015, 3701), (4011, 5890)]
    #listb = [(1,2), (100,200), (4500,6000)]
    #result = check_overlap_coordinates(lista, listb)
    #print result

    for a in domain_a:
        lista = [tuple([a[2],a[3]])]
        for b in domain_b:
            listb = [tuple([b[2],b[3]])]
            result = check_overlap_coordinates(lista, listb)
            #If there is overlapping
            if result is not False:
                #If first coordinate has greater score and are in same protein/gene
                if a[-1] > b[-1] and a[0] == b[0]:                    
                    domain_list.append(a)
                    #print a
                    #break
                #elif a[-1] < b[-1] and a[0] == b[0]:
                    #domain_list.append(b)

          
    return domain_list

def check_overlap_coordinates(lista, listb):
    a = 0
    b = 0
    found = False
    while a < len(lista) and b < len(listb):
        result = check( lista[a] , listb[b] )
        if result < 0:
            a += 1
            continue
        if result > 0:
            b += 1
            continue
        # we found overlapping intervals
        found = True
        return (found, a, lista[a], b, listb[b] )
    return found

def check( (astart, aend) , (bstart, bend) ):
    if aend < bstart:
        return -1
    if bend < astart:
        return 1
    return 0

def remove_duplicates(l):
    return list(set(l))

def extract_sequence(putative_glycopeptides, fasta_file):    
    header = []
    for item in putative_glycopeptides:
        header = ">"+item[0]+"|"+item[1]
        #header.append(tmp)    
        with open(fasta_file, 'r+') as f:
            while True:
                line = f.readline().rstrip()
                #print line            
                #print header
                if header in line:                    
                    seq = f.readline().rstrip()                    
                    generateFastaFile(line,seq)
                if len(line) == 0:
                    break    
            
def generateFastaFile(header,sequence):
    with open("Novel_BGC_glycopeptides2.fasta","a+") as file_fasta:
        file_fasta.write(header)      
        file_fasta.write('\n')
        file_fasta.write(sequence.upper())
        file_fasta.write('\n')
        #print 'Registry added to the multi-fasta!'
        return True

if __name__ == "__main__":
    #python search_glycopeptides.py Glycopeptides/temp.out x    
    "Run this file from here as a script"
    #Check if parameters are provided; if not, exit with explanation
    start_time = time.time() #Global variable
    if len(sys.argv) < 2:
        print """Please provide as parameters an output HMM file, and \
        an output file"""
        sys.exit()
    else:
        #python search_glycopeptides.py data/nrpspksdomains.hmm nrps/hmmsearch3_cutoff.out data/reduce_gene_clusters_nrps.fasta
        #Read command-line parameters
        hmm_file = sys.argv[1]
        hmm_output = sys.argv[2]        
        fasta_file = sys.argv[3]        
        lines = [line.rstrip('\n') for line in open(hmm_output)]
        dictionary = parse_hmmersearch(lines)        
        putative_glycopeptides = find_products(dictionary)        
        #print len(putative_glycopeptides)
        extract_sequence(putative_glycopeptides,fasta_file)
        
        print("--- %s seconds in parse hmmer file ---" % (time.time() - start_time))
        
        
        
