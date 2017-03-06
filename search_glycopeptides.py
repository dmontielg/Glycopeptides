
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

def parse_hmmersearch(hmm_output):

    dictionary = {} #Accession as key, sequence as value
    for entry in hmm_output:
        if "#" not in entry:
            item = []            
            line = entry.split()
            
            tmp_target = line[0].split("|")            
            target_name = line[0]
            
            domain = line[3]
            ievalue = line[12]
            score = float(line[13])      
            locus = tmp_target[-1]
            
            item.append(locus)
            item.append(domain)
            #item.append(ievalue)
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
    contador = 0
    #Order the dictionary by organism and cluster
    dictionary = collections.OrderedDict(sorted(dictionary.items()))
    for key in dictionary:          
        #Each gene cluster
        counts = {'AMP-binding':0,'Epimerization':0,'Cglyc':0,'X':0}
        ##loop all the items in the cluster/key
        for item in dictionary[key]:            
            if item[1] in counts:
                counts[item[1]] += 1      
        #Condition for each gene cluster
        if counts["AMP-binding"] > 4 and counts["Epimerization"]>=1\
            and counts["X"] >= 1 and counts["Cglyc"] >= 1:
            
            ############################################
            ## #print key Still the cluster...
            ## Starts from here
            ############################################                        
            ##Loop again the whole dictionary that fullfill the criteria
            # Not used
            c_domain = {'Condensation_DCL':0,'Condensation_LCL':0,'Condensation_Dual':0,'Condensation_Starter':0}            
            #print key
            x = []
            c = []
            g = []
            e = []
            #print key
            for array_items in dictionary[key]:
                #Contains the protein/domain/score
                if "X" in array_items:
                    #print array_items                    
                    x.append(array_items)
                if "Condensation" in array_items[1]:                    
                    #print array_items                
                    c.append(array_items[2])
                if "Cglyc" in array_items:
                    #print array_items
                    g.append(array_items)                    
                if "Epimerization" in array_items:
                    #print array_items
                    e.append(array_items)
                    
            ######################
            ## X domain rule 
            #####################
            ax = True
            for i in range(len(x)):
                value_x =  x[i][2]   
                if ax:
                    for j in range(len(c)):
                        value_c = c[j]
                        #if at least one x in c is minor breaks, rule does not work
                        if value_x < value_c:
                            ax = False
                            #print value_x
                            #print value_c                                                        
            if ax:
                print key                            
            ######################
            ## Cglyc domain rule 
            #####################
            ag = True
            for i in range(len(g)):
                value_g =  g[i][2]   
                if ag:
                    for j in range(len(c)):
                        value_c = c[j]
                        #if at least one x in c is minor breaks, rule does not work
                        if value_g < value_c:
                            ag = False
                            #print value_x
                            #print value_c                                                        
            if ag:
                print key
            ######################
            ## Epimerization domain rule 
            #####################            
            count = len(e)                       
            for i in range(len(e)):
                value_e = e[i][2]
                for j in range(len(c)):
                    value_c = c[j]
                    if value_e < value_c:                          
                            count -= 1            
                            break         
            if count > 0:
                print key
                contador += 1                
    print contador
    
    #i and x are the most important variables so far...
    ###To do:
    ## Use the combinations of the new key to somewhat compare scores for the domains
    ## and see if they can be considered as the ones responsible for glycopeptides
    ## also check the genes number of genes or so...
    
    return True


def remove_duplicates(l):
    return list(set(l))

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
        #Read command-line parameters
        hmm_output = sys.argv[1]
        outfile = sys.argv[2]        
        lines = [line.rstrip('\n') for line in open(hmm_output)]
        dictionary = parse_hmmersearch(lines)        
        new_dictionary = find_products(dictionary)        
        #np.save(outfile, new_dictionary) 
        # Load
        #read_dictionary = np.load(outfile).item()        
        #print len(read_dictionary.keys())
        print("--- %s seconds in parse hmmer file ---" % (time.time() - start_time))
        
        
        
