#!/usr/bin/env python
#-- Copyright (c) 2016 Xiaowen Lu
#-- Bioinformatics Group @ Wageningen University  --
#-- modified by D. Montiel, J. Navarro (2017)
#-- This is the module to calculate the distance between each CLP pair

import sys
from munkres import Munkres
import math
import itertools
import pandas as pd
import numpy as np

def get_distance(DM, compd1, compd2, index1, index2):
    DM_input = open(DM, 'r').read().split('\n')
    index = [0,5]
    domain_name = []
    for l in DM_input:
        if l != "":
            line = l.strip()
            domain = line.split(",")[0].split("|")[0] + "|" + line.split(",")[0].split("|")[-1]
            #domain = "|".join([l.split(",")[0].split("|")[i] for i in index])
            domain_name.append(domain)

        
    #sys.exit()

    domain_index = dict(zip(domain_name, range(0,len(domain_name))))
    A = ['%s|%s' % t for t in zip([compd1]*len(index1), index1)]
    B = ['%s|%s' % t for t in zip([compd2]*len(index2), index2)]

    output = []

    for i in range(len(A)):
        A_i = A[i]
        B_i = B[i]

        row_index = domain_index.get(A_i)
        col_index = domain_index.get(B_i)

        if row_index <= col_index:
            row = DM_input[col_index]
            col = row_index
        else:
            row = DM_input[row_index]
            col = col_index
        score_list = row.split(",")
        #score_list.pop()
     
        #print score_list
        score = 1-float(score_list[col+1])
        
        output.append(score)
    
    return output

def generate_BGCs_DMS(seq_simScore, cluster, pseudo_aa):
    input = open(cluster, 'r').read().split('>')[1:]
    cluster_id = [i.split('\n')[0] for i in input]
    pseudo_seq = [i.split('\n')[1] for i in input]
    cluster_seq = dict(zip(cluster_id, pseudo_seq))

    BGCs = {}
    """
    BGC[compound_name][specificity_letter] = [list of labels for the specificity]
    e.g. BGC["Complestatin"]["H"] == ["Complestatin_0", "Complestatin_2", ...]
    """
    
    # for every compound (header in pseudo-fasta)
    for id in cluster_seq.keys():
        cluster_i = {} # define a dictionary for every specificity_letter

        seq = cluster_seq[id] # sequence of specificities
        seq_dict = {}
        for i in range(len(seq)):
            seq_dict[id+'_'+str(i+1)] = seq[i] # e.g. seq_dict["Complestatin_0"] = H

        specificity_name = set(seq_dict.values())
        for d in specificity_name: # for every DIFFERENT letter in the sequence...
            key = d
            value = []
            for index, specificity in seq_dict.items():
                if specificity == d:
                    value.append(index)
            cluster_i[key] = value
        BGCs[id] = cluster_i
        # in the end a made-up label is made for every specificity in each BGC 
        # (which also contains the position of the specificity)



    DMS = {}
    """DMS -- dictionary of this structure: DMS = {'general_domain_name_x': { ('specific_domain_name_1',
    'specific_domain_name_2'): (sequence_identity, alignment_length), ... }   }
        - general_domain_name_x: as above
        - ('specific_domain_name_1', 'specific_domain_name_2'): pair of specific domains, sorted alphabetically
        - (sequence_identity, alignment_length): sequence identity and alignment length of the domain pair"""

    aa_list = [a for a in open(pseudo_aa, 'r').read().split('\n')[1:] if a != '']
    AA = [m.split('\t')[1] for m in aa_list]

    for aa in AA:
        domain_w_aa = []
        for i in BGCs.keys():
            cluster_i = BGCs[i]
            cluster_i_aa = cluster_i.keys()
            if aa in cluster_i_aa:
                domain_w_aa = domain_w_aa + cluster_i[aa]

        aa_pair_list = [tuple(sorted(list(i))) for i in list(itertools.combinations(domain_w_aa, 2))]
        aa_pair_dict = {}
        for p in aa_pair_list:
            compd1a = p[0].split('_')
            index1 = ['A'+ compd1a.pop()]
            compd1 = '_'.join(compd1a)
            compd2a = p[1].split('_')
            index2 = ['A'+ compd2a.pop()]
            compd2 = '_'.join(compd2a)
            
            compd2 = compd2.replace(" ","")
            compd1 = compd1.replace(" ","")
            
            #print seq_simScore
            #print compd1
            #print compd2
            #print index1
            #print index2
            
            score = get_distance(DM = seq_simScore, compd1 = compd1, compd2 = compd2, index1 = index1,                                  
            index2 = index2)[0]
            aa_pair_dict[p] = (score, 0)

        DMS[aa] = aa_pair_dict

    
    return BGCs, DMS, cluster_seq



def calculate_GK(A, B, nbhood): #nbhood = 3, can be changed
    # calculate the Goodman-Kruskal gamma index
    # Measures the synteny conservation (in this case, of substrate speficifities)
    #
    GK = 0.
    if len(set(A) & set(B)) > 1:
        pairsA = set( [(A[i],A[j]) for i in xrange(len(A)-1) for j in 
                       xrange(i+1, (i+nbhood if i+nbhood < len(A) else len(A)))] )
        pairsB = set( [(B[i],B[j]) for i in xrange(len(B)-1) for j in 
                       xrange(i+1, (i+nbhood if i+nbhood < len(B) else len(B)))] )
        allPairs = set(list(pairsA) + list(pairsB))
        Ns, Nr = 0.,0.
        for p in allPairs:
            if p in pairsA and p in pairsB: Ns += 1
            elif p in pairsA and tuple(p[::-1]) in pairsB: Nr += 1
            elif tuple(p[::-1]) in pairsA and p in pairsB: Nr += 1
            else: pass
        if (Nr + Ns) == 0:
            gamma = 0
        else:
            gamma = (Ns-Nr) / (Nr+Ns)
        GK = (1+gamma)/2.
    return GK


def cluster_distance(A, B, nbhood):
    #key is name of GC, values is list of specificity names
    clusterA = BGCs[A] # dictionary where keys are specificities, and values 
        # are lists of specificity labels that map to a specific sequence in 
        # the DMS variable
    clusterB = BGCs[B]

    specificities_A = set(clusterA.keys())
    specificities_B = set(clusterB.keys())
    intersect = specificities_A.intersection(specificities_B)
    
    # JACCARD INDEX
    Jaccard = len(intersect) / float(len(specificities_A.union(specificities_B)))


    #DDS: The difference in specificities' sequences between cluster
    #S: Max occurence of each specificity

    DDS,S = 0,0
    SumDistance = 0
    pair = ""

    # elements in either set but not in both: union - intersect
    not_intersect = specificities_A ^ specificities_B

    # sum as many copies of the unshared specificity there are, for each
    for unshared_specificity in not_intersect:
        dom_set = []
        try:
            dom_set = clusterA[unshared_specificity]
        except KeyError:
            dom_set = clusterB[unshared_specificity]

        DDS += len(dom_set) 
        S += len(dom_set)

    # compare sequence identity of shared specificities
    for shared_specificity in intersect:
        seta = clusterA[shared_specificity]
        setb = clusterB[shared_specificity]

        if len(seta+setb) == 2: #The specificity occurs only once in both clusters
            pair = tuple(sorted([seta[0],setb[0]]))

            try:
                SumDistance = 1-DMS[shared_specificity][pair][0]
                # print 'SumDistance1', SumDistance

            except KeyError:
                print "KeyError on", pair
                print(shared_specificity)
                sys.exit()

            S += max(len(seta),len(setb))
            DDS += SumDistance

        else: #The specificity occurs more than once in both clusters
            # accumulated_distance = 0

            DistanceMatrix = [[1 for col in range(len(setb))] for row in range(len(seta))]
            # print DistanceMatrix
            for domsa in range(len(seta)):
                for domsb in range(len(setb)):
                    pair = tuple(sorted([seta[domsa], setb[domsb]]))
                    try:
                        Similarity = DMS[shared_specificity][pair][0]
                        # print Similarity
                    except KeyError:
                        print "KeyError on (case 2)", pair
                        print(shared_specificity)

                    seq_dist = 1-Similarity
                    DistanceMatrix[domsa][domsb] = seq_dist


            #Only use the best scoring pairs
            Hungarian = Munkres()
            BestIndexes = Hungarian.compute(DistanceMatrix)
            accumulated_distance = sum([DistanceMatrix[bi[0]][bi[1]] for bi in BestIndexes])
            SumDistance = (abs(len(seta)-len(setb)) + accumulated_distance)  #diff in abundance + sequence distance

            S += max(len(seta),len(setb))
            DDS += SumDistance

    DDS /= float(S)
    DDS = 1-DDS #transform into similarity


    #  calculate the Goodman-Kruskal gamma index
    A_pseudo_seq = list(cluster_seq[A])
    B_pseudo_seq = list(cluster_seq[B])
    Ar = [item for item in A_pseudo_seq]
    Ar.reverse()
    GK = max([calculate_GK(A_pseudo_seq, B_pseudo_seq, nbhood = nbhood), calculate_GK(Ar, B_pseudo_seq, nbhood = nbhood)])


    Distance = 1 - (Jaccardw * Jaccard) - (DDSw * DDS) - (GKw * GK)
    Similarity_score = (Jaccardw * Jaccard) + (DDSw * DDS) + (GKw * GK)
    # Similarity_score = 1 - DDS
    if Distance < 0:
        print "negative distance", Distance, "DDS", DDS, pair
        print "Probably a rounding issue"
        print "Distance is set to 0 for these clusters"
        Distance = 0
    print A, B, Jaccard, GK, DDS
    return Similarity_score


def generate_distance_matrix(cluster_list, scale, nbhood, outfile):

    Dist = [[0 for x in range(len(cluster_list))] for y in range(len(cluster_list))]
    df = pd.DataFrame(Dist, index = cluster_list, columns = cluster_list)
    pairs = list(itertools.combinations(cluster_list, 2))

    for p in pairs:
        sim_score = cluster_distance(A = p[1],  B=p[0], nbhood=nbhood)
        rowname = p[1]
        colname = p[0]
        df.ix[rowname, colname] = 1- sim_score/scale

    df.to_csv(outfile)

    return df

## run the script to calculate the distance

global Jaccardw
global GKw
global DDSw
Jaccardw = 0.5
GKw = 0.25
DDSw = 0.25

seq_simScore = "data/sequence_similarity_score.txt"
cluster = "data/glycopeptides_domain_pseudoSeq.fasta"   # fasta where sequences are aa specificities; 
pseudo_aa = "data/Clade_PhyloNode_Adomain_all_for_alignment.txt" 

outfile = "dist.txt"
tree_outfile = "tree_upgmma_all_clade.nwk"

BGCs, DMS, cluster_seq = generate_BGCs_DMS(seq_simScore = seq_simScore, cluster = cluster, pseudo_aa = pseudo_aa)

cluster_list = cluster_seq.keys()
dist_score_assembly_line = generate_distance_matrix(cluster_list, scale = 1, nbhood = 3, outfile = outfile)

#-- Plot the tree
import Bio.Phylo
from Bio.Phylo.TreeConstruction import _Matrix, _DistanceMatrix
from Bio.Phylo.TreeConstruction import _Matrix, _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
names = cluster_list
score = [s for s in open(outfile, 'r').read().split('\n')[1:] if s != '']
matrix = []
for i in range(len(score)):
    input_i = score[i].split(',')[1:(i+2)]
    input_i_int = [float(n) for n in input_i]
    matrix.append(input_i_int)
m = _DistanceMatrix(names, matrix)

constructor = DistanceTreeConstructor()
tree1 = constructor.upgma(m)
Bio.Phylo.draw_ascii(tree1)
Bio.Phylo.write(tree1, tree_outfile, 'newick')
