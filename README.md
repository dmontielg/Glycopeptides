# Glycopeptide antibiotics Prediction in non-ribosomal peptide-synthetase (NRPS)

### Installation
  Python 2.7 (Anaconda package recoomended)
  MUSCLE
  hmmsearch3
  RaxML 8

### Usage

python start.py

The following options will appear:
    Prediction of putative Glycopeptide Biosynthetic Gene Cluster
        1. Retrieve CDS from genbank
        2. Predict putative GBGC
        3. Extract AT domain from sequences
        4. Perform MSA with MUSCLE
        5. Maximum-Likelihood with RaxML

1. Folder where your putative glycopeptide gene clusters are (genbank) will extract the Coding Sequence and gene if there is.
2. Give the genes from the putative glycopeptide gene clusters and will try to predict if there is a potential glycopeptide gene cluster.
3. Uses hmmsearch for profile domain using Adenylation and Thiolation domains.
4. Receive a fasta file and use MUSCLE for Multiple Sequence Alingment.
5. Use an aligned file in Phylib format and execute RaxML for Maximum Likelihood of 100 bootstrap with the model WAG+I+G+F.
