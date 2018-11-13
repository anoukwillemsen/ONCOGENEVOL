#!/usr/bin/env python2.7
""" Read 4 seqs alignments, shuffles one of these groups and realigns """

## BUG: the total_tree_length() includes a nonexistent internal branch_length of 1.0 
from Bio import AlignIO
from Bio import Align
from Bio.Align import Applications
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import Phylo
import random     # shuffle
import subprocess # open pipe 
import sys        # check platform
import os

n_iters = 10

filename = "E6_AA_red3_selAll_all4Manatees_mafftali_newHead.fas"
align_ALL = AlignIO.read (filename, "fasta")
seqnames = [seq.id for seq in align_ALL]

idx = [seqnames.index(i) for i in ["FgPV1", "FlPV1", "PaPV1", "FcPV1", "CmPV1", "CcPV1"]]
align_fix1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # AvesTurtles
idx = [seqnames.index(i) for i in ["BpPV1", "HPV107", "RtiPV1", "FcaPV4", "VvPV1", "HPV133"]]
align_fix2 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # BetaGamma
idx = [seqnames.index(i) for i in ["EsPV3", "EhPV1", "SfPV1", "ElPV1", "CPV1", "PlpPV1"]]
align_fix3 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # Lambda
idx = [seqnames.index(i) for i in ["TmPV2", "TmPV3", "TmPV4", "TmPV1"]]
align_fix4 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # Manatees
idx = [seqnames.index(i) for i in ["BPV8", "CdPV2", "OaPV2", "OvPV1", "EcPV6", "EcPV4"]]
align_fix5 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # PerissoArtio
idx = [seqnames.index(i) for i in ["PphPV4", "MrPV1", "AgPV1", "HPV77", "HPV32", "HPV70"]]
align_var = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # Alpha

# DEBUG: small alignments
#align_fix = align_fix[:,:100]
#align_var = align_var[:,:100]

filesuffix = str(random.randint(1,1e10))

index = range (len (align_var[0])) # each alignment[] has record.seq and record.id (here is record.seq)
cline = Align.Applications.MuscleCommandline(clwstrict=False, maxhours=1.) # muscle command line 
outfile = open ("res_table_" + filesuffix + ".txt", "w", 0)

for iteration in range (n_iters):
    random.shuffle (index) # align_var will have columns changed 
    str_out = str(iteration) + " "    # output table row
    # create an alignment with zero length and same ids as align_var
    align_shuffled = Align.MultipleSeqAlignment([SeqRecord(Seq("",generic_dna), id=i.id) for i in align_var])
    # column is only one col (that's why the slice [,1:2] instead of [,1] )
    for column in [align_var[:,i:(i+1)] for i in index]:
        align_shuffled += column    # one-liner for this ? (maybe lambda for sum?)
    
    # we could also have used alignment.extend()
    align_input = Align.MultipleSeqAlignment([i for i in align_shuffled] + [i for i in align_fix1] + [i for i in align_fix2] + [i for i in align_fix3] + [i for i in align_fix4] + [i for i in align_fix5])

    # ML tree and branch lengths as Theobald
    SeqIO.write(align_input, "align" + filesuffix + ".phy", "phylip") 
    os.system ("phyml -i align" + filesuffix + ".phy -d aa -m LG -o tlr -b 0 -c 1") # run phyML
    tree = Phylo.read("align" + filesuffix + ".phy_phyml_tree", "newick") # read ML tree
    str_out += str(tree.total_branch_length()) + " "            # store brlen to print as output
    for line in open("align" + filesuffix + ".phy_phyml_stats"):          # search for lnL value
        if ". Log-likelihood:" in line:
            str_out += str(line.split(" ")[-1].strip()) + " "           # store last column  
            break

    child = subprocess.Popen(str(cline),  # already starts muscle, which is now waiting input from stdin
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
#                stderr=subprocess.PIPE,
                shell=(sys.platform!="win32"))

    SeqIO.write(align_input, child.stdin, "fasta") # after this muscle is still waiting, so...
    child.stdin.close() # ... we must close the handle by hand, which will then make muscle start calculations
    alignment = AlignIO.read(child.stdout, "fasta")  # read from stdout as a fasta file

    # ML tree and branch lengths on properly aligned sequences
    SeqIO.write(alignment, "align" + filesuffix + ".phy", "phylip") 
    os.system ("phyml -i align" + filesuffix + ".phy -d aa -m LG -o tlr -b 0 -c 1") # run phyML
    tree = Phylo.read("align" + filesuffix + ".phy_phyml_tree", "newick") # read ML tree
    str_out += str(tree.total_branch_length()) + " "            # store brlen to print as output
    for line in open("align" + filesuffix + ".phy_phyml_stats"):    # search for lnL value
        if ". Log-likelihood:" in line:
            str_out += str(line.split(" ")[-1].strip()) + " "           # store last column  
            break

    str_out += str (len (alignment[0])) 
    outfile.write (str_out + "\n")

outfile.close()
