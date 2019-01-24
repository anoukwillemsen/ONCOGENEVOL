#!/usr/bin/env python2.7
""" calc ML of best model for reshuffled and realigned subsets """ 

from Bio import AlignIO
from Bio import Align
from Bio import SeqIO
from Bio import Phylo
from Bio.Align import Applications
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from StringIO import StringIO
import random     # shuffle
import subprocess # open pipe 
import sys        # check platform
import os
import re         # regular expressions, perl-style

n_iters = 100
suffix = str(random.randint(1,1e10))

def phyml_ml (treefilename):
    """ calculates ML tree and branch lengths - estimates gamma, pinvar, and assumes model eq. freqs """
    os.system ("phyml -i mc.paml -d aa -m LG -b 0 -c 4 -v e" + treefilename) # run phyML
    tree = Phylo.read("mc.paml_phyml_tree", "newick") # read ML tree
    strout = str(tree.total_branch_length() - 1.) + " "        # store brlen to print as output
    for line in open("mc.paml_phyml_stats"):          # search for lnL value
        if ". Log-likelihood:" in line:
            lnl = str(line.split(" ")[-1].strip())      # store last column using string.split() 
        if "- Gamma shape parameter:" in line:
            gamma = str(line.split(" ")[-1].strip()) + " "  # store last column using string.split()  

    return strout, gamma, lnl

def align (alignment):
    """ muscle alignment using pipes (not files) """
    cline = Align.Applications.MuscleCommandline(clwstrict=False, maxhours=1.) # muscle command line 
    child = subprocess.Popen(str(cline), 
                    stdin=subprocess.PIPE, 
                    stdout=subprocess.PIPE, 
#                    stderr=subprocess.PIPE, 
                    shell=(sys.platform!="win32"))

    SeqIO.write(alignment, child.stdin, "fasta") # after this muscle is still waiting, so...
    child.stdin.close() # ... we must close the handle by hand, which will then make muscle start calculations
    return AlignIO.read(child.stdout, "fasta")  # read from stdout as a fasta file

def write_to_phylip (align, filename):
    phy = open (filename, "w", 0)
    phy.write (str (len (align)) + "  " + str (len (align[0])) + "\n")
    for rec in align:
        phy.write (str (rec.id) + "   " + str (rec.seq) + "\n")

    phy.close()


os.system ("cp -a ../000.template_tmp_dir tmp." + suffix)
os.chdir ("./tmp." + suffix)
outfile = open ("table_LnL.txt", "w", 0)

filename = "P1-max.ALL.fas"
align_ALL = AlignIO.read (filename, "fasta")
seqnames = [seq.id for seq in align_ALL]

idx = [seqnames.index(i) for i in ["CUT_HPV29_E5_BETA", "CUT_HPV81_E5_BETA", "GW_HPV91_E5_DELTA", "GW_HPV7_E5_DELTA", "GW_HPV11_E5_GAMMA_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "MUC_HPV70_E5_ALPHA1", "MUC_HPV18_E5_ALPHA1", "MUC_HPV34_E5_ALPHA2", "MUC_HPV16_E5_ALPHA2", "MUC_MfPV5_E5_EPSILON_DSETA", "MUC_MfPV11_E5_EPSILON_DSETA"]]
align_fix1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx]) #red
idx = [seqnames.index(i) for i in ["AaPV1", "BgPV1", "BPV13", "BPV14", "BPV1", "BPV2", "CcaPV1", "OaPV1", "OaPV2", "OvPV1", "RalPV1", "RtPV1"]]
align_var = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  #blue

colindex = range (len (align_var[0])) # each alignment[] has record.seq and record.id (here is record.seq)

for iteration in range (n_iters):
    random.shuffle (colindex) # align_var will have columns changed 
    # create an alignment with zero length and same ids as align_var
    align_shuffled = Align.MultipleSeqAlignment([SeqRecord(Seq("",generic_dna), id=i.id) for i in align_var])
    # column is only one col (that's why the slice [,1:2] instead of [,1] )
    for column in [align_var[:,i:(i+1)] for i in colindex]:
        align_shuffled += column    # one-liner for this ? (maybe lambda for sum?)
    # we could also have used alignment.extend()
    align_ALL = Align.MultipleSeqAlignment([i for i in align_shuffled] + [i for i in align_fix1])

    str_out = ""
    model_out = ""
    lnl_out = ""

    # align simulated data
    align_ALL = align (align_ALL)
    seqnames = [seq.id for seq in align_ALL]
    
    idx = [seqnames.index(i) for i in ["CUT_HPV29_E5_BETA", "CUT_HPV81_E5_BETA", "GW_HPV91_E5_DELTA", "GW_HPV7_E5_DELTA", "GW_HPV11_E5_GAMMA_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "MUC_HPV70_E5_ALPHA1", "MUC_HPV18_E5_ALPHA1", "MUC_HPV34_E5_ALPHA2", "MUC_HPV16_E5_ALPHA2", "MUC_MfPV5_E5_EPSILON_DSETA", "MUC_MfPV11_E5_EPSILON_DSETA"]]
    align_red = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_red, "mc.paml")
    strout_red, gamma, lnl_red = phyml_ml ("c50.PP.red.tree")
    str_out += strout_red
    model_out += gamma
    lnl_out += lnl_red + " "

#    idx = [seqnames.index(i) for i in ["SfPV1", "OcPV1"]]
#    align_yellow = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
#    write_to_phylip (align_yellow, "mc.paml")
#    strout_yellow, gamma, lnl_yellow = phyml_ml ("c50.PP.yellow.tree")
#    str_out += strout_yellow
#    model_out += gamma
#    lnl_out += lnl_yellow + " "

    idx = [seqnames.index(i) for i in ["AaPV1", "BgPV1", "BPV13", "BPV14", "BPV1", "BPV2", "CcaPV1", "OaPV1", "OaPV2", "OvPV1", "RalPV1", "RtPV1"]]
    align_blue = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_blue, "mc.paml")
    strout_blue, gamma, lnl_blue = phyml_ml ("c50.PP.blue.tree")
    str_out += strout_blue
    model_out += gamma 
    lnl_out += lnl_blue + " "

######################################

######################################
# ALL together
    write_to_phylip (align_ALL, "mc.paml")
    strout_ALL, gamma, lnl_ALL = phyml_ml ("c50.PP.ALL.tree")
    str_out += strout_ALL
    model_out += gamma
    lnl_out += lnl_ALL + " "   
######################################

    # output: ML tree lengths
    # outorder: treelen_red treelen_blue treelen_ALL lnl_red lnl_blue lnl_ALL ali_len

    outfile.write (str_out + lnl_out + str(len (align_ALL[0])) + "\n")

outfile.close()
os.system ("mv table_LnL.txt ../tableLnL_baliphy_" + suffix + ".txt")
os.chdir  ("../")
os.system ("rm -rf tmp." + suffix)
