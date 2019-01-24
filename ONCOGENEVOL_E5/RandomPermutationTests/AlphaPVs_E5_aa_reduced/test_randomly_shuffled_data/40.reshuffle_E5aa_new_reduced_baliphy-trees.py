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

idx = [seqnames.index(i) for i in ["MUC_HPV70_E5_ALPHA1", "MUC_HPV39_E5_ALPHA1", "MUC_HPV18_E5_ALPHA1", "MUC_HPV85_E5_ALPHA1"]]
align_fix1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx]) #alpha1
idx = [seqnames.index(i) for i in ["MUC_HPV34_E5_ALPHA2", "MUC_HPV16_E5_ALPHA2", "MUC_HPV31_E5_ALPHA2", "MUC_HPV58_E5_ALPHA2"]]
align_fix2 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  #alpha2
idx = [seqnames.index(i) for i in ["CUT_HPV29_E5_BETA", "CUT_HPV106_E5_BETA", "CUT_HPV57_E5_BETA", "CUT_HPV81_E5_BETA"]]
align_fix3 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  #beta
idx = [seqnames.index(i) for i in ["GW_HPV11_E5_GAMMA_DELTA", "GW_HPV6_E5_GAMMA_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "GW_PpPV1_E5_GAMMA_DELTA"]]
align_fix4 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  #gammadelta
idx = [seqnames.index(i) for i in ["GW_HPV7_E5_DELTA", "GW_HPV40_E5_DELTA", "GW_HPV43_E5_DELTA", "GW_HPV91_E5_DELTA"]]
align_fix5 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  #delta
idx = [seqnames.index(i) for i in ["MUC_MfPV5_E5_EPSILON_DSETA", "MUC_MfPV11_E5_EPSILON_DSETA", "MUC_PhPV1_E5_EPSILON_DSETA", "MUC_MfPV4_E5_EPSILON_DSETA"]]
align_var = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  #epsilonzeta


colindex = range (len (align_var[0])) # each alignment[] has record.seq and record.id (here is record.seq)

for iteration in range (n_iters):
    random.shuffle (colindex) # align_var will have columns changed 
    # create an alignment with zero length and same ids as align_var
    align_shuffled = Align.MultipleSeqAlignment([SeqRecord(Seq("",generic_dna), id=i.id) for i in align_var])
    # column is only one col (that's why the slice [,1:2] instead of [,1] )
    for column in [align_var[:,i:(i+1)] for i in colindex]:
        align_shuffled += column    # one-liner for this ? (maybe lambda for sum?)
    # we could also have used alignment.extend()
    align_ALL = Align.MultipleSeqAlignment([i for i in align_shuffled] + [i for i in align_fix1] + [i for i in align_fix2] + [i for i in align_fix3] + [i for i in align_fix4] + [i for i in align_fix5])

    str_out = ""
    model_out = ""
    lnl_out = ""

    # align simulated data
    align_ALL = align (align_ALL)
    seqnames = [seq.id for seq in align_ALL]
    
    idx = [seqnames.index(i) for i in ["MUC_HPV70_E5_ALPHA1", "MUC_HPV39_E5_ALPHA1", "MUC_HPV18_E5_ALPHA1", "MUC_HPV85_E5_ALPHA1"]]
    align_alpha1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_alpha1, "mc.paml")
    strout_alpha1, gamma, lnl_alpha1 = phyml_ml ("c50.PP.alpha1.tree")
    str_out += strout_alpha1
    model_out += gamma
    lnl_out += lnl_alpha1 + " "

    idx = [seqnames.index(i) for i in ["MUC_HPV34_E5_ALPHA2", "MUC_HPV16_E5_ALPHA2", "MUC_HPV31_E5_ALPHA2", "MUC_HPV58_E5_ALPHA2"]]
    align_alpha2 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_alpha2, "mc.paml")
    strout_alpha2, gamma, lnl_alpha2 = phyml_ml ("c50.PP.alpha2.tree")
    str_out += strout_alpha2
    model_out += gamma
    lnl_out += lnl_alpha2 + " "

    idx = [seqnames.index(i) for i in ["CUT_HPV29_E5_BETA", "CUT_HPV106_E5_BETA", "CUT_HPV57_E5_BETA", "CUT_HPV81_E5_BETA"]]
    align_beta = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_beta, "mc.paml")
    strout_beta, gamma, lnl_beta = phyml_ml ("c50.PP.beta.tree")
    str_out += strout_beta
    model_out += gamma
    lnl_out += lnl_beta + " "

    idx = [seqnames.index(i) for i in ["GW_HPV11_E5_GAMMA_DELTA", "GW_HPV6_E5_GAMMA_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "GW_PpPV1_E5_GAMMA_DELTA"]]
    align_gammadelta = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_gammadelta, "mc.paml")
    strout_gammadelta, gamma, lnl_gammadelta = phyml_ml ("c50.PP.gammadelta.tree")
    str_out += strout_gammadelta
    model_out += gamma 
    lnl_out += lnl_gammadelta + " "

    idx = [seqnames.index(i) for i in ["GW_HPV7_E5_DELTA", "GW_HPV40_E5_DELTA", "GW_HPV43_E5_DELTA", "GW_HPV91_E5_DELTA"]]
    align_delta = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_delta, "mc.paml")
    strout_delta, gamma, lnl_delta = phyml_ml ("c50.PP.delta.tree")
    str_out += strout_delta
    model_out += gamma 
    lnl_out += lnl_delta + " "


    idx = [seqnames.index(i) for i in ["MUC_MfPV5_E5_EPSILON_DSETA", "MUC_MfPV11_E5_EPSILON_DSETA", "MUC_PhPV1_E5_EPSILON_DSETA", "MUC_MfPV4_E5_EPSILON_DSETA"]]
    align_epsilonzeta = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_epsilonzeta, "mc.paml")
    strout_epsilonzeta, gamma, lnl_epsilonzeta = phyml_ml ("c50.PP.epsilonzeta.tree")
    str_out += strout_epsilonzeta
    model_out += gamma 
    lnl_out += lnl_epsilonzeta + " "

######################################
# alternative IO hypotheses
#   alpha1-alpha2
    idx = [seqnames.index(i) for i in ["MUC_HPV70_E5_ALPHA1", "MUC_HPV39_E5_ALPHA1", "MUC_HPV18_E5_ALPHA1", "MUC_HPV85_E5_ALPHA1", "MUC_HPV34_E5_ALPHA2", "MUC_HPV16_E5_ALPHA2", "MUC_HPV31_E5_ALPHA2", "MUC_HPV58_E5_ALPHA2"]]
    align_H1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx]) 
    write_to_phylip (align_H1, "mc.paml")
    strout_H1, gamma, lnl_H1 = phyml_ml ("c50.PP.alpha1_alpha2.tree")
    str_out += strout_H1
    model_out += gamma
    lnl_out += lnl_H1 + " "

#   gammadelta-delta
    idx = [seqnames.index(i) for i in ["GW_HPV11_E5_GAMMA_DELTA", "GW_HPV6_E5_GAMMA_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "GW_PpPV1_E5_GAMMA_DELTA", "GW_HPV7_E5_DELTA", "GW_HPV40_E5_DELTA", "GW_HPV43_E5_DELTA", "GW_HPV91_E5_DELTA"]]
    align_H2 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_H2, "mc.paml")
    strout_H2, gamma, lnl_H2 = phyml_ml ("c50.PP.gammadelta_delta.tree")
    str_out += strout_H2
    model_out += gamma
    lnl_out += lnl_H2 + " "

#   alpha1-alpha2-gammadelta-delta
    idx = [seqnames.index(i) for i in ["MUC_HPV70_E5_ALPHA1", "MUC_HPV39_E5_ALPHA1", "MUC_HPV18_E5_ALPHA1", "MUC_HPV85_E5_ALPHA1", "MUC_HPV34_E5_ALPHA2", "MUC_HPV16_E5_ALPHA2", "MUC_HPV31_E5_ALPHA2", "MUC_HPV58_E5_ALPHA2", "GW_HPV11_E5_GAMMA_DELTA", "GW_HPV6_E5_GAMMA_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "GW_PpPV1_E5_GAMMA_DELTA", "GW_HPV7_E5_DELTA", "GW_HPV40_E5_DELTA", "GW_HPV43_E5_DELTA", "GW_HPV91_E5_DELTA"]]
    align_H3 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_H3, "mc.paml")
    strout_H3, gamma, lnl_H3 = phyml_ml ("c50.PP.alpha1_alpha2_gammdelta_delta.tree")
    str_out += strout_H3
    model_out += gamma
    lnl_out += lnl_H3 + " "

#   alpha1-alpha2-gammadelta-delta-epsilonzeta
    idx = [seqnames.index(i) for i in ["MUC_HPV70_E5_ALPHA1", "MUC_HPV39_E5_ALPHA1", "MUC_HPV18_E5_ALPHA1", "MUC_HPV85_E5_ALPHA1", "MUC_HPV34_E5_ALPHA2", "MUC_HPV16_E5_ALPHA2", "MUC_HPV31_E5_ALPHA2", "MUC_HPV58_E5_ALPHA2", "GW_HPV11_E5_GAMMA_DELTA", "GW_HPV6_E5_GAMMA_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "GW_PpPV1_E5_GAMMA_DELTA", "GW_HPV7_E5_DELTA", "GW_HPV40_E5_DELTA", "GW_HPV43_E5_DELTA", "GW_HPV91_E5_DELTA", "MUC_MfPV5_E5_EPSILON_DSETA", "MUC_MfPV11_E5_EPSILON_DSETA", "MUC_PhPV1_E5_EPSILON_DSETA", "MUC_MfPV4_E5_EPSILON_DSETA"]]
    align_H4 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_H4, "mc.paml")
    strout_H4, gamma, lnl_H4 = phyml_ml ("c50.PP.alpha1_alpha2_gammadelta_delta_epsilonzeta.tree")
    str_out += strout_H4
    model_out += gamma
    lnl_out += lnl_H4 + " "

#   alpha1-alpha2-epsilonzeta
    idx = [seqnames.index(i) for i in ["MUC_HPV70_E5_ALPHA1", "MUC_HPV39_E5_ALPHA1", "MUC_HPV18_E5_ALPHA1", "MUC_HPV85_E5_ALPHA1", "MUC_HPV34_E5_ALPHA2", "MUC_HPV16_E5_ALPHA2", "MUC_HPV31_E5_ALPHA2", "MUC_HPV58_E5_ALPHA2", "MUC_MfPV5_E5_EPSILON_DSETA", "MUC_MfPV11_E5_EPSILON_DSETA", "MUC_PhPV1_E5_EPSILON_DSETA", "MUC_MfPV4_E5_EPSILON_DSETA"]]
    align_H5 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx]) 
    write_to_phylip (align_H5, "mc.paml")
    strout_H5, gamma, lnl_H5 = phyml_ml ("c50.PP.alpha1_alpha2_epsilonzeta.tree")
    str_out += strout_H5
    model_out += gamma
    lnl_out += lnl_H5 + " "

######################################
# ALL together
    write_to_phylip (align_ALL, "mc.paml")
    strout_ALL, gamma, lnl_ALL = phyml_ml ("c50.PP.ALL.tree")
    str_out += strout_ALL
    model_out += gamma
    lnl_out += lnl_ALL + " "   
######################################

    # output: ML tree lengths
    # outorder: treelen_alpha1 treelen_alpha2 treelen_beta treelen_gammadelta treelen_delta treelen_epsilon_dseta treelen_part_H1 treelen_part_H2 treelen_part_H3 treelen_part_H4 treelen_part_H5 treelen_ALL lnl_alpha1 lnl_alpha2 lnl_beta lnl_gammadelta treelen_delta lnl_epsilon_dseta lnl_part_H1 lnl_part_H2 lnl_part_H3 lnl_part_H4 lnl_part_H5 lnl_ALL ali_len

    outfile.write (str_out + lnl_out + str(len (align_ALL[0])) + "\n")

outfile.close()
os.system ("mv table_LnL.txt ../tableLnL_baliphy_" + suffix + ".txt")
os.chdir  ("../")
os.system ("rm -rf tmp." + suffix)
