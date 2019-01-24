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

filename = "original_red-blue.phy"
align_ALL = AlignIO.read (filename, "phylip-relaxed")
seqnames = [seq.id for seq in align_ALL]

idx = [seqnames.index(i) for i in ["MUC_HPV16_E5_ALPHA", "MUC_HPV18_E5_ALPHA", "MUC_HPV31_E5_ALPHA", "MUC_HPV33_E5_ALPHA", "MUC_HPV34_E5_ALPHA", "MUC_HPV35_E5_ALPHA", "MUC_HPV39_E5_ALPHA", "MUC_HPV45_E5_ALPHA", "MUC_HPV52_E5_ALPHA", "MUC_HPV58_E5_ALPHA", "MUC_HPV59_E5_ALPHA", "MUC_HPV67_E5_ALPHA", "MUC_HPV68_E5_ALPHA", "MUC_HPV70_E5_ALPHA", "MUC_HPV73_E5_ALPHA", "MUC_HPV85_E5_ALPHA", "MUC_HPV97_E5_ALPHA", "CUT_CgPV1_E5_BETA", "CUT_HPV2_E5_BETA", "CUT_HPV3_E5_BETA", "CUT_HPV10_E5_BETA", "CUT_HPV27_E5_BETA", "CUT_HPV28_E5_BETA", "CUT_HPV29_E5_BETA", "CUT_HPV57_E5_BETA", "CUT_HPV61_E5_BETA", "CUT_HPV62_E5_BETA", "CUT_HPV71_E5_BETA", "CUT_HPV72_E5_BETA", "CUT_HPV77_E5_BETA", "CUT_HPV78_E5_BETA", "CUT_HPV81_E5_BETA", "CUT_HPV83_E5_BETA", "CUT_HPV84_E5_BETA", "CUT_HPV86_E5_BETA", "CUT_HPV87_E5_BETA", "CUT_HPV89_E5_BETA", "CUT_HPV90_E5_BETA", "CUT_HPV94_E5_BETA", "CUT_HPV102_E5_BETA", "CUT_HPV106_E5_BETA", "CUT_HPV114_E5_BETA", "CUT_HPV117_E5_BETA", "CUT_HPV125_E5_BETA", "CUT_HPV160_E5_BETA", "MUC_MfPV3_E5_EPSILON_DSETA", "MUC_MfPV4_E5_EPSILON_DSETA", "MUC_MfPV5_E5_EPSILON_DSETA", "MUC_MfPV6_E5_EPSILON_DSETA", "MUC_MfPV7_E5_EPSILON_DSETA", "MUC_MfPV8_E5_EPSILON_DSETA", "MUC_MfPV9_E5_EPSILON_DSETA", "MUC_MfPV10_E5_EPSILON_DSETA", "MUC_MfPV11_E5_EPSILON_DSETA", "MUC_MfuPV1_E5_EPSILON_DSETA", "MUC_MmPV1_E5_EPSILON_DSETA", "MUC_PhPV1_E5_EPSILON_DSETA", "GW_HPV6_E5_GAMMA_DELTA", "GW_HPV7_E5_DELTA", "GW_HPV11_E5_GAMMA_DELTA", "GW_HPV13_E5_GAMMA_DELTA", "GW_HPV40_E5_DELTA", "GW_HPV43_E5_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "GW_HPV54_E5", "GW_HPV74_E5_GAMMA_DELTA", "GW_HPV91_E5_DELTA", "GW_PpPV1_E5_GAMMA_DELTA", "GW_CchPV1_E5_GAMMA_DELTA"]]
align_fix1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx]) #red
#idx = [seqnames.index(i) for i in ["SfPV1", "OcPV1"]]
#align_fix2 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  #yellow
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
    
    idx = [seqnames.index(i) for i in ["MUC_HPV16_E5_ALPHA", "MUC_HPV18_E5_ALPHA", "MUC_HPV31_E5_ALPHA", "MUC_HPV33_E5_ALPHA", "MUC_HPV34_E5_ALPHA", "MUC_HPV35_E5_ALPHA", "MUC_HPV39_E5_ALPHA", "MUC_HPV45_E5_ALPHA", "MUC_HPV52_E5_ALPHA", "MUC_HPV58_E5_ALPHA", "MUC_HPV59_E5_ALPHA", "MUC_HPV67_E5_ALPHA", "MUC_HPV68_E5_ALPHA", "MUC_HPV70_E5_ALPHA", "MUC_HPV73_E5_ALPHA", "MUC_HPV85_E5_ALPHA", "MUC_HPV97_E5_ALPHA", "CUT_CgPV1_E5_BETA", "CUT_HPV2_E5_BETA", "CUT_HPV3_E5_BETA", "CUT_HPV10_E5_BETA", "CUT_HPV27_E5_BETA", "CUT_HPV28_E5_BETA", "CUT_HPV29_E5_BETA", "CUT_HPV57_E5_BETA", "CUT_HPV61_E5_BETA", "CUT_HPV62_E5_BETA", "CUT_HPV71_E5_BETA", "CUT_HPV72_E5_BETA", "CUT_HPV77_E5_BETA", "CUT_HPV78_E5_BETA", "CUT_HPV81_E5_BETA", "CUT_HPV83_E5_BETA", "CUT_HPV84_E5_BETA", "CUT_HPV86_E5_BETA", "CUT_HPV87_E5_BETA", "CUT_HPV89_E5_BETA", "CUT_HPV90_E5_BETA", "CUT_HPV94_E5_BETA", "CUT_HPV102_E5_BETA", "CUT_HPV106_E5_BETA", "CUT_HPV114_E5_BETA", "CUT_HPV117_E5_BETA", "CUT_HPV125_E5_BETA", "CUT_HPV160_E5_BETA", "MUC_MfPV3_E5_EPSILON_DSETA", "MUC_MfPV4_E5_EPSILON_DSETA", "MUC_MfPV5_E5_EPSILON_DSETA", "MUC_MfPV6_E5_EPSILON_DSETA", "MUC_MfPV7_E5_EPSILON_DSETA", "MUC_MfPV8_E5_EPSILON_DSETA", "MUC_MfPV9_E5_EPSILON_DSETA", "MUC_MfPV10_E5_EPSILON_DSETA", "MUC_MfPV11_E5_EPSILON_DSETA", "MUC_MfuPV1_E5_EPSILON_DSETA", "MUC_MmPV1_E5_EPSILON_DSETA", "MUC_PhPV1_E5_EPSILON_DSETA", "GW_HPV6_E5_GAMMA_DELTA", "GW_HPV7_E5_DELTA", "GW_HPV11_E5_GAMMA_DELTA", "GW_HPV13_E5_GAMMA_DELTA", "GW_HPV40_E5_DELTA", "GW_HPV43_E5_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "GW_HPV54_E5", "GW_HPV74_E5_GAMMA_DELTA", "GW_HPV91_E5_DELTA", "GW_PpPV1_E5_GAMMA_DELTA", "GW_CchPV1_E5_GAMMA_DELTA"]]
    align_red = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_red, "mc.paml")
    strout_red, gamma, lnl_red = phyml_ml ("original_red.phy_phyml_tree")
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
    strout_blue, gamma, lnl_blue = phyml_ml ("original_blue.phy_phyml_tree")
    str_out += strout_blue
    model_out += gamma 
    lnl_out += lnl_blue + " "

######################################
## alternative IO hypotheses
##   red-blue
#    idx = [seqnames.index(i) for i in ["MUC_HPV16_E5_ALPHA", "MUC_HPV18_E5_ALPHA", "MUC_HPV31_E5_ALPHA", "MUC_HPV33_E5_ALPHA", "MUC_HPV34_E5_ALPHA", "MUC_HPV35_E5_ALPHA", "MUC_HPV39_E5_ALPHA", "MUC_HPV45_E5_ALPHA", "MUC_HPV52_E5_ALPHA", "MUC_HPV58_E5_ALPHA", "MUC_HPV59_E5_ALPHA", "MUC_HPV67_E5_ALPHA", "MUC_HPV68_E5_ALPHA", "MUC_HPV70_E5_ALPHA", "MUC_HPV73_E5_ALPHA", "MUC_HPV85_E5_ALPHA", "MUC_HPV97_E5_ALPHA", "CUT_CgPV1_E5_BETA", "CUT_HPV2_E5_BETA", "CUT_HPV3_E5_BETA", "CUT_HPV10_E5_BETA", "CUT_HPV27_E5_BETA", "CUT_HPV28_E5_BETA", "CUT_HPV29_E5_BETA", "CUT_HPV57_E5_BETA", "CUT_HPV61_E5_BETA", "CUT_HPV62_E5_BETA", "CUT_HPV71_E5_BETA", "CUT_HPV72_E5_BETA", "CUT_HPV77_E5_BETA", "CUT_HPV78_E5_BETA", "CUT_HPV81_E5_BETA", "CUT_HPV83_E5_BETA", "CUT_HPV84_E5_BETA", "CUT_HPV86_E5_BETA", "CUT_HPV87_E5_BETA", "CUT_HPV89_E5_BETA", "CUT_HPV90_E5_BETA", "CUT_HPV94_E5_BETA", "CUT_HPV102_E5_BETA", "CUT_HPV106_E5_BETA", "CUT_HPV114_E5_BETA", "CUT_HPV117_E5_BETA", "CUT_HPV125_E5_BETA", "CUT_HPV160_E5_BETA", "MUC_MfPV3_E5_EPSILON_DSETA", "MUC_MfPV4_E5_EPSILON_DSETA", "MUC_MfPV5_E5_EPSILON_DSETA", "MUC_MfPV6_E5_EPSILON_DSETA", "MUC_MfPV7_E5_EPSILON_DSETA", "MUC_MfPV8_E5_EPSILON_DSETA", "MUC_MfPV9_E5_EPSILON_DSETA", "MUC_MfPV10_E5_EPSILON_DSETA", "MUC_MfPV11_E5_EPSILON_DSETA", "MUC_MfuPV1_E5_EPSILON_DSETA", "MUC_MmPV1_E5_EPSILON_DSETA", "MUC_PhPV1_E5_EPSILON_DSETA", "GW_HPV6_E5_GAMMA_DELTA", "GW_HPV7_E5_DELTA", "GW_HPV11_E5_GAMMA_DELTA", "GW_HPV13_E5_GAMMA_DELTA", "GW_HPV40_E5_DELTA", "GW_HPV43_E5_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "GW_HPV54_E5", "GW_HPV74_E5_GAMMA_DELTA", "GW_HPV91_E5_DELTA", "GW_PpPV1_E5_GAMMA_DELTA", "GW_CchPV1_E5_GAMMA_DELTA", "AaPV1", "BgPV1", "BPV13", "BPV14", "BPV1", "BPV2", "CcaPV1", "OaPV1", "OaPV2", "OvPV1", "RalPV1", "RtPV1"]]
#    align_H1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx]) 
#    write_to_phylip (align_H1, "mc.paml")
#    strout_H1, gamma, lnl_H1 = phyml_ml ("c50.PP.red-blue.tree")
#    str_out += strout_H1
#    model_out += gamma
#    lnl_out += lnl_H1 + " "

##   blue-yellow
#    idx = [seqnames.index(i) for i in ["AaPV1", "BgPV1", "BPV13", "BPV14", "BPV1", "BPV2", "CcaPV1", "OaPV1", "OaPV2", "OvPV1", "RalPV1", "RtPV1", "SfPV1", "OcPV1"]]
#    align_H2 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
#    write_to_phylip (align_H2, "mc.paml")
#    strout_H2, gamma, lnl_H2 = phyml_ml ("c50.PP.blue-yellow.tree")
#    str_out += strout_H2
#    model_out += gamma
#    lnl_out += lnl_H2 + " "

######################################
# ALL together
    write_to_phylip (align_ALL, "mc.paml")
    strout_ALL, gamma, lnl_ALL = phyml_ml ("original_red-blue.phy_phyml_tree")
    str_out += strout_ALL
    model_out += gamma
    lnl_out += lnl_ALL + " "   
######################################

    # output: ML tree lengths
    # outorder: treelen_red treelen_blue treelen_ALL lnl_red lnl_blue lnl_ALL ali_len

    outfile.write (str_out + lnl_out + str(len (align_ALL[0])) + "\n")

outfile.close()
os.system ("mv table_LnL.txt ../tableLnL_phyml_" + suffix + ".txt")
os.chdir  ("../")
os.system ("rm -rf tmp." + suffix)
