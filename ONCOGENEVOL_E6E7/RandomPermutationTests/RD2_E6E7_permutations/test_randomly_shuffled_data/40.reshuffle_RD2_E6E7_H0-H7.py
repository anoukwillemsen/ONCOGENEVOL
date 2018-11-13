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

filename = "original_ALL.phy"
align_ALL = AlignIO.read (filename, "phylip-relaxed")
seqnames = [seq.id for seq in align_ALL]

idx = [seqnames.index(i) for i in ["CcPV1", "CmPV1", "FcPV1", "FgPV1", "FlPV1", "PaPV1", "PePV1"]]
align_fix1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # AvesTurtles
idx = [seqnames.index(i) for i in ["ChPV1", "EePV1", "FcaPV3", "HPV123", "MfPV2", "VvPV1"]]
align_fix2 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # BetaGamma
idx = [seqnames.index(i) for i in ["CPV6", "FcaPV1", "MscPV2", "OcPV1", "RfPV1", "TePV1"]]
align_fix3 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # Lambda
idx = [seqnames.index(i) for i in ["TmPV2", "TmPV3", "TmPV1", "TmPV4"]]
align_fix4 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # Manatees
idx = [seqnames.index(i) for i in ["CdPV1", "CePV2", "EcPV1", "EcPV2", "OaPV1", "OvPV1"]]
align_fix5 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # PerissoArtio
idx = [seqnames.index(i) for i in ["HPV42", "HPV59", "HPV78", "SscPV3", "SsPV1", "UmPV1"]]
align_var = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  # Alpha

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
    
    idx = [seqnames.index(i) for i in ["HPV42", "HPV59", "HPV78", "SscPV3", "SsPV1", "UmPV1"]]
    align_Alpha = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_Alpha, "mc.paml")
    strout_Alpha, gamma, lnl_Alpha = phyml_ml ("original_Alpha.phy_phyml_tree")
    str_out += strout_Alpha
    model_out += gamma
    lnl_out += lnl_Alpha + " "

    idx = [seqnames.index(i) for i in ["CcPV1", "CmPV1", "FcPV1", "FgPV1", "FlPV1", "PaPV1", "PePV1"]]
    align_AvesTurtles = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_AvesTurtles, "mc.paml")
    strout_AvesTurtles, gamma, lnl_AvesTurtles = phyml_ml ("original_AvesTurtles.phy_phyml_tree")
    str_out += strout_AvesTurtles
    model_out += gamma
    lnl_out += lnl_AvesTurtles + " "

    idx = [seqnames.index(i) for i in ["ChPV1", "EePV1", "FcaPV3", "HPV123", "MfPV2", "VvPV1"]]
    align_BetaGamma = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_BetaGamma, "mc.paml")
    strout_BetaGamma, gamma, lnl_BetaGamma = phyml_ml ("original_BetaGamma.phy_phyml_tree")
    str_out += strout_BetaGamma
    model_out += gamma 
    lnl_out += lnl_BetaGamma + " "

    idx = [seqnames.index(i) for i in ["CPV6", "FcaPV1", "MscPV2", "OcPV1", "RfPV1", "TePV1"]]
    align_Lambda = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_Lambda, "mc.paml")
    strout_Lambda, gamma, lnl_Lambda = phyml_ml ("original_Lambda.phy_phyml_tree")
    str_out += strout_Lambda
    model_out += gamma
    lnl_out += lnl_Lambda + " "

    idx = [seqnames.index(i) for i in ["TmPV2", "TmPV3", "TmPV1", "TmPV4"]]
    align_Manatees = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_Manatees, "mc.paml")
    strout_all4Manatees, gamma, lnl_all4Manatees = phyml_ml ("original_Manatees.phy_phyml_tree")
    str_out += strout_all4Manatees
    model_out += gamma
    lnl_out += lnl_all4Manatees + " "

    idx = [seqnames.index(i) for i in ["CdPV1", "CePV2", "EcPV1", "EcPV2", "OaPV1", "OvPV1"]]
    align_PerissoArtio = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_PerissoArtio, "mc.paml")
    strout_PerissoArtio, gamma, lnl_PerissoArtio = phyml_ml ("original_PerissoArtio.phy_phyml_tree")
    str_out += strout_PerissoArtio
    model_out += gamma
    lnl_out += lnl_PerissoArtio + " "

######################################
# alternative IO hypotheses
# blue-yellow-red-black-green
    idx = [seqnames.index(i) for i in ["CdPV1", "CePV2", "EcPV1", "EcPV2", "OaPV1", "OvPV1", "CPV6", "FcaPV1", "MscPV2", "OcPV1", "RfPV1", "TePV1", "HPV42", "HPV59", "HPV78", "SscPV3", "SsPV1", "UmPV1", "TmPV2", "TmPV3", "TmPV1", "TmPV4", "ChPV1", "EePV1", "FcaPV3", "HPV123", "MfPV2", "VvPV1"]]
    align_H1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx]) 
    write_to_phylip (align_H1, "mc.paml")
    strout_H1, gamma, lnl_H1 = phyml_ml ("original_H1.phy_phyml_tree")
    str_out += strout_H1
    model_out += gamma
    lnl_out += lnl_H1 + " "

# yellow-red-black-green
    idx = [seqnames.index(i) for i in ["CPV6", "FcaPV1", "MscPV2", "OcPV1", "RfPV1", "TePV1", "HPV42", "HPV59", "HPV78", "SscPV3", "SsPV1", "UmPV1", "TmPV2", "TmPV3", "TmPV1", "TmPV4", "ChPV1", "EePV1", "FcaPV3", "HPV123", "MfPV2", "VvPV1"]]
    align_H2 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_H2, "mc.paml")
    strout_H2, gamma, lnl_H2 = phyml_ml ("original_H2.phy_phyml_tree")
    str_out += strout_H2
    model_out += gamma
    lnl_out += lnl_H2 + " "

# red-black-green
    idx = [seqnames.index(i) for i in ["HPV42", "HPV59", "HPV78", "SscPV3", "SsPV1", "UmPV1", "TmPV2", "TmPV3", "TmPV1", "TmPV4", "ChPV1", "EePV1", "FcaPV3", "HPV123", "MfPV2", "VvPV1"]]
    align_H3 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_H3, "mc.paml")
    strout_H3, gamma, lnl_H3 = phyml_ml ("original_H3.phy_phyml_tree")
    str_out += strout_H3
    model_out += gamma
    lnl_out += lnl_H3 + " "

# black-green
    idx = [seqnames.index(i) for i in ["TmPV2", "TmPV3", "TmPV1", "TmPV4", "ChPV1", "EePV1", "FcaPV3", "HPV123", "MfPV2", "VvPV1"]]
    align_H4 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_H4, "mc.paml")
    strout_H4, gamma, lnl_H4 = phyml_ml ("original_H4.phy_phyml_tree")
    str_out += strout_H4
    model_out += gamma
    lnl_out += lnl_H4 + " "

# red-black
    idx = [seqnames.index(i) for i in ["HPV42", "HPV59", "HPV78", "SscPV3", "SsPV1", "UmPV1", "TmPV2", "TmPV3", "TmPV1", "TmPV4"]]
    align_H5_1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_H5_1, "mc.paml")
    strout_H5_1, gamma, lnl_H5_1 = phyml_ml ("original_H5_1.phy_phyml_tree")
    str_out += strout_H5_1
    model_out += gamma
    lnl_out += lnl_H5_1 + " "

# green-yellow
    idx = [seqnames.index(i) for i in ["ChPV1", "EePV1", "FcaPV3", "HPV123", "MfPV2", "VvPV1", "CPV6", "FcaPV1", "MscPV2", "OcPV1", "RfPV1", "TePV1"]]
    align_H5_2 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])
    write_to_phylip (align_H5_2, "mc.paml")
    strout_H5_2, gamma, lnl_H5_2 = phyml_ml ("original_H5_2.phy_phyml_tree")
    str_out += strout_H5_2
    model_out += gamma
    lnl_out += lnl_H5_2 + " "
######################################
# ALL together
    write_to_phylip (align_ALL, "mc.paml")
    strout_ALL, gamma, lnl_ALL = phyml_ml ("original_ALL.phy_phyml_tree")
    str_out += strout_ALL
    model_out += gamma
    lnl_out += lnl_ALL + " "   
######################################

    # output: ML tree lengths
    # outorder: treelen_red	treelen_gray	treelen_green	treelen_yellow	treelen_black	treelen_blue	treelen_part_H1	treelen_part_H2	treelen_part_H3	treelen_part_H4	treelen_part_H5_1	treelen_part_H5_2	treelen_ALL	lnl_red	lnl_gray	lnl_green	lnl_yellow	lnl_black	lnl_blue	lnl_part_H1	lnl_part_H2	lnl_part_H3	lnl_part_H4	lnl_part_H5_1	lnl_part_H5_2	lnl_ALL	ali_len

    outfile.write (str_out + lnl_out + str(len (align_ALL[0])) + "\n")

outfile.close()
os.system ("mv table_LnL.txt ../tableLnL_" + suffix + ".txt")
os.chdir  ("../")
os.system ("rm -rf tmp." + suffix)
