#!/usr/bin/env python2.7
""" Read 4 seqs alignments, and calculates likelihoods without reshuffling """

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

namesuffix = ""

def phyml_ml (alignment, suffix):
  """ calculates ML tree and branch lengths  """
  namesuffix = str(random.randint(1,1e10)) 
  SeqIO.write(alignment, "align" + namesuffix + ".phy", "phylip-relaxed") # alignment in phylip format 
  os.system ("phyml -i align" + namesuffix + ".phy -d aa -m LG -o tlr -b 0 -c 1") # run phyML
  tree = Phylo.read("align" + namesuffix + ".phy_phyml_tree", "newick") # read ML tree
  strout = str(tree.total_branch_length()) + " "            # store brlen to print as output
  for line in open("align" + namesuffix + ".phy_phyml_stats"):          # search for lnL value
    if ". Log-likelihood:" in line:
      strout += str(line.split(" ")[-1].strip()) + " "           # store last column  
      break

  return strout

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


filename = "E5_red-blue-yellow_clade_aa_mafftAln.fas"
align_ALL = AlignIO.read (filename, "fasta")
seqnames = [seq.id for seq in align_ALL]

idx = [seqnames.index(i) for i in ["MUC_HPV16_E5_ALPHA", "MUC_HPV18_E5_ALPHA", "MUC_HPV31_E5_ALPHA", "MUC_HPV33_E5_ALPHA", "MUC_HPV34_E5_ALPHA", "MUC_HPV35_E5_ALPHA", "MUC_HPV39_E5_ALPHA", "MUC_HPV45_E5_ALPHA", "MUC_HPV52_E5_ALPHA", "MUC_HPV58_E5_ALPHA", "MUC_HPV59_E5_ALPHA", "MUC_HPV67_E5_ALPHA", "MUC_HPV68_E5_ALPHA", "MUC_HPV70_E5_ALPHA", "MUC_HPV73_E5_ALPHA", "MUC_HPV85_E5_ALPHA", "MUC_HPV97_E5_ALPHA", "CUT_CgPV1_E5_BETA", "CUT_HPV2_E5_BETA", "CUT_HPV3_E5_BETA", "CUT_HPV10_E5_BETA", "CUT_HPV27_E5_BETA", "CUT_HPV28_E5_BETA", "CUT_HPV29_E5_BETA", "CUT_HPV57_E5_BETA", "CUT_HPV61_E5_BETA", "CUT_HPV62_E5_BETA", "CUT_HPV71_E5_BETA", "CUT_HPV72_E5_BETA", "CUT_HPV77_E5_BETA", "CUT_HPV78_E5_BETA", "CUT_HPV81_E5_BETA", "CUT_HPV83_E5_BETA", "CUT_HPV84_E5_BETA", "CUT_HPV86_E5_BETA", "CUT_HPV87_E5_BETA", "CUT_HPV89_E5_BETA", "CUT_HPV90_E5_BETA", "CUT_HPV94_E5_BETA", "CUT_HPV102_E5_BETA", "CUT_HPV106_E5_BETA", "CUT_HPV114_E5_BETA", "CUT_HPV117_E5_BETA", "CUT_HPV125_E5_BETA", "CUT_HPV160_E5_BETA", "MUC_MfPV3_E5_EPSILON_DSETA", "MUC_MfPV4_E5_EPSILON_DSETA", "MUC_MfPV5_E5_EPSILON_DSETA", "MUC_MfPV6_E5_EPSILON_DSETA", "MUC_MfPV7_E5_EPSILON_DSETA", "MUC_MfPV8_E5_EPSILON_DSETA", "MUC_MfPV9_E5_EPSILON_DSETA", "MUC_MfPV10_E5_EPSILON_DSETA", "MUC_MfPV11_E5_EPSILON_DSETA", "MUC_MfuPV1_E5_EPSILON_DSETA", "MUC_MmPV1_E5_EPSILON_DSETA", "MUC_PhPV1_E5_EPSILON_DSETA", "GW_HPV6_E5_GAMMA_DELTA", "GW_HPV7_E5_DELTA", "GW_HPV11_E5_GAMMA_DELTA", "GW_HPV13_E5_GAMMA_DELTA", "GW_HPV40_E5_DELTA", "GW_HPV43_E5_DELTA", "GW_HPV44_E5_GAMMA_DELTA", "GW_HPV54_E5", "GW_HPV74_E5_GAMMA_DELTA", "GW_HPV91_E5_DELTA", "GW_PpPV1_E5_GAMMA_DELTA", "GW_CchPV1_E5_GAMMA_DELTA"]]
align_fix1 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx]) #red
#idx = [seqnames.index(i) for i in ["SfPV1", "OcPV1"]]
#align_fix2 = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  #yellow
idx = [seqnames.index(i) for i in ["AaPV1", "BgPV1", "BPV13", "BPV14", "BPV1", "BPV2", "CcaPV1", "OaPV1", "OaPV2", "OvPV1", "RalPV1", "RtPV1"]]
align_var = Align.MultipleSeqAlignment([align_ALL[i] for i in idx])  #blue

# DEBUG: small alignments
#align_fix = align_fix[:,:500]
#align_var = align_var[:,:500]

filesuffix = "original"

outfile = open ("table-" + filesuffix + ".txt", "w", 0)

# we could also have used alignment.extend()
align_input = Align.MultipleSeqAlignment([i for i in align_var] + [i for i in align_fix1])

##   red-blue
#align_H1 = Align.MultipleSeqAlignment([i for i in align_fix1] + [i for i in align_var])

##   blue-yellow
#align_H2 = Align.MultipleSeqAlignment([i for i in align_var] + [i for i in align_fix2])


outfile.write ("original alignment: " + phyml_ml (align_input, namesuffix) + str (len (align_input[0])) + "\n")
align_input = align (align_input)
outfile.write ("realignd alignment: " + phyml_ml (align_input, namesuffix) + str (len (align_input[0])) + "\n")

outfile.write ("original red: " + phyml_ml (align_fix1, namesuffix) + str (len (align_fix1[0])) + "\n")
align_fix1 = align (align_fix1)
outfile.write ("realignd red: " + phyml_ml (align_fix1, namesuffix) + str (len (align_fix1[0])) + "\n")

#outfile.write ("original yellow: " + phyml_ml (align_fix2, namesuffix) + str (len (align_fix2[0])) + "\n")
#align_fix2 = align (align_fix2)
#outfile.write ("realignd yellow: " + phyml_ml (align_fix2, namesuffix) + str (len (align_fix2[0])) + "\n")

outfile.write ("original blue: " + phyml_ml (align_var, namesuffix) + str (len (align_var[0])) + "\n")
align_var = align (align_var)
outfile.write ("realignd blue: " + phyml_ml (align_var, namesuffix) + str (len (align_var[0])) + "\n")

#outfile.write ("original red-blue: " + phyml_ml (align_H1, namesuffix) + str (len (align_H1[0])) + "\n")
#align_H1 = align (align_H1)
#outfile.write ("realignd red-blue: " + phyml_ml (align_H1, namesuffix) + str (len (align_H1[0])) + "\n")

#outfile.write ("original blue-yellow: " + phyml_ml (align_H2, namesuffix) + str (len (align_H2[0])) + "\n")
#align_H2 = align (align_H2)
#outfile.write ("realignd blue-yellow: " + phyml_ml (align_H2, namesuffix) + str (len (align_H2[0])) + "\n")


outfile.close()
