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


filename = "E5_ALL_aa_reduced_mafftAln.fas"
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

# DEBUG: small alignments
#align_fix = align_fix[:,:500]
#align_var = align_var[:,:500]

filesuffix = "original"

outfile = open ("table-" + filesuffix + ".txt", "w", 0)

# we could also have used alignment.extend()
align_input = Align.MultipleSeqAlignment([i for i in align_var] + [i for i in align_fix1] + [i for i in align_fix2] + [i for i in align_fix3] + [i for i in align_fix4] + [i for i in align_fix5])

#   alpha1-alpha2
align_H1 = Align.MultipleSeqAlignment([i for i in align_fix1] + [i for i in align_fix2])

#   gammadelta-delta
align_H2 = Align.MultipleSeqAlignment([i for i in align_fix4] + [i for i in align_fix5])

#   alpha1-alpha2-gammadelta-delta
align_H3 = Align.MultipleSeqAlignment([i for i in align_fix1] + [i for i in align_fix2] + [i for i in align_fix4] + [i for i in align_fix5])

#   alpha1-alpha2-gammadelta-delta-epsilonzeta
align_H4 = Align.MultipleSeqAlignment([i for i in align_fix1] + [i for i in align_fix2] + [i for i in align_fix4] + [i for i in align_fix5] + [i for i in align_var])

#   alpha1-alpha2-epsilonzeta
align_H5 = Align.MultipleSeqAlignment([i for i in align_fix1] + [i for i in align_fix2] + [i for i in align_var])


outfile.write ("original alignment: " + phyml_ml (align_input, namesuffix) + str (len (align_input[0])) + "\n")
align_input = align (align_input)
outfile.write ("realignd alignment: " + phyml_ml (align_input, namesuffix) + str (len (align_input[0])) + "\n")

outfile.write ("original alpha1: " + phyml_ml (align_fix1, namesuffix) + str (len (align_fix1[0])) + "\n")
align_fix1 = align (align_fix1)
outfile.write ("realignd alpha1: " + phyml_ml (align_fix1, namesuffix) + str (len (align_fix1[0])) + "\n")

outfile.write ("original alpha2: " + phyml_ml (align_fix2, namesuffix) + str (len (align_fix2[0])) + "\n")
align_fix2 = align (align_fix2)
outfile.write ("realignd alpha2: " + phyml_ml (align_fix2, namesuffix) + str (len (align_fix2[0])) + "\n")

outfile.write ("original beta: " + phyml_ml (align_fix3, namesuffix) + str (len (align_fix3[0])) + "\n")
align_fix3 = align (align_fix3)
outfile.write ("realignd beta: " + phyml_ml (align_fix3, namesuffix) + str (len (align_fix3[0])) + "\n")

outfile.write ("original gammadelta: " + phyml_ml (align_fix4, namesuffix) + str (len (align_fix4[0])) + "\n")
align_fix4 = align (align_fix4)
outfile.write ("realignd gammadelta: " + phyml_ml (align_fix4, namesuffix) + str (len (align_fix4[0])) + "\n")

outfile.write ("original delta: " + phyml_ml (align_fix5, namesuffix) + str (len (align_fix5[0])) + "\n")
align_fix5 = align (align_fix5)
outfile.write ("realignd delta: " + phyml_ml (align_fix5, namesuffix) + str (len (align_fix5[0])) + "\n")

outfile.write ("original epsilonzeta: " + phyml_ml (align_var, namesuffix) + str (len (align_var[0])) + "\n")
align_var = align (align_var)
outfile.write ("realignd epsilonzeta: " + phyml_ml (align_var, namesuffix) + str (len (align_var[0])) + "\n")

outfile.write ("original alpha1-alpha2: " + phyml_ml (align_H1, namesuffix) + str (len (align_H1[0])) + "\n")
align_H1 = align (align_H1)
outfile.write ("realignd alpha1-alpha2: " + phyml_ml (align_H1, namesuffix) + str (len (align_H1[0])) + "\n")

outfile.write ("original gammadelta-delta: " + phyml_ml (align_H2, namesuffix) + str (len (align_H2[0])) + "\n")
align_H2 = align (align_H2)
outfile.write ("realignd gammadelta-delta: " + phyml_ml (align_H2, namesuffix) + str (len (align_H2[0])) + "\n")

outfile.write ("original alpha1-alpha2-gammadelta-delta: " + phyml_ml (align_H3, namesuffix) + str (len (align_H3[0])) + "\n")
align_H3 = align (align_H3)
outfile.write ("realignd alpha1-alpha2-gammadelta-delta: " + phyml_ml (align_H3, namesuffix) + str (len (align_H3[0])) + "\n")

outfile.write ("original alpha1-alpha2-gammadelta-delta-epsilonzeta: " + phyml_ml (align_H4, namesuffix) + str (len (align_H4[0])) + "\n")
align_H4 = align (align_H4)
outfile.write ("realignd alpha1-alpha2-gammadelta-delta-epsilonzeta: " + phyml_ml (align_H4, namesuffix) + str (len (align_H4[0])) + "\n")

outfile.write ("original alpha1-alpha2-epsilonzeta: " + phyml_ml (align_H5, namesuffix) + str (len (align_H5[0])) + "\n")
align_H5 = align (align_H5)
outfile.write ("realignd alpha1-alpha2-epsilonzeta: " + phyml_ml (align_H5, namesuffix) + str (len (align_H5[0])) + "\n")

outfile.close()