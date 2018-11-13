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
  SeqIO.write(alignment, "align" + namesuffix + ".phy", "phylip") # alignment in phylip format 
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
#align_fix = align_fix[:,:500]
#align_var = align_var[:,:500]

filesuffix = "original"

outfile = open ("table-" + filesuffix + ".txt", "w", 0)

# we could also have used alignment.extend()
align_input = Align.MultipleSeqAlignment([i for i in align_var] + [i for i in align_fix1] + [i for i in align_fix2] + [i for i in align_fix3] + [i for i in align_fix4] + [i for i in align_fix5])

#   blue-yellow-red-black-green
align_H1 = Align.MultipleSeqAlignment([i for i in align_var] + [i for i in align_fix2] + [i for i in align_fix3] + [i for i in align_fix4] + [i for i in align_fix5])

#   yellow-red-black-green
align_H2 = Align.MultipleSeqAlignment([i for i in align_var] + [i for i in align_fix2] + [i for i in align_fix3] + [i for i in align_fix4])

#   red-black-green
align_H3 = Align.MultipleSeqAlignment([i for i in align_var] + [i for i in align_fix2] + [i for i in align_fix4])

#   black-green
align_H4 = Align.MultipleSeqAlignment([i for i in align_fix2] + [i for i in align_fix4])

#   red-black
align_H5_1 = Align.MultipleSeqAlignment([i for i in align_var] + [i for i in align_fix4])

#   green-yellow
align_H5_2 = Align.MultipleSeqAlignment([i for i in align_fix2] + [i for i in align_fix3])

outfile.write ("original alignment: " + phyml_ml (align_input, namesuffix) + str (len (align_input[0])) + "\n")
align_input = align (align_input)
outfile.write ("realignd alignment: " + phyml_ml (align_input, namesuffix) + str (len (align_input[0])) + "\n")

outfile.write ("original AvesTurtles: " + phyml_ml (align_fix1, namesuffix) + str (len (align_fix1[0])) + "\n")
align_fix1 = align (align_fix1)
outfile.write ("realignd AvesTurtles: " + phyml_ml (align_fix1, namesuffix) + str (len (align_fix1[0])) + "\n")

outfile.write ("original BetaGamma: " + phyml_ml (align_fix2, namesuffix) + str (len (align_fix2[0])) + "\n")
align_fix2 = align (align_fix2)
outfile.write ("realignd BetaGamma: " + phyml_ml (align_fix2, namesuffix) + str (len (align_fix2[0])) + "\n")

outfile.write ("original Lambda: " + phyml_ml (align_fix3, namesuffix) + str (len (align_fix3[0])) + "\n")
align_fix3 = align (align_fix3)
outfile.write ("realignd Lambda: " + phyml_ml (align_fix3, namesuffix) + str (len (align_fix3[0])) + "\n")

outfile.write ("original Manatees: " + phyml_ml (align_fix4, namesuffix) + str (len (align_fix4[0])) + "\n")
align_fix4 = align (align_fix4)
outfile.write ("realignd Manatees: " + phyml_ml (align_fix4, namesuffix) + str (len (align_fix4[0])) + "\n")

outfile.write ("original PerissoArtio: " + phyml_ml (align_fix5, namesuffix) + str (len (align_fix5[0])) + "\n")
align_fix5 = align (align_fix5)
outfile.write ("realignd PerissoArtio: " + phyml_ml (align_fix5, namesuffix) + str (len (align_fix5[0])) + "\n")

outfile.write ("original Alpha: " + phyml_ml (align_var, namesuffix) + str (len (align_var[0])) + "\n")
align_var = align (align_var)
outfile.write ("realignd Alpha: " + phyml_ml (align_var, namesuffix) + str (len (align_var[0])) + "\n")

outfile.write ("original H1: " + phyml_ml (align_H1, namesuffix) + str (len (align_H1[0])) + "\n")
align_H1 = align (align_H1)
outfile.write ("realignd H1: " + phyml_ml (align_H1, namesuffix) + str (len (align_H1[0])) + "\n")

outfile.write ("original H2: " + phyml_ml (align_H2, namesuffix) + str (len (align_H2[0])) + "\n")
align_H2 = align (align_H2)
outfile.write ("realignd H2: " + phyml_ml (align_H2, namesuffix) + str (len (align_H2[0])) + "\n")

outfile.write ("original H3: " + phyml_ml (align_H3, namesuffix) + str (len (align_H3[0])) + "\n")
align_H3 = align (align_H3)
outfile.write ("realignd H3: " + phyml_ml (align_H3, namesuffix) + str (len (align_H3[0])) + "\n")

outfile.write ("original H4: " + phyml_ml (align_H4, namesuffix) + str (len (align_H4[0])) + "\n")
align_H4 = align (align_H4)
outfile.write ("realignd H4: " + phyml_ml (align_H4, namesuffix) + str (len (align_H4[0])) + "\n")

outfile.write ("original H5_1: " + phyml_ml (align_H5_1, namesuffix) + str (len (align_H5_1[0])) + "\n")
align_H5_1 = align (align_H5_1)
outfile.write ("realignd H5_1: " + phyml_ml (align_H5_1, namesuffix) + str (len (align_H5_1[0])) + "\n")

outfile.write ("original H5_2: " + phyml_ml (align_H5_2, namesuffix) + str (len (align_H5_2[0])) + "\n")
align_H5_2 = align (align_H5_2)
outfile.write ("realignd H5_2: " + phyml_ml (align_H5_2, namesuffix) + str (len (align_H5_2[0])) + "\n")

outfile.close()
