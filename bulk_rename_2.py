# Scrip to remove specific species from alignments.
# Written by Brogan Harris - 15/08/2019

# Imports and Libraries required #
import os
import sys
import glob
from Bio import SeqIO
from Bio.Seq import Seq


count = 0 

# Read in the alignments
for alingment in glob.glob("*fa*"):
    new_file = str(count) + ".fasta"
    os.rename(alingment, new_file)
    count += 1
    

