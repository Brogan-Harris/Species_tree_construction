# Script to constuct super matrix for phylogenetic analysis.
# Replaces all missing genes with graps
# Outputs Supermatrix alignment and warning for genes with high gaps.
# Written by Brogan Harris - 05/08/2019

# Imports and Libraries required #
import os
import sys
import glob
import time
import subprocess
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import matplotlib.pyplot as plt

# Print dependencies
print(""" Dependencies: Biopython, Matplotlib, Pandas, Seaborn
          Command line args: List of species
          Format:  All files suffix *.aln, all species "Arabidopsis_thaliana"
      """)

# Functions
def count_AA_occurence(sequence, amino_acids_to_track='ACDEFGHIKLMNPQRSTVWY-'):

    """ function to count amino acide occurences """

    amino_acids_per_train = 0
    amino_acids = dict.fromkeys(amino_acids_to_track, 0)

    for char in sequence:
        if char in amino_acids:
            amino_acids_per_train += 1
            amino_acids[char] += 1

    AA_percentages = {k:(v*100.0)/amino_acids_per_train for k,
                   v in amino_acids.items()}
    return AA_percentages

# Read in list of species, and convert to a dictionary.
species_file = sys.argv[1]
species_dict = defaultdict(list)
with open(species_file) as f:
    for l in f:
        line = l.rstrip()
        species = line.replace(" ", "_")
        species_dict[species] = []

# Store sequence and description
all_seqs = {}
all_descriptions = {}

# Read in the alignments
for alingment in glob.glob("*.aln"):
    inhandle = SeqIO.parse(alingment, "fasta")

    present_species = set()
    all_species = set()

    for record in inhandle:
        # Add tp dictionary / remove duplicates
        all_seqs[record.description] = str(record.seq)
        all_descriptions[record.description] = str(record.seq)

        # Reformat description to be the same as search key
        description = str(record.description)
        sequence = str(record.seq)
        seq_length = len(sequence)
        description_search = description.replace(" ", "_")

        # Search description to see if species if present
        for key, val in species_dict.items():
            key_search = key
            all_species.add(key_search)
            if key_search in description_search:
                species_dict[key_search].append(sequence)
                present_species.add(key_search)

    # Fill in gaps for all species with no sequences.
    not_present = all_species - present_species
    for key_search in not_present:
        species_dict[key_search].append("-" * seq_length)

# Concatenate all seqiences and write to file
Concatenated_seqs = {}
for key,val in species_dict.items():
    concat = "".join(val)
    Concatenated_seqs[key] = concat

# Name file super matrix + current time and write to fasta file
timestr = time.strftime("%d%m%y-%H%M")
file_name = "Concatenate_" + timestr + "_.fasta"
ofile = open(file_name, "w")
for key, val in Concatenated_seqs.items():
    ofile.write(">" + key + "\n" + val + "\n")
ofile.close()

# Convert concatenation to Panda database.
Species = []
Sequences = []
for key, val in Concatenated_seqs.items():
    Species.append(key)
    Sequences.append(val)
data = {'Species': Species, 'Sequences': Sequences}
df_seq = pd.DataFrame.from_dict(data)

# Loop through sequence and count amino acid
Species_AA_content = {}
for species, seq in Concatenated_seqs.items():
    AA_content = count_AA_occurence(seq)
    Species_AA_content[species] = AA_content

# Convert Amino acid dictionary to panda dataframe and orient same as df_seq
df_AA = pd.DataFrame.from_dict(Species_AA_content)
df_AA = df_AA.transpose()
df_AA = df_AA.reset_index()
df_AA.rename(columns={ df_AA.columns[0]: "Species" }, inplace = True)

# Merge the two databases by species
df_concat = df_seq.merge(df_AA, on = "Species")
df_concat.to_csv("df_concat.txt")
print(df_concat)

# Make graphs of amino acid content and gap
sns.set(rc={'figure.figsize':(10,20)})
Gap_content = sns.barplot(x="Species", y="-", data=df_concat)
Gap_content.set_title("Concatenation Gap Content")
Gap_content.set_xlabel("Species")
Gap_content.set_ylabel("% of gaps in alignment")
figure = Gap_content.get_figure()
figure.savefig('Concatenation_gap_content.png', dpi=400)

# Make graph of problem species
df_problem_species = df_concat[(df_concat["-"] >= 80)]
if len(df_problem_species) >= 1:
    ax = sns.barplot(x="Species", y="-", data=df_problem_species)
    ax.set_title("High Gap Content species")
    ax.set_xlabel("Species")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
    ax.set_ylabel("% of gaps in alignment")
    ax.tick_params(labelsize=10)
    figure = ax.get_figure()
    figure.savefig('High_gap_content.png', dpi=400)


# Make new df without the problem species and convert to list.
df_conservative = df_concat[df_concat["-"] <= 80]
print(df_conservative)
conservative_species = list(df_conservative.Species)

# Write conservative species concatenation to a fasta file.
file_name = "Conservative_concatenate_" + timestr + "_.fasta"
ofile = open(file_name, "w")
for species in conservative_species:
    ofile.write(">" + species + "\n" + Concatenated_seqs[species] + "\n")
ofile.close()

# Write new species list
file_name = "New_species_list_" + timestr + ".txt"
ofile = open(file_name, "w")
for species in conservative_species:
    ofile.write(species + "\n")
ofile.close()

# Write list of species removed form the alignment
problem_species = list(df_problem_species.Species)
file_name = "Species_removed_from_concat_" + timestr + ".txt"
ofile = open(file_name, "w")
for species in problem_species:
    ofile.write(species + "\n")
ofile.close()

# Sort files
os.mkdir("Results")
os.mkdir("Concatenations")
os.mkdir("Original_trimmed_orthogroups")
os.mkdir("Super_matrix_py_script")

my_commands = ["mv *.txt Results",
               "mv *.png Results",
               "mv *.fasta Concatenations",
               "mv *.aln Original_trimmed_orthogroups",
               "mv *.py Super_matrix_py_script"]

for my_command in my_commands:
    output  = subprocess.check_output(my_command, shell=True,
                                  stderr=subprocess.STDOUT)

# Write summary text file
file_name = "Summary_" + timestr + ".txt"
ofile = open(file_name, "w")
ofile.write("Summary of Supermatrix.py output" + "\n" + \
            "Length of alignment = " + str(seq_length) + "\n" + \
            "All results, concatenations and dfs have been sorted" + "\n")
ofile.close()
