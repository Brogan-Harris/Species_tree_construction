# Imports and Libraries required #
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO
import pandas as pd
from collections import defaultdict


# Read in the species file
species_file = sys.argv[1]
species_dict = defaultdict(list)
with open(species_file) as f:
    for l in f:
        line = l.rstrip()
        species = line.replace(" ", "_")
        species_dict[species] = []

# Save all species.
species_set = set()
with open(species_file) as f:
        for l in f:
            species = l.rstrip()
            species_set.add(species)
print("Total number of species under analysis: ", len(species_set))


# Calulate the total number of orthogroups
orthogroups = []
total_orthogroups = 0
for orthogroup in glob.glob("*.fasta"):
    orthogroups.append(orthogroup)
    total_orthogroups += 1
print("Total Orthogroups:",total_orthogroups)

# Store failed Orthogroups
failed_orthogroups = []

# Make directory for succesfuly orthogroups
#os.mkdir("Successful_orthogroups")

# Loop through all the orthogroups
orthogroup_count = 0
for orthogroup in orthogroups:

    # Monitor progress
    orthogroup_count += 1
    print("Orthogroup under analysis: ", orthogroup)
    print("Orthologs analysed: ", orthogroup_count, "/", total_orthogroups)

    # Make a dictionary with all species/
    species_count = {}
    for species in species_set:
        species_count[species] = 0

    #pass in orthogroup sequence datadata
    inhandle = SeqIO.parse(orthogroup, "fasta")

    # Count the occurance of each sequence in both the ale and orignal groups.
    original_seq = []
    ale_seq = []

    ALE_orthogroup = {}
    species_seq = {}
    original_count = 0
    ale_count = 0

    original_desc = []
    new_desc = []


    # Loop through each record and convert it to an ALE record.
    for record in inhandle:
        original_count += 1
        record_description = record.description.replace(" ", "_")
        original_desc.append(record_description)
        seq = record.seq
        original_seq.append(seq)
        for species in species_set:
            if species in record_description:
                count = species_count[species]
                ale_count +=1
                species_count[species] += 1

                species_ale = species.replace("_", ".")
                ALE_reformat = species_ale + "_" + str(count+1)
                ALE_orthogroup[ALE_reformat] = seq
                ale_seq.append(seq)
                new_desc.append(ALE_reformat)

            else:
                continue
                
    # Back checks
    if len(original_seq) != len(ale_seq):
        print("Error missing sequences or species from :", orthogroup)
        failed_orthogroups.append(orthogroup)
    else:
        #Write new orthogroup
        file_name = "ALE_" + orthogroup
        ofile = open(file_name , "w")
        for key, val in ALE_orthogroup.items():
            ofile.write(">" + key + "\n" + str(val) + "\n")
        ofile.close()
    print("\n\n")

# Write failed orthogroups to file
ofile = open("failed_orthogroups.txt", "w")
for orthogroup in failed_orthogroups:
    ofile.write(orthogroup + "\n")
ofile.close()
print("Failed: ",len(failed_orthogroups),"/",total_orthogroups)

# Move all failed orthogroups to a new directory
os.mkdir("Failed_orthogroups")
for orthogroup in failed_orthogroups:
    my_command = "mv *" + orthogroup + "* Failed_orthogroups"
    output  = subprocess.check_output(my_command, shell=True,
                                      stderr=subprocess.STDOUT)
