# Imports and Libraries required #
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO
import pandas as pd
#------------------------------------------------------------------------------#
# Read in list of species, and convert to a dictionary.
species_file = sys.argv[1]

# Number of files in directory
total_orthogroups = 0
for orthogroup in glob.glob("*.fa"):
    total_orthogroups += 1

# Save all species.
species_set = set()
with open(species_file) as f:
        for l in f:
            species = l.rstrip()
            species_set.add(species)
print("Total number of species under analysis: ", len(species_set))

# Make dictionary for species represtation
species_rep = {}

# Read in orthogroups, all orthogroups must be labelled ".fa"
orthogroup_count = 0

# Loop though all orthogroups in the folder.
for orthogroup in glob.glob("*.fa"):
    orthogroup_count += 1
    print("Orthogroup under analysis: ", orthogroup)
    print("Orthologs analysed: ", orthogroup_count, "/", total_orthogroups)
    inhandle = SeqIO.parse(orthogroup, "fasta")

    # Count all species in the orthogroup
    species_count = {}
    for species in species_set:
        species_count[species] = 0

    all_descriptions = {}
    for record in inhandle:
        record_description = record.description.replace(" ", "_")
        all_descriptions[record_description] = str(record.seq)
        #print(record_description)
        for species in species_set:
            if species in record_description:
                species_count[species] += 1

    # Identify all single copy and multi copy genes
    sc_orthgroup_species = set()
    mc_orthgroup_species = set()
    nc_orthgroup_species = set()
    removed_species = set()

    for key, val in species_count.items():
        if val == 1:
            sc_orthgroup_species.add(key)
        elif val > 1:
            mc_orthgroup_species.add(key)
            removed_species.add(key)
        else:
            nc_orthgroup_species.add(key)
            removed_species.add(key)

    # Write all single copy orthologs to a new fasta file
    ofile = open(orthogroup[:9] + "_sc.faa", "w")
    for description, seq in all_descriptions.items():
        for species in sc_orthgroup_species:
            if species in description:
                ofile.write(">" + description + "\n" + seq + "\n")
    ofile.close()

    # Not all species in each group and type or removal
    Total_species = len(species_set)
    Species_with_sc = len(sc_orthgroup_species)
    Species_with_mc = len(mc_orthgroup_species)
    Species_with_nc = len(nc_orthgroup_species)
    #print(len(removed_species))

    # Back checks
    sequences = 0
    sc_orthogroup = SeqIO.parse(orthogroup[:9] +"_sc.faa", "fasta")
    for sequence in sc_orthogroup:
        sequences += 1
    if sequences == Species_with_sc:
        print("sc orthofile correctly printed")
    else:
        print("Error with code")

    # Calculate species representation
    species_representation = (Species_with_sc/Total_species)*100
    o_file = orthogroup[:9] + "_sc.faa"
    #print("Orthogroup: ", o_file, "has ",species_representation, "species_rep")
    species_rep[o_file] = species_representation


# Move all orthogroups with >66% representation to a seperate folder
os.mkdir("Orthogroups_more_80pc_srep")
os.mkdir("Orthogroups_less_80pc_srep")

ofile = open("species_rep.csv", "w")
for key,val in species_rep.items():
    ofile.write(key + "," + str(val) + "\n")
    if val >= 80:
        my_command = "mv " + key + " Orthogroups_more_80pc_srep"
        output  = subprocess.check_output(my_command, shell=True,
                                      stderr=subprocess.STDOUT)
        print(key, " moved to >80 file")
    if val < 80:
        my_command = "mv " + key + " Orthogroups_less_80pc_srep"
        output  = subprocess.check_output(my_command, shell=True,
                                      stderr=subprocess.STDOUT)
        print(key, " moved to <80 file")

# Move original_orthogroups and results into a file
os.mkdir("Original_orthogroups")
os.mkdir("Results")

my_command = "mv *.fa Original_orthogroups"
output  = subprocess.check_output(my_command, shell=True,
                              stderr=subprocess.STDOUT)
my_command = "mv *.csv Results"
output  = subprocess.check_output(my_command, shell=True,
                              stderr=subprocess.STDOUT)

# Print summary
print("\n\n***** All files have been moved to correct folders ****")
print("\n**** Thanks for using PRem2 *****")
