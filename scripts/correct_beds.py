#!/usr/bin/env python3

#Filename: correct_beds.py
#Author: Alex Truchon
#Date: 2025-04-16
#Description: This script takes a bedmethyl file that is output from modkit and checks with the reference genome to only include "true" nucleotides that match.

# Read in the reference genome file in fasta format
infile1 = open('../AaV_2021.fna')
read1 = infile1.readlines()

# Convert the genome into a single string
genome = ''
for line in read1:
	if line.startswith('>'):
		continue
	else:
		line = line.strip()
		genome += line

# Ask for the bedmethyl file as input and the type of methylation identified in the bedmethyl file (either 6mA or 5mC)
file = input("Input the uncorrected bedmethyl file: ")
methyl = input("Type of methylation (A or C): ")

# Read in the input bedmethyl file
infile = open(file, 'r')
read = infile.readlines()
# Open an outfile for the corrected bedmethyl
outfname = file.replace('.bed', '.corrected.bed')
outfile = open(outfname, 'w')
for line in read:
	line = line.strip()
	position = int(line.split('\t')[1])	# Identify the genomic location of the modification
	strand = line.split('\t')[5]		# Identify whether the modification is on the plus or minus strand
	if strand == '+':
		if methyl == 'C':				# Confirm the corresponding position on the reference is an identical nucleotide
			if genome[position] == 'C':
				outfile.write('%s\n' % (line))
		else:						# Write only the correct nucleotides to the outfile
			if genome[position] == 'A':
				outfile.write('%s\n' % (line))
	if strand == '-':
		if methyl == 'C':				# Confirm the corresponding position on the reference is the opposite nucleotide
			if genome[position] == 'G':
				outfile.write('%s\n' % (line))
		else:						# Write only the correct nucleotides to the outfile
			if genome[position] == 'T':
				outfile.write('%s\n' % (line))
