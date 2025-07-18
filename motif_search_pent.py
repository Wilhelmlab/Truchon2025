#!/usr/bin/env python3

#Filename: motif_search_pent.py
#Author: Alex Truchon
#Date: 2025-04-16
#Description: This script reads a bedmethyl file and outputs a csv table which contains the methylation frequencies of all possible pentamers
# which can surround an adenine

# Input the corrected bedmethyl table
file = input("Input the name of the corrected bedmethyl file: ")

# Read in the bedmethyl table
infile = open(file, 'r')
read = infile.readlines()

motifs = {}

# Read in the reference genome file
genome_in = open('../AaV_2021.fna')
genome_read = genome_in.readlines()
genome = ''

# Convert the genome file to a string
for line in genome_read:
	line = line.strip()
	if line.startswith('>'):
		continue
	else:
		genome += line

# Record the length of the genome
total = len(genome)
count = 0

# Open the bedmethyl file
for line in read:
	line = line.strip()
	coord = int(line.split('\t')[1])		# Define the location and the methylation frequency
	score = float(line.split('\t')[10])
	if score > -0.01:				# Use all scores
		count = count + 1
		if coord < 2:							# Identify the motif located within the first couple of base pairs
			ID = ('>seq%s\t%s\t%s' % (str(count), coord, score))
			start = total - (2 - coord)
			end = coord + 3
			motifs[ID] = genome[start:total] + genome[0:end]
		else:
			if coord > total - 2:						# Identify all remaining motifs throughout the genome
				ID = ('>seq%s\t%s\t%s' % (str(count), coord, score))
				start = coord - 2
				end = 3 - (total - coord)
				motifs[ID] = genome[start:total] + genome[0:end]
			else:
				ID = ('>seq%s\t%s\t%s' % (str(count), coord, score))
				start = coord - 2
				end = coord + 3
				motifs[ID] = genome[start:end]

# Generate dictionaries for scores and counts
score_dict = {}
count_dict = {}

# Create a translation table
fwd = 'ACGT'
rev = 'TGCA'
tab = str.maketrans(fwd, rev)

# For all the motifs that were identified, add the motif to the new dictionaries
for x in motifs:
	score = float(x.split('\t')[2])
	motif = motifs[x]
	if motif[2] == 'T':				# Consider the reverse motifs
		motif = motif.translate(tab)[::-1]
	if motif in score_dict:
		count_dict[motif] += 1			# Add scores and counts to existing motifs
		score_dict[motif] += score
	else:
		count_dict[motif] = 1			# Append the dictionary with motifs not yet included
		score_dict[motif] = score

ends = {}
starts = {}

# Open an outfile
outfile = open('pentamer-table.csv', 'w')

# Iterate through the score dictionary
for x in score_dict:
	end = x[2:5]		# Extract the first 3 nucleotides and last 3 nucleotides from each pentamer
	start = x[0:3]
	if end not in ends:
		ends[end] = ''			# Add to end and start dictionaries
		outfile.write(',%s' % (end))
	starts[start] = ''
del ends["AA"]

# Write each motif to a csv that has terminal nucleotides as columns and beginning nucleotides as rows
for x in starts:
	outfile.write('\n%s' % (x))
	for y in ends:
		motif = x[0:2] + y
		outfile.write(',%s' % (score_dict[motif]/count_dict[motif]))
