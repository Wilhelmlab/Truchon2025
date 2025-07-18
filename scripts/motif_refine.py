#!/usr/bin/env python3

#Filename: motif_refine.py
#Author: Alex Truchon
#Date: 2025-04-16
#Description: This script reports genome wide statistics for a submitted motif

import re

# Import the corrected bedmethyl file, the motif of interest, and the location of methylation on the motif
file = input("Input the corrected bedmethyl file: ")
m = input("Input the motif of interest: ")
p = int(input("Input the position of methylation in the motif of interest: "))
seq = m

# Create a translation dictionary for reverse motifs
# Note: Additional IUPAC codes may need to be added for more degenerate motifs
fwd = "CATGYN"
rev = "GTACRN"
transtab = str.maketrans(fwd, rev)

# Translate the motif into a reverse motif
rm = m[::-1]
rm = rm.translate(transtab)

# Read in the reference genome file
infile = open(file, 'r')
read = infile.readlines()
gin = open('../AaV_2021.fna')
gread = gin.readlines()

# Convert the genome sequence into a string
genome = ''
for line in gread:
	line = line.strip()
	if line.startswith('>'):
		continue
	else:
		genome += line

# Create a list for all scores
scores = []

# Determine the length of the motif
mlen = len(m)

# Replace degenerate nucleotide codes with regex derivatives
# Note: For more degenerate motifs similar lines may need to be added
if 'N' in m:
	m = m.replace('N','[CTAG]')
	rm = rm.replace('N', '[CTAG]')
if 'Y' in m:
	m = m.replace('Y', '[CT]')
	rm = rm.replace('R', '[GA]')

# Read in the bedmethyl file and search for the motif in question within it
for line in read:
	line = line.strip()
	site = int(line.split('\t')[2])
	motif = genome[site-p:site+mlen-p]		# Pull the region surrounding the adenine
	motifrev = genome[site-mlen+p-1:site+p-1]
	score = float(line.split('\t')[10])		# Take the methylation frequency of the adenine
	h = re.search(m, motif)
	if h != None:
		scores.append(score)			# Append the matching forward and reverse motifs to the score list
	else:
		h = re.search(rm, motifrev)
		if h != None:
			scores.append(score)

# Print results
print ('Motif of interest: %s' % (seq))
print ('Total sites of motif %s: %s' % (seq, len(scores)))
total = 0
count = 0

for i in scores:
	total = total + i
	if i > 69.99:
		count = count + 1

print ('Total methylated sites of motif %s: %s' % (seq, count))
print ('Methylation frequency of motif %s sites: %s' % (seq, total/len(scores)))
