#!/usr/bin/env python3

#Filename: motif_locs.py
#Author: Alex Truchon
#Date: 2025-04-16
#Description: This script takes a submitted motif and identifies its locations throughout a reference genome while collecting
# the methylation frequencies for each site

import re

# Input the corrected bedmethyl file, the motif you are searching for, and the position in the motif that is targeted for methylation
# For instance, if GATC was the motif in question for adenine methylation, p would equal 2
file = input("Input the corrected bedmethyl file: ")
m = input("Input the motif in question: ")
p = int(input("Input the location in the motif that is targeted for methylation: "))

# Create a translation dictionary for reverse motifs
# Note: Additional IUPAC codes may need to be added for more degenerate motifs
fwd = "CATGYN"
rev = "GTACRN"
transtab = str.maketrans(fwd, rev)

# Translate the motif into a reverse motif
rm = m[::-1]
rm = rm.translate(transtab)

# Read in the reference genome file and the bedmethyl file
infile = open(file, 'r')
read = infile.readlines()
gin = open('../AaV_2021.fna')
gread = gin.readlines()
genome = ''

# Create an outfile for motif scores
outfname = 'motif_scores_%s.csv' % (m)
outfile = open(outfname, 'w')

# Add the reference genome to a string
for line in gread:
	line = line.strip()
	if line.startswith('>'):
		continue
	else:
		genome += line

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

# Record the regex expression that is being searched across the genome
print ('Motif searched: %s' % (m))
print ('Reverse motif searched: %s' % (rm))

# Open the bedmethyl file
for line in read:
	line = line.strip()
	site = int(line.split('\t')[2])			# Identify genomic positions of the methylation sites
	motif = genome[site-p:site+mlen-p]
	motifrev = genome[site-mlen+p-1:site+p-1]	# Pull the correct region from the genome which would correspond to the motif in question
	score = float(line.split('\t')[10])		# Determine the score
	h = re.search(m, motif)
	if h != None:							# Print sites that correspond to the forward motif
		outfile.write('AaV_2021,%s,%s\n' % (site, score))
	else:
		h = re.search(rm, motifrev)				# Print sites that correspond to the reverse motif
		if h != None:
			outfile.write('AaV_2021,%s,%s\n' % (site, score))

