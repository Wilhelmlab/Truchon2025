#!/usr/bin/env python3

#Filename: motif_recip.py
#Author: Alex Truchon
#Date: 2025-04-16
#Description: This script takes as input a reciprocal motif and calculates statistics for complete, hemi, and unmethylated sites in the genome

import re

# Import the corrected bedmethyl file, the motif of interest, and the location of methylation on the motif
file = input("Input the corrected bedmethyl file: ")
m = input("Input the motif of interest: ")
p = int(input("Input the position of methylation in the motif of interest: "))

# Read in the corrected bedmethyl file
infile = open(file, 'r')
read = infile.readlines()
gin = open('../AaV_2021.fna')
gread = gin.readlines()
genome = ''

# Convert the genome to a string
for line in gread:
	line = line.strip()
	if line.startswith('>'):
		continue
	else:
		genome += line

# Make a dictionary for methylation frequencies
score_dict = {}
for line in read:
	line = line.strip()
	site = int(line.split('\t')[2])
	score = float(line.split('\t')[10])
	score_dict[site] = score

# Create count integers for full, hemi, and unmethylated sites
all = 0
non = 0
hemi = 0
full = 0

# Calculate numbers necessary for defining position
mlen = len(m)
rp = mlen - p + 1
diff = p - rp

# Replace degenerate nucleotides with the correct regular expression
if 'N' in m:
	m = m.replace('N', '[CTAG]')

# Bring up each site in the score dictionary
for x in score_dict:
	site = x
	score = score_dict[site]			# Search for each motif in the genome to check if it corresponds to the site
	motif = genome[site-p:site+mlen-p]
	h = re.search(m, motif)
	if h != None:
		all = all + 1
		rscore = score_dict[site - diff]
		if score >= 70:				# Count if the score and the reverse score are fully methylated
			if rscore >= 70:
				full = full + 1
			else:
				hemi = hemi + 1
		else:
			if rscore >= 70:
				hemi = hemi + 1
			else:
				non = non + 1

total = non + hemi + full
if all == total:
	print ("All motifs are accounted for.")
else:
	print ("Motifs are missing based on the regex search")

# Print results
print ("Proportion fully methylated: %s" % (full/total))
print ("Proportion hemi-methylated: %s" % (hemi/total))
print ("Proportion unmethylated: %s" % (non/total))
