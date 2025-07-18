#!/usr/bin/env python3

#Filename: nucl_freq.py
#Author: Alex Truchon
#Date: 2025-04-16
#Description: This script uses a minimum methylation frequency as a reference to identify the abundance of nucltides in
# specific positions upstream and downstream surrounding a methylated adenine

# Input the corrected bedmethyl file and an integer to use as reference for characterizing methylated sites
file = input("Input the corrected bedmethyl file: ")
lowest = float(input("Input the minimum score for assigning methylation: "))
p = int(input("Input the position in question surrounding methylated sites: "))

# Read in the bedmethyl file
infile = open(file, 'r')
read = infile.readlines()

# Read in the reference genome file
genome_in = open('../AaV_2021.fna')
genome_read = genome_in.readlines()
genome = ''

# Convert the reference genome file to a string
for line in genome_read:
	line = line.strip()
	if line.startswith('>'):
		continue
	else:
		genome += line

# Calculate the length of the genome
total = len(genome)

# Generate dictionaries for the high scoring methylation sites
dict80minus = {}
dict80plus = {}
dictall = {}
for line in read:
	line = line.strip()
	coord = int(line.split('\t')[1])	# Locate the coordinates and scores of the respective sites
	score = float(line.split('\t')[10])
	strand = line.split('\t')[5]
	dictall[coord] = str(score) + '_' + strand	# Add all to the collective dictionary
	if score > lowest:
		if strand == '-':
			dict80minus[coord] = score	# Add high scoring methylation sites to respective dicitonary
		else:
			dict80plus[coord] = score

# Load the position in question
pos = p

Acountu = 0
Ccountu = 0
Tcountu = 0
Gcountu = 0
Acountd = 0
Ccountd = 0
Tcountd = 0
Gcountd = 0

# Count the number of each nucleotide surrounding high scoring positions
for x in dict80plus:
	up = x - pos
	down = x + pos
	nuc_up = genome[up]
	nuc_down = genome[down]
	if nuc_up == 'A':
		Acountu = Acountu + 1
	if nuc_up == 'C':
		Ccountu = Ccountu + 1
	if nuc_up == 'T':
		Tcountu = Tcountu + 1
	if nuc_up == 'G':
		Gcountu = Gcountu + 1
	if nuc_down == 'A':
		Acountd = Acountd + 1
	if nuc_down == 'C':
		Ccountd = Ccountd + 1
	if nuc_down == 'G':
		Gcountd = Gcountd + 1
	if nuc_down == 'T':
		Tcountd = Tcountd + 1

for x in dict80minus:
	up = x + pos
	down = x - pos
	nuc_up = genome[up]
	nuc_down = genome[down]
	if nuc_up == 'A':
		Tcountu = Tcountu + 1
	if nuc_up == 'C':
		Gcountu = Gcountu + 1
	if nuc_up == 'T':
		Acountu = Acountu + 1
	if nuc_up == 'G':
		Ccountu = Ccountu + 1
	if nuc_down == 'A':
		Tcountd = Tcountd + 1
	if nuc_down == 'C':
		Gcountd = Gcountd + 1
	if nuc_down == 'T':
		Acountd = Acountd + 1
	if nuc_down == 'G':
		Ccountd = Ccountd + 1

print (pos, 'Upstream')
print ('A', Acountu)
print ('T', Tcountu)
print ('C', Ccountu)
print ('G', Gcountu, '\n')
print (pos, 'Downstream')
print ('A', Acountd)
print ('T', Tcountd)
print ('C', Ccountd)
print ('G', Gcountd, '\n')

print ('Percent adenine at positive', pos, ':', Acountd/(Acountd + Tcountd + Ccountd + Gcountd))
print ('Percent thymine at positive', pos, ':', Tcountd/(Acountd + Tcountd + Ccountd + Gcountd))
print ('Percent cytosine at positive', pos, ':', Ccountd/(Acountd + Tcountd + Ccountd + Gcountd))
print ('Percent guanine at positive', pos, ':', Gcountd/(Acountd + Tcountd + Ccountd + Gcountd), '\n')
print ('Percent adenine at negative', pos, ':', Acountu/(Acountu + Tcountu + Ccountu + Gcountu))
print ('Percent thymine at negative', pos, ':', Tcountu/(Acountu + Tcountu + Ccountu + Gcountu))
print ('Percent cytosine at negative', pos, ':', Ccountu/(Acountu + Tcountu + Ccountu + Gcountu))
print ('Percent guanine at negative', pos, ':', Gcountu/(Acountu + Tcountu + Ccountu + Gcountu))
