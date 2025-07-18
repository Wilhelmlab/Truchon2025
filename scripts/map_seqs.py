#!/usr/bin/env python3

#Filename: map_seqs.py
#Author: Alex Truchon
#Date: 2025-04-16
#Description: This script takes a set methylation frequency and calculates z-scores for the average methylation frequency for a
# 5000 bp region iterating through the genome by a 1000 bp sliding window

import statistics

# Read in the reference genome file in fasta format
infile1 = open('../AaV_2021.fna')
read1 = infile1.readlines()

# Convert the reference genome into a string
genome = ''
for line in read1:
	if line.startswith('>'):
		continue
	else:
		line = line.strip()
		genome += line

# Calculate the length of the genome
total = len(genome)

# Input the corrected bedmethyl file and an integer to use as reference for characterizing methylated sites
file = input("Input the corrected bedmethyl file: ")
lowest = float(input("Input the minimum score for assigning methylation: "))
type = input("Input the type of methylation (6mA or 5mC): ")

# Assign the correct nucleotides to the type of methylation
if type == '6mA':
	nuc = 'A'
	opp = 'T'
if type == '5mC':
	nuc = 'C'
	opp = 'G'

# Read in the corrected bedmethyl file
infile = open(file, 'r')
read = infile.readlines()

# Create a dictionary for highly methylated sites
hits_dict = {}
count = 0
for line in read:
	line = line.strip()
	position = int(line.split('\t')[1])		# Identify the position, scores, and strands of methylated sites
	strand = line.split('\t')[5]
	score = float(line.split('\t')[10])
	if score > lowest:
		count = count + 1					# Count the number of sites
		hits_dict[str(position) + '_' + strand] = score		# Add the described sites to a dictionary of methylated sites

avg = count/total

# Print statistics about respective genomic methylation
print ('The genome is %s kilobases long.\n' % (total/1000))
print ('There are %s methylated adenines every 5000 bases.\n' % (avg * 5000))

# Begin integers for interating through the reference genome
start = 0
end = 5000
fwd_total_count = 0
rev_total_count = 0
fwd_80_count = 0
rev_80_count = 0

# Create a new dictionary to describe the frequency of methylation throughout the genome
freq_dict = {}

# Generate a range that iterates 5000 base pairs every 1000 base pair
for i in range(start, total, 1000):
	fwd_range_count = 0
	rev_range_count = 0
	fwd_count = 0
	rev_count = 0
	if i + 5000 < total:			# Consider all regions before the terminal end of the genome
		reg = genome[i:i+5000]
	if i + 5000 > total:
		reg = genome[i:total] + genome[0:(5000-(total-i))]	# Create a 5000 bp region which transverses the end and start of the genome
	for n in reg:
		if n == nuc:					# Count the possible sites for methylation on the forward strand
			fwd_count = fwd_count + 1
			fwd_total_count = fwd_total_count + 1
		if n == opp:					# Count the possible sites for methylation on the reverse strand
			rev_count = rev_count + 1
			rev_total_count = rev_total_count + 1
	for x in hits_dict:
		pos = int(x.split('_')[0])			# Identify the sites from the previously established hit dictionary
		strand = x.split('_')[1]
		if i + 5000 < total:				# Determine whether the defined sites are located in the 5000 bp region
			if pos > i and pos < i + 5000:
				if strand == '+':
					fwd_range_count = fwd_range_count + 1	# Count the sites in the region on the forward strand
					fwd_80_count = fwd_80_count + 1
				else:
					rev_range_count = rev_range_count + 1	# Count the sites in the region on the reverse strand
					rev_80_count = rev_80_count + 1
		else:
			if pos > i or pos < 5000 - (total-i):			# Perform the same function for the final terminal region
				if strand == '+':				# of the reference genome
					fwd_range_count = fwd_range_count + 1
					fwd_80_count = fwd_80_count + 1
				else:
					rev_range_count = rev_range_count + 1
					rev_80_count = rev_80_count + 1
	fwdmethper = fwd_range_count/fwd_count
	revmethper = rev_range_count/rev_count			# Count the methylated sites per total sites in the range
	if i + 5000 < total:
		freq_dict[str(i) + '_' + str(i + 5000)] = str(fwdmethper) + '_' + str(revmethper)	# Add the region and the frequency of methylation
	else:												# to the frequency dictionary
		freq_dict[str(i) + '_' + str(5000 - (total-i))] = str(fwdmethper) + '_' + str(revmethper)

# Print total methylation statistics across the entire genome
print ('%s%% of the adenines on the forward strand are highly methylated. %s%% of the adenines on the negative strand are highly methylated.' % (str(int((fwd_80_count/fwd_total_count)*100000)/1000), str(int((rev_80_count/rev_total_count)*100000)/1000)))

# Create lists for summing the frequencies found in the frequency dictionary
count = 0
fwdsum = []
revsum = []

# Open the frequency dictionary
for x in freq_dict:
	scores = freq_dict[x]
	count = count + 1		# Identify the start, end, forward strand methylations
	start = x.split('_')[0]		# and reverse strand methylations
	end = x.split('_')[1]
	fwdmethper = float(scores.split('_')[0])
	revmethper = float(scores.split('_')[1])
	fwdsum.append(fwdmethper)
	revsum.append(revmethper)	# Append the methylation stats to each list

# Calculate mean and standard deviation for each set of methylation stats
fwdmean = statistics.mean(fwdsum)
fwdsd = statistics.stdev(fwdsum)
revmean = statistics.mean(revsum)
revsd = statistics.stdev(revsum)

# Create an outfile for the methylation map
outfname = type + "_methylation_map.txt"
outfile = open(outfname, 'w')
outfile.write('Start\tEnd\tFwd_freq\tRev_freq\tFwd_z-score\tRev_z-score\n')

# Open the frequency dictionary
for x in freq_dict:
	scores = freq_dict[x]
	start = x.split('_')[0]
	end = x.split('_')[1]
	fwdmethper = float(scores.split('_')[0])
	revmethper = float(scores.split('_')[1])
	fwd_z = (fwdmethper - fwdmean)/fwdsd		# Calculate Z-scores for forward and reverse regions using the mean and sd
	rev_z = (revmethper - revmean)/revsd
	outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (start, end, fwdmethper, revmethper, fwd_z, rev_z))	# Print the results to file
