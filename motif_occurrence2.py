#!/usr/bin/env python3

#Filename: motif_occurrence.py
#Author: Alex Truchon
#Date: 2025-04-16
#Description: This script takes a minimum methylation frequency and finds the average methylation frequency of all adenines surrounding the methylated adenines

# Input the corrected bedmethyl file, the minimum methylation frequency, and whether you are observing upstream or downstream of a methylation site
file = input("Input the corrected bedmethyl file: ")
mf = float(input("Input the minimum methylation frequency for a site to be considered methylated: "))
s = input("Input whether the interest is upstream (-) or downstream (+): ")

# Provide names to upstream variables and integer values
if s == '-':
	stream = -1
	stream1 = "Upstream"
if s == '+':
	stream = 1
	stream1 = "Downstream"

# Read in the bedmethyl file
infile = open(file, 'r')
read = infile.readlines()

# Create dictionaries for high plus strand sites, high minus strand sites, and all sites
dict80minus = {}
dict80plus = {}
dictall = {}
for line in read:
	line = line.strip()
	coord = int(line.split('\t')[1])		# Pull the coordinate for each site
	score = float(line.split('\t')[10])		# Pull the methylfreq score for each site
	strand = line.split('\t')[5]			# Pull the strand for each site
	dictall[coord] = str(score) + '_' + strand	# Add all sites to the complete dictionary
	if score > mf:
		if strand == '-':			# Add the high scores to the respective dictionary
			dict80minus[coord] = score
		else:
			dict80plus[coord] = score


# Create counts for the genomic positions relative to the methylated site
count1 = 0
count2 = 0
count3 = 0
count4 = 0
count5 = 0
count6 = 0
count7 = 0
count8 = 0
count9 = 0
count10 = 0
first = 0
second = 0
third = 0
fourth = 0
fifth = 0
sixth = 0
seventh = 0
eighth = 0
ninth = 0
tenth = 0

# Create counts for the opposite strand
count1m = 0
count2m = 0
count3m = 0
count4m = 0
count5m = 0
count6m = 0
count7m = 0
count8m = 0
count9m = 0
count10m = 0
firstm = 0
secondm = 0
thirdm = 0
fourthm = 0
fifthm = 0
sixthm = 0
seventhm = 0
eighthm = 0
ninthm = 0
tenthm = 0

# First add the scores from the plus strand dictionary
for x in dict80plus:
	site1 = x + (1*stream)	# Convert the site location to either upstream or downstream
	if site1 in dictall:
		score = float(dictall[site1].split('_')[0])	# Search for the respective site in the complete dictionary
		strand = dictall[site1].split('_')[1]
		if strand == '+':				# Add the score and count for the respective strand
			count1 = count1 + 1
			first = first + score
		else:
			count1m = count1m + 1
			firstm = firstm + score
	site2 = x + (2*stream)					# Repeat for the next ten positions following the methylated site
	if site2 in dictall:
		score = float(dictall[site2].split('_')[0])
		strand = dictall[site2].split('_')[1]
		if strand == '+':
			count2 = count2 +1
			second = second + score
		else:
			count2m = count2m + 1
			secondm = secondm + score
	site3 = x + (3*stream)
	if site3 in dictall:
		score = float(dictall[site3].split('_')[0])
		strand = dictall[site3].split('_')[1]
		if strand == '+':
			count3 = count3 + 1
			third = third + score
		else:
			count3m = count3m + 1
			thirdm = thirdm + score
	site4 = x + (4*stream)
	if site4 in dictall:
		score = float(dictall[site4].split('_')[0])
		strand = dictall[site4].split('_')[1]
		if strand == '+':
			count4 = count4 + 1
			fourth = fourth + score
		else:
			count4m = count4m + 1
			fourthm = fourthm + score
	site5 = x + (5*stream)
	if site5 in dictall:
		score = float(dictall[site5].split('_')[0])
		strand = dictall[site5].split('_')[1]
		if strand == '+':
			count5 = count5 + 1
			fifth = fifth + score
		else:
			count5m = count5m + 1
			fifthm = fifthm + score
	site6 = x + (6*stream)
	if site6 in dictall:
		score = float(dictall[site6].split('_')[0])
		strand = dictall[site6].split('_')[1]
		if strand == '+':
			count6 = count6 + 1
			sixth = sixth + score
		else:
			count6m = count6m + 1
			sixthm = sixthm + score
	site7 = x + (7*stream)
	if site7 in dictall:
		score = float(dictall[site7].split('_')[0])
		strand = dictall[site7].split('_')[1]
		if strand == '+':
			count7 = count7 + 1
			seventh = seventh + score
		else:
			count7m = count7m + 1
			seventhm = seventhm + score
	site8 = x + (8*stream)
	if site8 in dictall:
		score = float(dictall[site8].split('_')[0])
		strand = dictall[site8].split('_')[1]
		if strand == '+':
			count8 = count8 + 1
			eighth = eighth + score
		else:
			count8m = count8m + 1
			eighthm = eighthm + score
	site9 = x + (9*stream)
	if site9 in dictall:
		score = float(dictall[site9].split('_')[0])
		strand = dictall[site9].split('_')[1]
		if strand == '+':
			count9 = count9 + 1
			ninth = ninth + score
		else:
			count9m = count9m + 1
			ninthm = ninthm + score
	site10 = x + (10*stream)
	if site10 in dictall:
		score = float(dictall[site10].split('_')[0])
		strand = dictall[site10].split('_')[1]
		if strand == '+':
			count10 = count10 + 1
			tenth = tenth + score
		else:
			count10m = count10m + 1
			tenthm = tenthm + score

# Add the scores from the sites that are located on the reverse strand
for x in dict80minus:
	site1 = x - (1*stream)		# Convert the site so it is on the correct stream
	if site1 in dictall:
		score = float(dictall[site1].split('_')[0])		# Search for the site in the complete dictionary
		strand = dictall[site1].split('_')[1]
		if strand == '-':
			count1 = count1 + 1				# Add the score and count to the respective position counts
			first = first + score
		else:
			count1m = count1m + 1
			firstm = firstm + score
	site2 = x - (2*stream)						# Repeat for the next ten sites
	if site2 in dictall:
		score = float(dictall[site2].split('_')[0])
		strand = dictall[site2].split('_')[1]
		if strand == '-':
			count2 = count2 + 1
			second = second + score
		else:
			count2m = count2m + 1
			secondm = secondm + score
	site3 = x - (3*stream)
	if site3 in dictall:
		score = float(dictall[site3].split('_')[0])
		strand = dictall[site3].split('_')[1]
		if strand == '-':
			count3 = count3 + 1
			third = third + score
		else:
			count3m = count3m + 1
			thirdm = thirdm + score
	site4 = x - (4*stream)
	if site4 in dictall:
		score = float(dictall[site4].split('_')[0])
		strand = dictall[site4].split('_')[1]
		if strand == '-':
			count4 = count4 + 1
			fourth = fourth + score
		else:
			count4m = count4m + 1
			fourthm = fourthm + score
	site5 = x - (5*stream)
	if site5 in dictall:
		score = float(dictall[site5].split('_')[0])
		strand = dictall[site5].split('_')[1]
		if strand == '-':
			count5 = count5 + 1
			fifth = fifth + score
		else:
			count5m = count5m + 1
			fifthm = fifthm + score
	site6 = x - (6*stream)
	if site6 in dictall:
		score = float(dictall[site6].split('_')[0])
		strand = dictall[site6].split('_')[1]
		if strand == '-':
			count6 = count6 + 1
			sixth = sixth + score
		else:
			count6m = count6m + 1
			sixthm = sixthm + score
	site7 = x - (7*stream)
	if site7 in dictall:
		score = float(dictall[site7].split('_')[0])
		strand = dictall[site7].split('_')[1]
		if strand == '-':
			count7 = count7 + 1
			seventh = seventh + score
		else:
			count7m = count7m + 1
			seventhm = seventhm + score
	site8 = x - (8*stream)
	if site8 in dictall:
		score = float(dictall[site8].split('_')[0])
		strand = dictall[site8].split('_')[1]
		if strand == '-':
			count8 = count8 + 1
			eighth = eighth + score
		else:
			count8m = count8m + 1
			eighthm = eighthm + score
	site9 = x - (9*stream)
	if site9 in dictall:
		score = float(dictall[site9].split('_')[0])
		strand = dictall[site9].split('_')[1]
		if strand == '-':
			count9 = count9 + 1
			ninth = ninth + score
		else:
			count9m = count9m + 1
			ninthm = ninthm + score
	site10 = x - (10*stream)
	if site10 in dictall:
		score = float(dictall[site10].split('_')[0])
		strand = dictall[site10].split('_')[1]
		if strand == '-':
			count10 = count10 + 1
			tenth = tenth + score
		else:
			count10m = count10m + 1
			tenthm = tenthm + score

# Print the results, i.e. the average methylation for each site surrounding a methylated adenine
print ("%s; Positive Strand" % (stream1))

print (first/count1, count1)
print (second/count2, count2)
print (third/count3, count3)
print (fourth/count4, count4)
print (fifth/count5, count5)
print (sixth/count6, count6)
print (seventh/count7, count7)
print (eighth/count8, count8)
print (ninth/count9, count9)
print (tenth/count10, count10)

print ("%s; Negative Strand" % (stream1))

print (firstm/count1m, count1m)
print (secondm/count2m, count2m)
print (thirdm/count3m, count3m)
print (fourthm/count4m, count4m)
print (fifthm/count5m, count5m)
print (sixthm/count6m, count6m)
print (seventhm/count7m, count7m)
print (eighthm/count8m, count8m)
print (ninthm/count9m, count9m)
print (tenthm/count10m, count10m)

