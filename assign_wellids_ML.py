import os,sys,re

well_ids = {}
well_ids2 = {}

wells = []
for j in range(1,13):
	for i in ['A','B','C','D','E','F','G','H']:
		wells.append("%s\t%s" % (i,j))

fhi = open(sys.argv[1])

#print wells
index = 0
for barcode in fhi:
	split = barcode.split()
	well_ids[split[0]] = wells[index]
	index += 1

fhi.close()

fhi = open(sys.argv[2])

for line in fhi:
	split = line.split('_')
	species = split[0]
	barcode = split[2].split('-')[0]
	pos = well_ids[barcode][0]
	if species == "human":
		if pos == "A" or pos == "B":
			colour = "HeLa"
		elif pos == "C" or pos == "D":
			colour = "HeLa"
		elif pos == "E" or pos == "F":
			colour = "HAP1"
		else:
			colour = "HAP1"
	else:
		if pos == "A" or pos == "B":
			colour = "Patski"
		elif pos == "C" or pos == "D":
			colour = "Patski"
		elif pos == "E" or pos =="F":
			colour = "MEF"
		else:
			colour = "MEF"
	print "%s\t%s" % (line.rstrip('\n'), colour)	

fhi.close()
