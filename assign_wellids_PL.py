import os,sys,re

well_ids = {}

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
	if pos == "A" or pos == "B":
		colour = "HeLa"
	elif pos == "C" or pos == "D":
		colour = "HAP1"
	elif pos == "E" or pos == "F":
		colour = "Patski"
	else:
		colour = "MEF"
	print "%s\t%s" % (line.rstrip('\n'), colour)	

fhi.close()
