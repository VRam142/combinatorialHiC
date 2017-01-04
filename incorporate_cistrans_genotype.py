import os,sys,re

fhi = open(sys.argv[1]) #cis-trans

ratios = {}

for line in fhi:
	split = line.split()
	barcode = split[0]
	ratios[barcode] = split[1]

fhi.close()


fhi = open(sys.argv[2]) #percentages file

for line in fhi:
        split = line.split()
        barcode = "%s-%s" % (split[6], split[7])
	if barcode not in ratios:
		ratios[barcode] = 0
        print "%s\t%s\t%s" % (line.rstrip('\n'), ratios[barcode], 0) 
