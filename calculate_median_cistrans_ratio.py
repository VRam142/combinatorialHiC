import os,sys,re
from collections import Counter

fhi = open(sys.argv[1])

barcodes = {}

for line in fhi:
	split = line.split('\t')
	barcode = split[1]
	if barcode not in barcodes:
		barcodes[barcode] = Counter()
	if split[2] == ">20kb" or split[2] == "<1kb" or split[2] == "1kb-20kb":
#	if split[2] == ">20kb":
		barcodes[barcode][">20kb"] += float(split[3]) / float(split[4])
	if split[2] == "Interchromosomal":
		barcodes[barcode][split[2]] = float(split[3]) / float(split[4])

for barcode in barcodes:
	cistrans_ratio = barcodes[barcode][">20kb"] / barcodes[barcode]["Interchromosomal"]
	print "%s\t%s" % (barcode, cistrans_ratio)
	
