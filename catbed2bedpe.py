import os,sys,re

def uncouple_reads(fhi, barcodes, fho, stats):
	'''Take a file handle of a sorted bed file, as well as a list of name -- barcode associations, and print out the associations.\
		Write out the statistics for mismatches to stats, an fho'''
	matched = 0
	unmatched = 0
	linkage = {} #Dict to associate names with barcodes and barcode status
	unlinked = {} #Dict to store unmatched barcode names
	for line in barcodes: #iterate through the barcode file
		name, bc1, bc2, bcterm = line.split()
		name_correct = name.split("/")[0]
		if bc1 == bc2: #Check to see that both halves of the barcode match
			barcode = bc1 #And store the barcode
			matched += 1
			linkage[name_correct] = [barcode, bcterm]
		else:
			barcode = "%s-%s" % (bc1, bc2)
			unlinked[name_correct] = [barcode, bcterm]
			unmatched += 1
	print >> stats, "<H3>Read Pairs with Matched Barcodes:%s</H3>\n<H3>Read Pairs with Unmatched Barcodes:%s</H3>" % (matched, unmatched)
	read_lb = [1,2,3,4,5,6]
	for line in fhi:
		chrom, fend, rend, name, mapq, strand, chrom_frag, fend_frag, rend_frag, frag_id, space, frag_strand, dist = line.split()
		name_correct = name.split("/")[0]
		if name_correct in linkage:
			if read_lb[3]  == name_correct:
				if chrom == read_lb[0]:
					if int(fend) == int(read_lb[1]) or int(rend) == int(read_lb[2]):
						print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %\
							(chrom, fend, rend, read_lb[0], read_lb[1], read_lb[2], name, mapq, read_lb[4], "+", "-", linkage[name_correct][0]\
								, linkage[name_correct][1], frag_id, dist, read_lb[9], read_lb[12])
					elif int(fend) < int(read_lb[2]):
						print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"  %\
							(chrom, fend, rend, read_lb[0], read_lb[1], read_lb[2], name, mapq, read_lb[4], strand, read_lb[5], linkage[name_correct][0]\
								, linkage[name_correct][1], frag_id, dist, read_lb[9], read_lb[12])
					else:
						print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"  %\
							(read_lb[0], read_lb[1], read_lb[2], chrom, fend, rend, name, read_lb[4], mapq, read_lb[5], strand, linkage[name_correct][0]\
								, linkage[name_correct][1], read_lb[9], read_lb[12], frag_id, dist)
				else:
					if chrom < read_lb[0]:
						print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"  %\
							(chrom, fend, rend, read_lb[0], read_lb[1], read_lb[2], name, mapq, read_lb[4], strand, read_lb[5], linkage[name_correct][0]\
								, linkage[name_correct][1], frag_id, dist, read_lb[9], read_lb[12])
					else:
						print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"  %\
							(read_lb[0], read_lb[1], read_lb[2], chrom, fend, rend, name, read_lb[4], mapq, read_lb[5], strand, linkage[name_correct][0]\
								, linkage[name_correct][1], read_lb[9], read_lb[12], frag_id, dist)
		elif name_correct in unlinked:
			if read_lb[3]  == name_correct:
				if chrom == read_lb[0]:
					if int(fend) == int(read_lb[1]) or int(rend) == int(read_lb[2]):
						print >> fho, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %\
							(chrom, fend, rend, read_lb[0], read_lb[1], read_lb[2], name, mapq, read_lb[4], "+", "-", unlinked[name_correct][0]\
								, unlinked[name_correct][1], frag_id, dist, read_lb[9], read_lb[12])
					elif int(fend) < int(read_lb[2]):
						print >> fho, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"  %\
							(chrom, fend, rend, read_lb[0], read_lb[1], read_lb[2], name, mapq, read_lb[4], strand, read_lb[5], unlinked[name_correct][0]\
								, unlinked[name_correct][1], frag_id, dist, read_lb[9], read_lb[12])
					else:
						print >> fho, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"  %\
							(read_lb[0], read_lb[1], read_lb[2], chrom, fend, rend, name, read_lb[4], mapq, read_lb[5], strand, unlinked[name_correct][0]\
								, unlinked[name_correct][1], read_lb[9], read_lb[12], frag_id, dist)
				else:
					if chrom < read_lb[0]:
						print >> fho, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"  %\
							(chrom, fend, rend, read_lb[0], read_lb[1], read_lb[2], name, mapq, read_lb[4], strand, read_lb[5], unlinked[name_correct][0]\
								, unlinked[name_correct][1], frag_id, dist, read_lb[9], read_lb[12])
					else:
						print >> fho, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"  %\
							(read_lb[0], read_lb[1], read_lb[2], chrom, fend, rend, name, read_lb[4], mapq, read_lb[5], strand, unlinked[name_correct][0]\
								, unlinked[name_correct][1], read_lb[9], read_lb[12], frag_id, dist)
		read_lb = [chrom, fend, rend, name_correct, mapq, strand, chrom_frag, fend_frag, rend_frag, frag_id, space, frag_strand, dist]

def main():
	fhi = open(sys.argv[1])
	barcodes = open(sys.argv[2])
	mismatch_out = open(sys.argv[3], 'w')
	outfile = open(sys.argv[4], 'w')
	uncouple_reads(fhi, barcodes, mismatch_out, outfile)
	fhi.close()
	barcodes.close()
	mismatch_out.close()
	outfile.close()

if __name__ == "__main__":
	main()
