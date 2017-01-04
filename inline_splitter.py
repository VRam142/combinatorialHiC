import os,sys,re
import itertools as it
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Levenshtein import hamming
import gzip as gz

def checkHamming(barcodes,barcode):
	'''Given a list of barcodes, check that the given barcode is within edit\
		distance 2 to any of the list of barcodes'''
	for bc in barcodes:
		match = False
		hd = hamming(barcode, bc)
		if hd <= 2:
			match = True
			barcode = bc
			break
	return (match, barcode)

def split_fastqs(r1, r2, r1_o, r2_o, barcodes):
	'''Given paired-end fastq data, split reads based off of an inline 10 (or 11) mer barcode'''
	mismatch = 0
	not_found_R1 = 0
	not_found_R2 = 0
	reads = 0
	fqr1 = FastqGeneralIterator(r1)
	fqr2 = FastqGeneralIterator(r2)
	seqzip = it.izip(fqr1, fqr2) #Zip up the two iterators for expediency
	for pairs in seqzip:
		title1, seq1, qual1 = pairs[0]
		title2, seq2, qual2 = pairs[1]
		barcode1 = seq1[:8] #Just look in read 1--barcodes SHOULD be the same on both ends of the molecule)
		barcode2 = seq2[:8] #Check barcode 2 as well; print out how many times they disagree
		test1 = checkHamming(barcodes, barcode1)
		test2 = checkHamming(barcodes, barcode2)
		if test1[0]:
			if test2[0]:
			#if the barcodes match, print out the trimmed / split reads to new files
				if test1[1] == test2[1]:
					print >> r1_o, "@%s&%s\n%s\n+\n%s" % (title1,test1[1], seq1[11:], qual1[11:])
					print >> r2_o, "@%s&%s\n%s\n+\n%s" % (title2,test2[1], seq2[11:], qual2[11:])
				else:
					mismatch += 1
			else: #If there isn't a match in R1
				not_found_R2 += 1
		elif test2[0]:
			not_found_R1 += 1
		else:
			not_found_R1 += 1
			not_found_R2 += 1
		reads += 1
	out_error0 = "<H3>Total number of reads:%s</H3>" % reads
	out_error1 = "<H3>Total number of barcode mismatches:%s</H3>" % mismatch
	out_error2 = "<H3>Total number of missed R1 barcodes:%s</H3>" % not_found_R1
	out_error3 = "<H3>Total number of missed R2 barcodes:%s</H3>" % not_found_R2
	sys.stderr.write(out_error0 + '\n' + out_error1 + '\n' + out_error2 + '\n' + out_error3 + '\n')

def main():
	r1 = gz.open(sys.argv[1]) #R1 filehandle
	r2 = gz.open(sys.argv[2]) #R2 filehandle
	barcode_fhi = open(sys.argv[3]) #Barcode filehandle
	#Barcodes
	barcodes = []
	for line in barcode_fhi:
		barcodes.append(line.split()[0])
	r1_o = open(sys.argv[4], 'w') #R1 filehandle out
	r2_o = open(sys.argv[5], 'w') #R2 filehandle out
	split_fastqs(r1, r2, r1_o, r2_o, barcodes) #Split some fastqs, yo
	#Close the filehandles
	r1.close()
	r2.close()
	r1_o.close()
	r2_o.close()
	barcode_fhi.close()

if __name__ == "__main__":
	main()
