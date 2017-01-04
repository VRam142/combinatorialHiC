#Script for analyzing sciHi-C fastqs.
import os,sys,re
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools as it
import gzip as gz
from Levenshtein import hamming

def addWildcards(string):
	foo= ""
	for i in range(len(string)):
        	foo+=string[:i] + '.' + string[i+1:] + "|"
	return foo[:-1]

def rc(seq):
	rcd = ""
	rcs = {"A":"T","C":"G","T":"A","G":"C", "N":"N"}
	for char in seq[::-1]:
		rcd+=rcs[char]
	return rcd

def checkHamming(barcodes,barcode):
	for bc in barcodes:
		match = False
		hd = hamming(barcode, bc)
		if hd <= 2:
			match = True
			barcode = bc
			break
	return (match, barcode)

#@profile
def associateBarcodes(fhi_r1, fhi_r2, fho_r1, fho_r2, internal_bc):
	'''Takes a Fastq filehandle in and filehandle out for each read; writes barcode-clipped fastqs to\
	filehandle out and returns a dictionary with barcode pairs associated with given read names'''
	seq_r1 = FastqGeneralIterator(fhi_r1)
	seq_r2 = FastqGeneralIterator(fhi_r2)
	seqzip = it.izip(seq_r1, seq_r2) #zip up the two fastq iterators to increase performance
	bc_seq = {}
	bc1s = [] #Container to store internal barcode IDs
	bridge = ''
	#Create a list of all possible internal barcodes (necessary
	#for hamming distance calculations)
	for line in internal_bc:
		bc = line.rstrip('\n')
        	strip = 'GATC' + bc + 'GAATTC'
        	bc1s.append(bc)
	bridge_re = re.compile('GATC[ACGT]{8}GAATTC')
    #We want to find all reads where
    #the internal adaptor is present
    #in either R1 or R2. Accepting that we're
    #throwing out some stuff because of mismatches,
    ##let's just use the re module to do this.
	for i in seqzip:
		title1, seq1, qual1 = i[0]
		title2, seq2, qual2 = i[1]
		new_title = title1.split()[0]
		m1 = bridge_re.search(seq1)
		m2 = bridge_re.search(seq2)
		terminal_bc = title1.split('&')[1]
		if m1:
			if m2:
				pos1 = m1.start()
				pos2 = m2.start()
				trim_seq1 = seq1[:pos1]
				trim_seq2 = seq2[:pos2]
				trim_qual1 = qual1[:pos1]
				trim_qual2 = qual2[:pos2]
				#Check for length of the read
				if len(trim_seq1) < 30 or len(trim_seq2) < 30:
					continue #and continue if either one is too short
				barcode_seq1 = seq1[pos1 + 4:pos1 + 12]
				barcode_seq2 = seq2[pos2 + 4:pos2 + 12]
				ham_bc1 = checkHamming(bc1s, barcode_seq1)
				ham_bc2 = checkHamming(bc1s, barcode_seq2)
				if ham_bc1[0]:
					if ham_bc2[0]:
						print >> fho_r1, "@%s\n%s\n+\n%s" % (new_title, trim_seq1, trim_qual1)
						print >> fho_r2, "@%s\n%s\n+\n%s" % (new_title, trim_seq2, trim_qual2)
						bc1 = ham_bc1[1]
						bc2 = ham_bc2[1]
						print "%s\t%s\t%s\t%s" % (new_title, bc1, bc2, terminal_bc)
			else:
				pos1 = m1.start()
				pos1_end = m1.end()
				trim_seq1 = seq1[:pos1]
				trim_qual1 = qual1[:pos1]
				#Check for length of the rea
				if len(trim_seq1) < 30:
					continue #and continue if either one is too short
				if len(seq1[pos1_end:]) < 12:
					continue
				barcode_seq2 = rc(seq1[pos1_end:pos1_end + 8])
				barcode_seq1 = seq1[pos1 + 4:pos1 + 12]
				ham_bc1 = checkHamming(bc1s, barcode_seq1)
				ham_bc2 = checkHamming(bc1s, barcode_seq2)
				if ham_bc1[0]:
                    			if ham_bc2[0]:
                    				print >> fho_r1, "@%s\n%s\n+\n%s" % (new_title, trim_seq1, trim_qual1)
						print >> fho_r2, "@%s\n%s\n+\n%s" % (new_title, seq2, qual2)
						bc1 = ham_bc1[1]
						bc2 = ham_bc2[1]
						print "%s\t%s\t%s\t%s" % (new_title, bc1, bc2, terminal_bc)
		else:
			if m2:
				pos2 = m2.start()
				pos2_end = m2.end()
				trim_seq2 = seq2[:pos2]
				trim_qual2 = qual2[:pos2]
				#Check for length of the read
#				print seq2[pos2_end:]
				if len(trim_seq2) < 30:
					continue #and continue if either one is too short
				if len(seq2[pos2_end:]) < 12:
					continue #and continue if either one is too short
				barcode_seq1 = rc(seq2[pos2_end:pos2_end + 8])
				barcode_seq2 = seq2[pos2 + 4:pos2 + 12]
				ham_bc1 = checkHamming(bc1s, barcode_seq1)
				ham_bc2 = checkHamming(bc1s, barcode_seq2)
				if ham_bc1[0]:
					if ham_bc2[0]:
						print >> fho_r1, "@%s\n%s\n+\n%s" % (new_title, seq1, qual1)
						print >> fho_r2, "@%s\n%s\n+\n%s" % (new_title, trim_seq2, trim_qual2)
        			    		bc1 = ham_bc1[1]
						bc2 = ham_bc2[1]
						print "%s\t%s\t%s\t%s" % (new_title, bc1, bc2, terminal_bc)

def main():
	fhi_bc = open(sys.argv[1])
	fhi_r1 = open(sys.argv[2])
	fhi_r2 = open(sys.argv[3])
	fho_r1 = open(sys.argv[4], 'w')
	fho_r2 = open(sys.argv[5], 'w')
	assoc = associateBarcodes(fhi_r1, fhi_r2, fho_r1, fho_r2, fhi_bc)
	fhi_bc.close()
	fhi_r1.close()
	fhi_r2.close()
	fho_r1.close()
	fho_r2.close()

if __name__ == "__main__":
	main()
