import os,sys,re

def define_valid_chroms(genome_file):
	valid_chroms = {}
	for line in genome_file:
		split = line.split()
		valid_chroms[split[0]] = True
	return valid_chroms

def filter_bed(bedpe, valid_chroms):
	for line in bedpe:
		split = line.split()
		chrid1 = split[0]
		chrid2 = split[3]
		if chrid1 in valid_chroms and chrid2 in valid_chroms:
			print line.rstrip('\n')

def main():
	genome_file = open(sys.argv[1])
	bedpe = open(sys.argv[2])
	valid_chroms = define_valid_chroms(genome_file)
	filter_bed(bedpe, valid_chroms)
	genome_file.close()
	bedpe.close()

if __name__ == "__main__":
	main()
