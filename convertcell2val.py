import os,sys,re
import numpy as np
import scipy as sp
from math import pow
from math import log10

def define_max_index(chromsizes, resolution):
    index = 0
    lines = chromsizes.readlines()
    chromid = {}
    for line in lines:
        chromname, length = line.split()
        for i in range(0,int(length),resolution):
            chromid[index] = chromname
            index += 1
    return index, chromid

def cell2val(matrix_list, index):
    pos_index = 0
    matpos = {}
    for i in range(0,index):
        for j in range(0,index):
            if i <= j:
                matpos[(i,j)] = pos_index
                pos_index += 1
    valmatrix = np.zeros(pos_index)
    for matrix in matrix_list:
        newvec = np.zeros(pos_index)
        matrix_open = open(matrix)
        for line in matrix_open:
            bin1, bin2, raw, norm, chrom1, chrom2 = line.split()
#            if chrom1 == chrom2: #Ignore intrachromosomal contacts for now; could also implement as a check to see how close bin1 and bin2 are
#                norm = np.NAN
            coord = (int(bin1), int(bin2))
            pos = matpos[coord]
	    if bin1 == bin2 or raw == 0:
                newvec[pos] = 0 #set to which value you want to compute cors for (root normalized coverage, or raw coverage)
            else:
		newvec[pos] = log10(float(raw))
        valmatrix = np.vstack((valmatrix, newvec))
        matrix_open.close()
    valmatrix = np.delete(valmatrix, 0, 0)
    return valmatrix

def main():
    matrices = open(sys.argv[1])    #positional argument 1 --> list of matrix filenames
    genome_size = open(sys.argv[2]) #positional argument 2 --> genome size file
    resolution = int(sys.argv[3])   #positional argument 3 --> resolution of matrices
    index, chromid = define_max_index(genome_size, resolution)
    matrix_list = []
    matrix_coverage = []
    for line in matrices:
	cell_identity = line.split()[1]
        coverage = line.split("_")[1]
	zipped = "%s\t%s" % (coverage, cell_identity)
        matrix_coverage.append(zipped)
        matrix_list.append(line.split()[0])
    for cov in matrix_coverage:
        print cov
    valmatrix = cell2val(matrix_list, index)
    np.savetxt("%s" % (sys.argv[4]), valmatrix, delimiter=',')
    matrices.close()

if __name__=="__main__":
    main()
