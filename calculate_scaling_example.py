import numpy as np
from pylab import *
from scipy import optimize
from numpy import log10
from sklearn import linear_model
from math import *
from math import pow
from random import shuffle
from collections import Counter

#Goal: Examine the scaling coefficients of power law decay in single cells. Here we're generating the null model by
#shuffling the barcode combos. This is a JANKY implementation that isn't optimized for memory usage or runtime, but
#it works.

def initialize_bins():
    '''Initialize log bins that will be used'''
    bins = []
    i = 5e3
    bins.append(i)
    index = 1
    while i < 1.6e8:
        i_new = int(i * pow(1.12, index))
        index += 1
        if i_new == i:
            i = i_new
            continue
        else:
            i = i_new
            bins.append(i)
    return bins

def scaling_histogram(bedpe, log_bins):
    '''JANKY way of binning counts in log bins. One should / could write a better / more generalizable implementation 
    that uses interval trees.'''
    histogram = {}
    total = Counter()
    for line in bedpe:
        split = line.split()
        barcode = "%s-%s" % (split[11], split[12])
        if barcode not in histogram:
            histogram[barcode] = Counter()
            for i in log_bins:
                histogram[barcode][i] = 0
        chrid1 = split[0]
        chrid2 = split[3]
        pos1 = (int(split[1]) + int(split[2])) / 2
        pos2 = (int(split[4]) + int(split[5])) / 2
        distance = abs(pos1 - pos2) #Consider the distance the difference between fragment midpoints
        if distance >= log_bins[-1]: continue #Ignore distances longer than the max distance
        if chrid1 != chrid2: continue #obviously, trans contacts are no good here
        frag1 = int(split[13].split('_')[-1])
        frag2 = int(split[15].split('_')[-1])
        if abs(frag1 - frag2) < 2: continue #Ignore fragments closer than 2 restriction sites away
        #Just...don't look:
        if log_bins[0] <= distance < log_bins[1]:
            histogram[barcode][log_bins[0]] += 1
        if log_bins[1] <= distance < log_bins[2]:
            histogram[barcode][log_bins[1]] += 1            
        if log_bins[2] <= distance < log_bins[3]:
            histogram[barcode][log_bins[2]] += 1
        if log_bins[3] <= distance < log_bins[4]:
            histogram[barcode][log_bins[3]] += 1  
        if log_bins[4] <= distance < log_bins[5]:
            histogram[barcode][log_bins[4]] += 1
        if log_bins[5] <= distance < log_bins[6]:
            histogram[barcode][log_bins[5]] += 1            
        if log_bins[6] <= distance < log_bins[7]:
            histogram[barcode][log_bins[6]] += 1
        if log_bins[7] <= distance < log_bins[8]:
            histogram[barcode][log_bins[7]] += 1    
        if log_bins[8] <= distance < log_bins[9]:
            histogram[barcode][log_bins[8]] += 1            
        if log_bins[9] <= distance < log_bins[10]:
            histogram[barcode][log_bins[9]] += 1
        if log_bins[10] <= distance < log_bins[11]:
            histogram[barcode][log_bins[10]] += 1
        if log_bins[11] <= distance < log_bins[12]:
            histogram[barcode][log_bins[11]] += 1            
        if log_bins[12] <= distance < log_bins[13]:
            histogram[barcode][log_bins[12]] += 1
        if log_bins[13] <= distance < log_bins[14]:
            histogram[barcode][log_bins[13]] += 1
    return histogram

def plot_scaling(counter_object, log_bins):
    '''Normalize contact probabilities and calculate slope between bins 6 and 9.'''
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y: (y - fitfunc(p, x))
    pinit = [1.0, -1.0]
    x = []
    y = []
    for i in range(len(log_bins) - 1):
        x.append(log_bins[i])
        value = float(counter_object[log_bins[i]]) / float(abs(log_bins[i] - log_bins[i + 1]))
        y.append(value)
    area = 0
    for i in range(len(x)):
        area += y[i]
        scaling_factor = float(1) / float(area)
    newx = np.array(x[:-1])
    newy = np.array(y[:-1])
    logx = newx
    logy = newy
    for i in range(len(newx)):
        logx[i] = log10(newx[i])
        logy[i] = log10(scaling_factor * newy[i])
    out = optimize.leastsq(errfunc, pinit, args=(logx[6:10], logy[6:10]), full_output=1)
    pfinal = out[0]   
    index = pfinal[1]
    amp = np.exp(pfinal[0])    
    logy_fit = []
    for i in logx:
        new = pfinal[0] + pfinal[1] * i
        logy_fit.append(new)
    return (index, logx, logy)

bedpe = open(sys.argv[1]) # Filtered BEDPE file that you want to process 
log_bins = initialize_bins() #Initialize these bins
valid_human_cells = open(sys.argv[2]) #Valid cells list
valid_barcodes = {} #Keep track of valid barcodes
total_barcodes = [] #Used to generate a list of barcodes for shuffling
valid_reads = [] #These lists
random_reads = [] #get REALLY BIG
for line in valid_human_cells:
    split = line.split() #split the file
    barcode = split[0].split('_')[2] #barcode
    coverage = int(split[0].split('_')[1]) #cellular coverage
    if coverage < 5000: continue #Coverage cutoff
    valid_barcodes[barcode] = split[1] #Can change this to whatever you'd like (right now it is the cell type label)

#More janky stuff: storing all of the reads in a list to generate the shuffled null
for line in bedpe:
        split = line.split()
        barcode = "%s-%s" % (split[11], split[12])
        if barcode not in valid_barcodes: continue
        else:
            valid_reads.append(line.rstrip('\n'))
            total_barcodes.append(barcode)

shuffle(total_barcodes) #shuffle the barcodes
i = 0

#The jank continues: assign the shuffled barcodes to the existing reads. Note 
#that this generates REALLY BIG LISTS
for read in valid_reads:
    split = read.split()
    barcode1 = total_barcodes[i].split('-')[0]
    barcode2 = total_barcodes[i].split('-')[1]
    new_read = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
        (split[0], split[1], split[2], split[3], split[4], split[5], split[6], split[7], split[8], split[9], split[10], barcode1, barcode2, split[13], split[14], split[15], split[16])
    random_reads.append(new_read)
    i+=1    

human_histogram = scaling_histogram(valid_reads, log_bins) #Calculate scaling for the actual data
random_histogram = scaling_histogram(random_reads, log_bins) #And the shuffled null
    
#print everything so that you can make pretty plots
#Format is xcoord<\t>ycoord<\t>barcode<\t>cell_type<\t>coverage<\t>slope<\t>observed/shuffled
for barcode in valid_barcodes:
    index, logx, logy = plot_scaling(human_histogram[barcode], log_bins) 
    rindex,rlogx,rlogy = plot_scaling(random_histogram[barcode], log_bins)
    coverage = sum(human_histogram[barcode].values())
    rcoverage = sum(random_histogram[barcode].values())
    for i in range(len(logx)):
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (logx[i], logy[i], barcode, valid_barcodes[barcode], coverage, index, "Observed")
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rlogx[i], rlogy[i], barcode, valid_barcodes[barcode], rcoverage, rindex, "Shuffled")

bedpe.close()
valid_human_cells.close()
