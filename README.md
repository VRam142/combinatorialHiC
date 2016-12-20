Processing single cell combinatorial indexed Hi-C (sciHi-C) data
================================================================

The following README documents code & analyses performed in (and are largely reproduced from the methods section therein):

Ramani V, Deng X, Qiu R, Gunderson KL, Steemers FJ, Noble WS, Disteche CM, Duan Z, Shendure J. Massively Multiplex Single-Cell Hi-C.
Nature Methods (2017)

*All scripts and binaries are provided as is, without any warrenty and for use at your own risk. This is not the release of a software package. We are only providing this information and code in addition to a description of methods for making it easier to reproduce our analyses. We are __not__ providing any support for these scripts.* 

Running the wrapper script
--------------------------
The wrapper script has a few software dependencies. They are: 
* SeqPrep
* samtools
* bedtools
* bowtie2
* python / BioPython

The wrapper script also requires a few reference files, some of which are provided here in this repository. These include:
* List of inner barcodes (provided)
* List of outer barcodes (provided)
* Combined bowtie2 reference for hg19 and mm10 (not provided)
* BED file of all DpnII cutsites in the combined hg19/mm10 reference (not provided)
* List of valid chromosome ids in format "chrid\tchr_length" (not provided, as this could conceivably change depending on the analysis)

Finally, the wrapper script takes raw, gzipped fastqs as input.

The wrapper script is run as follows:
```bash
#Run the sciHi-C analysis pipeline
bash scihic_pipe.sh inner_barcodes.txt [R1.fq.gz] [R2.fq.gz] outer_barcodes.txt [outfile_prefix]
```

By default, the wrapper script writes most files to local storage (/tmp/) to minimize network burden.

Barcode Association & Read Trimming
-----------------------------------
Illumina indexed reads were demultiplexed allowing 1 mismatch. Demultiplexed raw reads are provided @ dbGaP study 21085 and GEO GSE84290.

We first adaptor clip all reads using the SeqPrep utility. Then, to obtain round 2 (i.e. terminal) barcodes, we use a custom Python script to iterate through both mates, compare the first 8 bases of each read against the 96 known barcode sequences, and then assign barcodes to each mate using a Levenshtein distance cutoff of 2. Reads “split” in this way are output such that the first 11 bases of each read, which derivefrom the custom barcoded Y adaptors, are removed. Mates where either terminal barcode went unidentified, or where the terminal barcodes did not match, are discarded.

For each resulting “split” pair of reads, the two reads are then scanned using a custom Python script to find the common portion of the bridge adaptor sequence. The 8 bases immediately 5’ of this sequence are isolated and compared against the 96 known bridge adaptor barcodes, again using a Levenshtein distance cutoff of 2. There are cases where the entire bridge adaptor, including both barcodes flanking the ligation junction, is encountered in one mate, and not the other. To account for these cases, we also isolate the 8 bases flanking the 3’ end of the common bridge adaptor sequence (when it is encountered within a read), reverse complement it, and compare the resulting 8-mer against the 96 known bridge adaptor barcodes. Output reads are then clipped to remove the bridge adaptor and all 3’ sequence. Barcodes flanking the ligation junction should match; again, mates where barcodes do not match, or where a barcode is not found are discarded. 

The result of this processing module are three files: filtered reads 1 and 2, and an “associations” file—a tab-delimited file where the name of each read passing the above filters and their associated barcode combination are listed.

Read Alignment, Read Pairing, & Barcode Assocation
--------------------------------------------------
As is standard for Hi-C reads, the resulting processed and filtered reads 1 and 2 were aligned separately using bowtie2/2.2.3 to a Burrows-Wheeler Index of the concatenated mouse (mm10) and human (hg19) genomes. Individual SAM files were then converted to BED format and filtered for alignments with MAPQ >= 30 using a combination of samtools, bedtools, and awk. Using bedtools closest along with a BED file of all DpnII sites in both genomes (generated using HiC-Pro), the closest DpnII site to each read was determined, after which BED files were concatenated, sorted on read ID using UNIX sort, and then processed using a custom Python script to generate a BEDPE format file where 5’ mates always precede 3’ mates, and where a simple Python dictionary is used to associate barcode combinations contained in the “associations” file with each pair of reads.  Reads were then sorted by barcode, read 1 chromosome, start, end, read 2 chromosome, start, and end using UNIX sort, and deduplicated using a custom Python script on the following criteria: reads were considered to be PCR duplicates if they were associated with the same cellular index, and if they comprised a ligation between the same two restriction sites as defined using bedtools closest.

Cellular Demultiplexing & Quality Analysis
------------------------------------------
When demultiplexing cells, we run two custom Python scripts. First, we generate a “percentages” file that includes the species purity of each cellular index, the coverage of each index, and the number of times a particular restriction fragment is observed once, twice, thrice, and four times. We also include the cis:trans ratio described above, and, if applicable, the fraction of homozygous alternate HeLa alleles observed. We use these percentages files to filter BEDPE files (see below) and generate, at any desired resolution, single cell matrices in long format (i.e. BIN1-BIN2-COUNT), with only the “upper diagonal” of the matrix included to reduce storage footprint. These matrices are then converted to numpy matrices for visualization and further analysis.

Filtration of Cellular Indices
------------------------------
We applied several filters to our resulting cellular indices to arrive at the cells analyzed in this study. We first removed all cellular indices with fewer than 1000 unique reads. We next filtered out all indices where the cis:trans ratio was lower than 1. Finally, for all experiments we removed cellular indices where less than 95% of reads aligned uniquely to either the mouse (mm10) or human (hg19) genomes. For all human cells from HAP1 and HeLa S3 mixing experiments (Libraries ML1, ML2, PL1, and PL2) further filtration by genotype was performed. For each cellular index, we examined all reads overlapping with known alternate homozygous sites in the HeLa S3 genome and computed the fraction of sites where the alternate allele is observed. We then drew cutoffs to filter out all cells where this fraction fell between 56% and 99%. We employ this filtering step purely as an additional, conservative measure, and note that this is not strictly necessary. The clear separation of two populations in data derived from library ML4 (nocadazole arrest experiment), where no genotype filtration was performed, illustrates this.

PCA of sciHi-C Data
-------------------
Single-cell matrices at interchromosomal contact resolution (log10 of contact counts) and 10 Mb resolution (binarized; 0 if absent, 1 if present) were vectorized and concatenated using custom Python scripts. Concatenation was performed such that redundant entries of each contact matrix (i.e. Cij and Cji) were only represented once. Resulting matrices, where rows represent single-cells and columns represent observed contacts, were then decomposed using the PCA function in scikit-learn. For interchromosomal matrices, entries for intrachromosomal contacts (i.e. the diagonal) were set to 0. For 10 Mb intrachromosomal matrices, all interchromosomal contacts were ignored and all entries Cij where | i – j | < 3 were set to zero.

Calculation of Contact Probabilities in Single Cells
----------------------------------------------------
Methods to calculate the scaling probability within single cells were adapted from Imakaev, Fudenberg et al, Nature Methods and Sanborn, Rao et al, PNAS. A histogram of contact distances normalized by bin size was generated using logarithmically increasing bins (increasing by powers of 1.12n). We obtained the scaling coefficient by calculating the line of best fit for the log-log plot of this histogram between distances of 50 kb and 8 Mb. Shuffled controls were generated by randomly reassigning all cellular indices and repeating the above analysis; this importantly maintains the coverage distribution of the new set of simulated “single cells.”

