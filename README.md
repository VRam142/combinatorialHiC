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
* List of valid chromosome ids in format "chrid\tchr_length" (combo_hg19_mm10.genomesize, but this could conceivably change depending on the analysis)

One will need to specify some of the paths in the wrapper script, namely, the paths to the Python scripts / barcode files provided in this repo, and the paths to the reference files.

Finally, the wrapper script takes raw, gzipped fastqs as input.

The wrapper script is run as follows:
```bash
#Run the sciHi-C analysis pipeline
bash scihic_pipe.sh inner_barcodes.txt [R1.fq.gz] [R2.fq.gz] outer_barcodes.txt [outfile_prefix]
```

As written, the wrapper script writes most files to local storage (/tmp/) to minimize network burden.

Barcode Association & Read Trimming
-----------------------------------
Illumina indexed reads were demultiplexed allowing 1 mismatch. Demultiplexed raw reads are provided @ dbGaP study accession phs001269.v1.p1 and GEO GSE84920.

We first adaptor clip all reads using the SeqPrep utility. Then, to obtain round 2 (i.e. terminal) barcodes, we use a custom Python script to iterate through both mates, compare the first 8 bases of each read against the 96 known barcode sequences, and then assign barcodes to each mate using a Levenshtein distance cutoff of 2. Reads “split” in this way are output such that the first 11 bases of each read, which derivefrom the custom barcoded Y adaptors, are removed. Mates where either terminal barcode went unidentified, or where the terminal barcodes did not match, are discarded.

For each resulting “split” pair of reads, the two reads are then scanned using a custom Python script to find the common portion of the bridge adaptor sequence. The 8 bases immediately 5’ of this sequence are isolated and compared against the 96 known bridge adaptor barcodes, again using a Levenshtein distance cutoff of 2. There are cases where the entire bridge adaptor, including both barcodes flanking the ligation junction, is encountered in one mate, and not the other. To account for these cases, we also isolate the 8 bases flanking the 3’ end of the common bridge adaptor sequence (when it is encountered within a read), reverse complement it, and compare the resulting 8-mer against the 96 known bridge adaptor barcodes. Output reads are then clipped to remove the bridge adaptor and all 3’ sequence. Barcodes flanking the ligation junction should match; again, mates where barcodes do not match, or where a barcode is not found are discarded. 

The custom Python script to filter terminal barcodes can be run separately as:
```bash
#Separate out terminal barcodes from reads
python inline_splitter.py /tmp/$fq_r1.clipped /tmp/$fq_r2.clipped $barcodes \
	/tmp/$fq_r1.split /tmp/$fq_r2.split 2> $outdir/splitting_stats.html
```

The result of this processing module are three files: filtered reads 1 and 2, and an “associations” file—a tab-delimited file where the name of each read passing the above filters and their associated barcode combination are listed.

Read Alignment, Read Pairing, & Barcode Assocation
--------------------------------------------------
As is standard for Hi-C reads, the resulting processed and filtered reads 1 and 2 were aligned separately using bowtie2/2.2.3 to a Burrows-Wheeler Index of the concatenated mouse (mm10) and human (hg19) genomes. Individual SAM files were then converted to BED format and filtered for alignments with MAPQ >= 30 using a combination of samtools, bedtools, and awk. Using bedtools closest along with a BED file of all DpnII sites in both genomes (generated using HiC-Pro), the closest DpnII site to each read was determined, after which BED files were concatenated, sorted on read ID using UNIX sort, and then processed using a custom Python script to generate a BEDPE format file where 5’ mates always precede 3’ mates, and where a simple Python dictionary is used to associate barcode combinations contained in the “associations” file with each pair of reads.  Reads were then sorted by barcode, read 1 chromosome, start, end, read 2 chromosome, start, and end using UNIX sort, and deduplicated using a custom Python script on the following criteria: reads were considered to be PCR duplicates if they were associated with the same cellular index, and if they comprised a ligation between the same two restriction sites as defined using bedtools closest.

Cellular Demultiplexing & Quality Analysis
------------------------------------------
When demultiplexing cells, we run two custom Python scripts. First, we generate a “percentages” file that includes the species purity of each cellular index, the coverage of each index, and the number of times a particular restriction fragment is observed once, twice, thrice, and four times. We also include the cis:trans ratio described above, and, if applicable, the fraction of homozygous alternate HeLa alleles observed. We use these percentages files to filter BEDPE files (see below) and generate, at any desired resolution, single cell matrices in long format (i.e. BIN1-BIN2-COUNT), with only the “upper diagonal” of the matrix included to reduce storage footprint. These matrices can then be converted to numpy matrices for visualization and further analysis.

The percentages file, and the filtration criteria used in the paper (excluding HeLa genotype, which is not generalizable) can be generated separately by running the following (taken from the wrapper script):
```bash
#Generate PERCENTAGES file
python calculate_cell_distro_long.py $bc_assoc.bedpe.mapq0.deduped.filtered > \
	$bc_assoc.deduped.percentages 2> $bc_assoc.deduped.REoccurrences
#Calculate single-cell length distribution
python single_cell_length_distros.py $bc_assoc.deduped.percentages\
	 $bedpe.mapq0.deduped.filtered $bc_assoc > $bc_assoc.single_cell_lengths
#Calculate cis-trans ratios from these distributions
python calculate_median_cistrans_ratio.py\
	 $bc_assoc.single_cell_lengths > $bc_assoc.cistrans.txt
#Incorporate cistrans ratios in the percentages file
python incorporate_cistrans_genotype.py $bc_assoc.cistrans.txt $bc_assoc.deduped.percentages >\
	 $bc_assoc.deduped.percentages.filterable
```

Filtration of Cellular Indices
------------------------------
We applied several filters to our resulting cellular indices to arrive at the cells analyzed in this study. We first removed all cellular indices with fewer than 1000 unique reads. We next filtered out all indices where the cis:trans ratio was lower than 1. Finally, for all experiments we removed cellular indices where less than 95% of reads aligned uniquely to either the mouse (mm10) or human (hg19) genomes. For all human cells from HAP1 and HeLa S3 mixing experiments (Libraries ML1, ML2, PL1, and PL2) further filtration by genotype was performed. For each cellular index, we examined all reads overlapping with known alternate homozygous sites in the HeLa S3 genome and computed the fraction of sites where the alternate allele is observed. We then drew cutoffs to filter out all cells where this fraction fell between 56% and 99%. We employ this filtering step purely as an additional, conservative measure, and note that this is not strictly necessary. The clear separation of two populations in data derived from library ML4 (nocadazole arrest experiment), where no genotype filtration was performed, illustrates this.

All of the experiments presented in Ramani et al employed a "programmed barcoding" strategy, where specific cell types were seeded in specific wells for the first round of barcoding (prior to mixing all cells together for second round barcoding). This enables easy cell type validation when analyzing sciHi-C data. To "label" the cells, use the appropriate assign_wellids Python script (e.g. assign_wellids_ML.py; assign_wellids_PL.py) as follows. In this example, assign_wellids_ML.py is used:
```bash
#Assign cell type IDs to list of valid cells
python assign_wellids_ML.py inner_barcodes.txt [valid_cells.list] > [valid_cells.list].labeled
```
Where [valid_cells.list] is a list of valid cells (generated by the wrapper script).

PCA of sciHi-C Data
-------------------
Single-cell matrices at interchromosomal contact resolution (log10 of contact counts) and 10 Mb resolution (binarized; 0 if absent, 1 if present) were vectorized and concatenated using custom Python scripts. Concatenation was performed such that redundant entries of each contact matrix (i.e. Cij and Cji) were only represented once. Resulting matrices, where rows represent single-cells and columns represent observed contacts, were then decomposed using the PCA function in scikit-learn. For interchromosomal matrices, entries for intrachromosomal contacts (i.e. the diagonal) were set to 0. For 10 Mb intrachromosomal matrices, all interchromosomal contacts were ignored and all entries Cij where | i – j | < 3 were set to zero.

Code to perform PCA is not included in the wrapper script, but is included in this repository. It can be run as follows:
```bash
#Generate interchromosomal interaction matrices (rows are single-cells, columns are normalized interaction counts)
python convertcell2val.py [valid_human_cells.list.labeled] [hg19.genomesize] 500000000\
	 [outfile_name] > [outfile_name_for_coverages]
#Perform PCA (requires matplotlib, numpy, scipy) on matrices. This example script takes in 8 positional arguments:
#the coverage outfiles from above for ML1, ML2, PL1, PL2, and the matrix files for those same experiments. It performs PCA
#on the resulting combined matrix, and then plots the result using matplotlib. The script will also provide the first three 
#PC vectors in txt format to plot using other utilities in the file [outfile.txt].
python PCA_example.py [ML1_coverage] [ML2_coverage] [PL1_coverage] [PL2_coverage]\
	 [ML1_matrix] [ML2_matrix] [PL1_matrix] [PL2_matrix] > [outfile.txt]
```

Calculation of Contact Probabilities in Single Cells
----------------------------------------------------
Methods to calculate the scaling probability within single cells were adapted from Imakaev, Fudenberg et al, Nature Methods and Sanborn, Rao et al, PNAS. A histogram of contact distances normalized by bin size was generated using logarithmically increasing bins (increasing by powers of 1.12n). We obtained the scaling coefficient by calculating the line of best fit for the log-log plot of this histogram between distances of 50 kb and 8 Mb. Shuffled controls were generated by randomly reassigning all cellular indices and repeating the above analysis; this importantly maintains the coverage distribution of the new set of simulated “single cells.”

Code to perform this analysis is not included in the wrapper script, but an example is included in this repository. It can be run as follows:
```bash
#Calculate scaling probabilities for single cells. Takes a valid_cells.list and a filtered BEDPE file (both generated from wrapper, can be filtered as desired)
#and generates a file with single-cell contact probabilities as function of genomic distance, and calculated fit to log-log plot used in paper
python calculate_scaling_example.py [valid_cells.list] [filtered_BEDPE.file] > [scaling_output.txt]
```
