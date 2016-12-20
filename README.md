Processing single cell combinatorial indexed Hi-C (sciHi-C) data
================================================================

The following README documents code & analyses performed in (and are largely reproduced from the methods section therein):

Ramani V, Deng X, Qiu R, Gunderson KL, Steemers FJ, Noble WS, Disteche CM, Duan Z, Shendure J. Massively Multiplex Single-Cell Hi-C.
Nature Methods (2017)

*All scripts and binaries are provided as is, without any warrenty and for use at your own risk. This is not the release of a software package. We are only providing this information and code in addition to a description of methods for making it easier to reproduce our analyses. We are __not__ providing any support for these scripts.* 

Running the wrapper script
--------------------------
The wrapper script has a few software dependencies. They are: 
SeqPrep
samtools
bedtools
bowtie2
python / BioPython

The wrapper script also requires a few reference files, some of which are provided here in this repository. These include:
List of inner barcodes (provided)
List of outer barcodes (provided)
Combined bowtie2 reference for hg19 and mm10 (not provided)
BED file of all DpnII cutsites in the combined hg19/mm10 reference (not provided)
List of valid chromosome ids in format "chrid\tchr_length" (not provided, as this could conceivably change depending on the analysis)

Finally, the wrapper script takes raw, gzipped fastqs as input.

The wrapper script is run as follows:
```bash
#Run the sciHi-C analysis pipeline
bash scihic_pipe.sh inner_barcodes.txt [R1.fq.gz] [R2.fq.gz] outer_barcodes.txt [outfile_prefix]
```

By default, the wrapper script writes most files to local storage (/tmp/) to minimize network burden.

Barcode Association & Read Trimming
-----------------------------------
Illumina indexed reads were demultiplexed allowing 1 mismatch. Demultiplexed raw reads are provided @ dbGaP and GEO.

We first adaptor clip all reads using the SeqPrep utility. Then, to obtain round 2 (i.e. terminal) barcodes, we use a custom Python script to iterate through both mates, compare the first 8 bases of each read against the 96 known barcode sequences, and then assign barcodes to each mate using a Levenshtein distance cutoff of 2. Reads “split” in this way are output such that the first 11 bases of each read, which derivefrom the custom barcoded Y adaptors, are removed. Mates where either terminal barcode went unidentified, or where the terminal barcodes did not match, are discarded.

For each resulting “split” pair of reads, the two reads are then scanned using a custom Python script to find the common portion of the bridge adaptor sequence. The 8 bases immediately 5’ of this sequence are isolated and compared against the 96 known bridge adaptor barcodes, again using a Levenshtein distance cutoff of 2. There are cases where the entire bridge adaptor, including both barcodes flanking the ligation junction, is encountered in one mate, and not the other. To account for these cases, we also isolate the 8 bases flanking the 3’ end of the common bridge adaptor sequence (when it is encountered within a read), reverse complement it, and compare the resulting 8-mer against the 96 known bridge adaptor barcodes. Output reads are then clipped to remove the bridge adaptor and all 3’ sequence. Barcodes flanking the ligation junction should match; again, mates where barcodes do not match, or where a barcode is not found are discarded. 

The result of this processing module are three files: filtered reads 1 and 2, and an “associations” file—a tab-delimited file where the name of each read passing the above filters and their associated barcode combination are listed.

Read Alignment, Read Pairing, & Barcode Assocation
--------------------------------------------------
As is standard for Hi-C reads, the resulting processed and filtered reads 1 and 2 were aligned separately using bowtie2/2.2.3 to a Burrows-Wheeler Index of the concatenated mouse (mm10) and human (hg19) genomes. Individual SAM files were then converted to BED format and filtered for alignments with MAPQ >= 30 using a combination of samtools, bedtools, and awk. Using bedtools closest along with a BED file of all DpnII sites in both genomes (generated using HiC-Pro), the closest DpnII site to each read was determined, after which BED files were concatenated, sorted on read ID using UNIX sort, and then processed using a custom Python script to generate a BEDPE format file where 5’ mates always precede 3’ mates, and where a simple Python dictionary is used to associate barcode combinations contained in the “associations” file with each pair of reads.  Reads were then sorted by barcode, read 1 chromosome, start, end, read 2 chromosome, start, and end using UNIX sort, and deduplicated using a custom Python script on the following criteria: reads were considered to be PCR duplicates if they were associated with the same cellular index, and if they comprised a ligation between the same two restriction sites as defined using bedtools closest.

Cellular Demultiplexing & Quality Analysis
------------------------------------------


When demultiplexing cells, we run two custom Python scripts. First, we generate a “percentages” file that includes the species purity of each cellular index, the coverage of each index, and the number of times a particular restriction fragment is observed once, twice, thrice, and four times. We also include the cis:trans ratio described above, and, if applicable, the fraction of homozygous alternate HeLa alleles observed. We use these percentages files to filter BEDPE files (see below) and generate, at any desired resolution, single cell matrices in long format (i.e. BIN1-BIN2-COUNT), with only the “upper diagonal” of the matrix included to reduce storage footprint. These matrices are then converted to numpy matrices for visualization and further analysis.


-----------------------

For nuclesome peak calling, the L-WPS is locally adjusted to a running median of zero in 1 kb windows and smoothed using a Savitzky-Golay filter ([Savitzky and Golay, 1964](http://pubs.acs.org/doi/abs/10.1021/ac60214a047)) (window size 21, 2nd order polynomial). The L-WPS track is then segmented into above-zero regions (allowing up to 5 consecutive positions below zero). If the resulting region is 50-150 bp, we identify the median L-WPS value of that region and search for the maximum-sum contiguous window above the median. We report the start, end and center coordinates of this window as the “peak,” or local maximum of nucleosome protection. All calculations involving distances between peaks are based on these center coordinates. A score for each peak is determined as the distance between maximum value in the window and the average of the two adjacent L-WPS minima neighboring the region. If the identified region is 150-450 bp, we apply the same above median contiguous window approach, but only report those windows that are 50-150 bp. For score calculation of multiple windows derived from 150-450 bp regions, we set the neighboring minima to zero. We discard regions <50 bp or >450 bp.

Peak calling is implemented in `callPeaks.py` and expects WIG on STDIN:

```bash
# Note that EXTREGION should be larger than REGION (by the maximum read length, i.e. 180) to prevent skipping alignments
./samtools view -u -m 120 -M 180 BAMFILE.bam $EXTREGION | ./FilterUniqueBAM.py -p | ./extractReadStartsFromBAM2Wig.py -p -r $REGION -w 120 -c OFF -s OFF | ./callPeaks.py -s > calls.bed

# or from bigWig (http://hgdownload.cse.ucsc.edu/admin/exe/) and save as block-gzip compressed (http://www.htslib.org/doc/tabix.html) BED file:

bigWigToWig -chrom=chr1 -start=12000000 -end=13000000 ${SAMPLE}.bw /dev/stdout | ./callPeaks.py -s | bgzip -c > calls.bed.gz
```

Again note the variables `EXTREGION` and `REGION` above. When piping reads from a region, we might miss the forward read of a read pair and the provided scripts usually only extract information from the first read of a pair. Thus, `EXTREGION` should include additional bases around the region (i.e. 200bp or another value that guarantees that all read insert sizes are included; here 180bp would be sufficient; e.g. if REGION is 1:1000-2000, EXTREGION should be 1:800-2200).


Bed files can be converted to bigBed using the above mentioned [UCSC tools](http://hgdownload.cse.ucsc.edu/admin/exe/). We uploaded bigBed (bb) files of our nucleosome calls to GEO for [BH01](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71378/suppl/GSE71378_BH01.bb), [IH01](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71378/suppl/GSE71378_IH01.bb), [IH02](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71378/suppl/GSE71378_IH02.bb), [CH01](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71378/suppl/GSE71378_CH01.bb), and [CA01](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71378/suppl/GSE71378_CA01.bb).

Analysis of DHS sites
---------------------

DHS peaks for 349 primary tissue and cell line samples were downloaded from http://www.uwencode.org/proj/Science_Maurano_Humbert_et_al/data/all_fdr0.05_hot.tgz. Samples derived from fetal tissues (233 of 349) were removed from the analysis as they behaved inconsistently within tissue type, possibly because of unequal representation of cell types within each tissue sample. 116 DHS callsets from a variety of cell lineages were retained for analysis. For the midpoint of each DHS peak in a particular set, the nearest upstream and downstream calls in the CH01 callset were identified, and the distance between the centers of those two calls was calculated. The distribution of all such distances was visualized for each DHS peak callset using a smoothed density estimate calculated for distances between 0 and 500 bp.

See the folder `DHS` for more details and supporting files.  

Analysis of A/B compartments
----------------------------

We downloaded A/B compartment calls with 100 kb resolution from [GSE63525](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525), specifically the file `GSE63525_GM12878_subcompartments.bed.gz`. A/B segmentation calls were compared to peak spacing within equivalent 100 kb genomic windows.  In a separate analysis, peak density was compared to windowed GC content in 10 kb bins.

See the folder `peak_density` for more details and supporting files.

Annotation of TFBSs and genomic features
----------------------------------------

We downloaded clustered FIMO (motif-based) intervals ([Grant et al., 2011](http://www.ncbi.nlm.nih.gov/pubmed/21330290); [Maurano et al., 2012](http://www.ncbi.nlm.nih.gov/pubmed/22955828), http://www.uwencode.org/proj/Maurano_et_al_func_var/) defining a set of computationally predicted TFBSs (`hg19.taipale.1e-4.starch`). For a subset of clustered TFs (AP-2-2, AP-2, CTCF_Core-2, E2F-2, EBF1, Ebox-CACCTG, Ebox, ESR1, ETS, MAFK, MEF2A-2, MEF2A, MYC-MAX, PAX5-2, RUNX2, RUNX-AML, STAF-2, TCF-LEF, YY1), we retained only predicted TFBSs that overlap with ENCODE ChIP-seq peaks (TfbsClusteredV3 set downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/). Details on download and basic processing of these data sets are available in `data/Maurano_et_al_func_var/README` and `data/ENCODE_TfbsClusteredV3/README`.

Genomic coordinates of transcription start sites, transcription end sites, start codons, and splice donor and acceptor sites were obtained from Ensembl Build version 75 (ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz, details are available in `data/Ensembl_v75/README`). 

Filtering active CTCF sites
---------------------------

CTCF sites first included clustered FIMO binding site predictions (described above). This set was intersected with ENCODE ChIP-seq peaks (TfbsClusteredV3, described above), and then further intersected with a set of CTCF binding sites experimentally observed to be active across 19 tissues (Wang et al., 2012; details on the download of this data set are available in `data/Wang_et_al_CTCF/README`), to produce three increasingly stringent sets. For each CTCF site, distances between each of 20 flanking nucleosomes  (10 upstream and 10 downstream) were calculated. 

The mean S-WPS and L-WPS at each position relative to the center of the CTCF binding motifs were calculated within bins defined by spacing between -1 and +1 nucleosomes (>160 bp, 161-200 bp, 201-230 bp, 231-270 bp, 271-420 bp, 421-460 bp, and >460 bp). In figures, CTCF binding sites are shifted such that the zero coordinate on the x-axis is the center of its 52 bp binding footprint ([Ong and Corces, 2014](http://www.ncbi.nlm.nih.gov/pubmed/24614316)). 

Overlaying WPS at aligned genomic features
------------------------------------------

WPS values for a specific data set (in most cases CH01) and the corresponding simulation were extracted for each position in a 5 kb window around the start coordinate of each TFBS, and was aggregated within each TF cluster. The mean WPS of the first and last 500 bp (which is predominantly flat and represents a mean offset) of the 5 kb window was subtracted from the original WPS at each position. For L-WPS only, a sliding window mean is calculated using a 200 bp window and subtracted from the original signal. Finally, the corrected WPS profile for the simulation is subtracted from the corrected WPS profile of the data set to correct for signal that is a product of fragment length and ligation bias. This final profile is plotted and termed the Adjusted WPS. 

We provide Python and R scripts for extracting and plotting WPS overlays from block-gzip compressed WIG files. Such files are generated by `extractReadStartsFromBAM2Wig.py` or `normalizeWPSwigs.py` and need to be tabix-indexed to retrieve specific genomic regions using the tabix routines of [pysam](http://pysam.readthedocs.org/en/latest/api.html).

###Motif-based TF binding site predictions

We outline how S-WPS and L-WPS plots were generated from clustered FIMO (motif-based) sites in a README in `WPS_overlays/Maurano_et_al_TFclusters`.

###Motif-based TF binding sites overlapped with ENCODE ChIP-seq peaks

Please find a README in `WPS_overlays/Maurano_et_al_plus_ChIP` which outlines how S-WPS and L-WPS plots were generated for a subset of clustered FIMO (motif-based) sites overlapped with ENCODE ChIP-seq peaks of AP-2-2, AP-2, CTCF_Core-2, E2F-2, EBF1, Ebox-CACCTG, Ebox, ESR1, ETS, MAFK, MEF2A-2, MEF2A, MYC-MAX, PAX5-2, RUNX2, RUNX-AML, STAF-2, TCF-LEF, and YY1.

###Active CTCF from ChIP-seq in 19 cell-lines

`WPS_overlays/CTCF_19CellLines` has a README file which outlines how S-WPS and L-WPS plots were generated for this subset of CTCF sites.

###Adjusted WPS around genic features 

A README describing how the adjusted WPS was extracted and overlayed around TSS, TSE, splice acceptor/donor, and start/stop codon can be found in `WPS_overlays/genicFeatures`.

###Adjusted WPS around genic features in five NB-4 expression bins

The file `WPS_overlays/genicFeatures_by_expression/README` describes how adjusted WPS was extracted and overlayed around TSS, TSE, splice acceptor/donor, and start/stop codon for the five NB-4 expression bins and how the genes in each bin were defined.

Gene expression analysis
------------------------

###Human Protein Atlas

FPKM gene expression (GE) values measured for 20,344 Ensembl gene identifiers in 44 human cell lines and 32 primary tissues by the Human Protein Atlas ([Uhlén et al., 2015](http://www.ncbi.nlm.nih.gov/pubmed/25613900)) was downloaded from http://www.proteinatlas.org/download/rna.csv.zip. Genes with 3 or more non-zero expression values were retained (n=19,378 genes). The GE data set is provided with one decimal precision for the FPKM values. Thus, a zero GE value (0.0) indicates expression in the interval [0, 0.05) Unless otherwise noted, we set the minimum GE value to 0.04 FPKM before log2-transformation..

###Fast Fourier transformation (FFT) and smoothing of trajectories

We use parameters to smooth (3 bp Daniell smoother; moving average giving half weight to the end values) and de-trend the data (i.e. subtract the mean of the series and remove a linear trend). A recursive time series filter implemented in R was used to remove high frequency variation from trajectories. 24 filter frequencies (1/seq(5,100,4)) were used, and the first 24 values of the trajectory were taken as init values. The 24-value shift in the resulting trajectories was corrected by repeating the last 24 values of the trajectory.

###FFT intensity correlation with expression

L-WPS was used to calculate periodograms of genomic regions using Fast Fourier Transform (FFT, spec.pgram in R) with frequencies between 1/500 and 1/100 bases. Intensity values for the 120-280 bp frequency range were determined from smooth FFT periodograms. S-shaped Pearson correlation between GE values and FFT intensities was observed around the major inter-nucleosome distance peak, along with a pronounced negative correlation in the 193-199 bp frequency range. The mean intensity in this frequency range was correlated with the average intensity with log2-transformed GE values for downstream analysis. 

###Running the analysis

```bash

zcat transcriptAnno-GRCh37.75.tsv.gz | awk 'BEGIN{ FS="\t"; OFS="\t" }{ if ($5 == "+") { print $1,$2,$3-10000,$3,$5 } else { print $1,$2,$4,$4+10000,$5 } }' > transcriptAnno-GRCh37.75.upstream.tsv                                                   
zcat transcriptAnno-GRCh37.75.tsv.gz | awk 'BEGIN{ FS="\t"; OFS="\t" }{ if ($5 == "+") { print $1,$2,$4,$4+10000,$5 } else { print $1,$2,$3-10000,$3,$5 } }' > transcriptAnno-GRCh37.75.downstream.tsv                                                 
zcat transcriptAnno-GRCh37.75.tsv.gz | awk 'BEGIN{ FS="\t"; OFS="\t" }{ if ($5 == "+") { print $1,$2,$3-1,$3-1+10000,$5 } else { print $1,$2,$4-1-10000,$4-1,$5 } }' > transcriptAnno-GRCh37.75.body.tsv

# Extract counts for a sample
mkdir -p body/$SAMPLE
./extractReadStartsFromBAM_Region_WPS.py --minInsert=120 --maxInsert=180 -i transcriptAnno-GRCh37.75.body.tsv -o 'body/$SAMPLE/block_%s.tsv.gz' BAMFILE.bam

# Run FFT & convert to summary tbl
mkdir -p /tmp/body/$SAMPLE/fft
( cd body/$SAMPLE/counts; ls block_*.tsv.gz ) | xargs -n 500 Rscript fft_path.R body/$SAMPLE/counts /tmp/body/$SAMPLE/fft
mkdir -p body/fft_summaries/
./convert_files.py -a transcriptAnno-GRCh37.75.body.tsv -t /tmp/ -r  -p body -i $SAMPLE
rm -fR /tmp/body/$SAMPLE/fft

# Correlate the intensities with the expression data and generate PDFs in R
R --vanilla --quiet < plots.R
```

Samtools version used
---------------------

Please note that samtools binary is included with these scripts. Among other things, this samtools binary allows filtering reads based on insert size/read length. This is an early version of the samtools branch released on https://github.com/mpieva/samtools-patched. We have not tested newer versions, but assume that they will work with our scripts too.
