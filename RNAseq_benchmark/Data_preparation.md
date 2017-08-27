# Preparing simulated datasets for RNA benchmark
Several scripts from RNAseqEval project were used for our benchmark of RNA mappers on third generation sequencing reads. Mapers were tested on a set of real and synthetic datasets. To simulate synthetic datasets, we used PBSIM version 1.0.3, downloaded from https://code.google.com/archive/p/pbsim/.

Synthetic datasets were created from the following organisms:
-	Saccharomyces cerevisiae S288 (baker’s yeast)
-	Drosophila melanogaster r6 (fruit fly)
-	Homo Sapiens GRCh38.p7 (human)
Reference genomes for all organisms were downloaded from http://www.ncbi.nlm.nih.gov.

PBSIM is intended to be used as a genomic reads simulator, taking as input a reference sequence and a set of simulation parameters (e.g., coverage, read length, error profile). To generate RNA-seq reads, PBSIM was applied to a set of transcripts generated from a particular genome using the gene annotations downloaded from https://genome.ucsc.edu/cgi-bin/hgTables. To make the datasets as realistic as possible, real datasets were analyzed and used to determine simulation parameters. Real gene expression datasets were used to select a set of transcripts for simulation (downloaded from http://bowtie-bio.sourceforge.net/recount/; core (human), nagalakshmi (yeast) and modencodefly (fruit fly) datasets were used).

## Simulated data preparation
Simulated datasets were generated using the following workflow:
1.	Analyze real datasets to determine error profiles.
2.	Filter annotations (keep only primary assembly information) and unify chromosome names.
3.	Separate annotations for genes with one isoform and genes with alternative splicing, keeping up to 3 isoforms randomly for each gene with alternative splicing.
4.	Generate a transcriptome from processed annotations and a reference genome.
5.	Analyze gene expression data and determine gene coverage histogram (see Figure).
6.	Approximate gene coverage histogram with 3 points to determine coverage and number of genes in simulated dataset (see Figure). Scale coverages proportionally down to make a smaller dataset, more suitable for testing.
7.	Extract 6 subsets of sequences from generated transcriptome, 3 for genes with single splicing and 3 for genes with alternative splicing. Each set contains a number of transcripts corresponding to the number of genes from a previous step.
8.	Using PBSIM, simulate reads on each generated subset of transcriptome, using coverages determined in step 6 and error profiles determined in step 1.
9.	Combine generated reads into a single generated dataset.

### 1. Error profile of real datasets
Real datasets were analyzed to determine error profiles. PacBio datasets were simulated using PacBio ROI isoseq error profile. Due to the lack of ONT MinION RNA reads, ONT MinION dataset was simulated using DNA error profile. Error profiles were determined by mapping the reads to the reference using GraphMap (https://github.com/isovic/graphmap), and by running errorrates.py script from https://github.com/isovic/samscripts. The script take SAM file with alignments and a reference as input, and calculated error rates by examining CIGAR strings for alignments and by comparing the corresponding bases in the read to the ones in the reference. Since PacBio error profile was determined from RNA reads, the reference was a transcriptome (_S.Cerevisiae, D. Melanogaster_ and _H. Sapiens_), while for ONT MinION reads, the reference was a genome (_D. Melanogaster_).

Usage example for errorrates.py:

    errorrates.py base d_melanogaster_transcriptome.fa pacBio_ROI_dataset.sam

### 2. Fintering annotations any unifying chromosome names



### 3. Separating annotations into groups

For simplicity, we rounded the coverage and number of genes from each transcriptome subset. For example, the Table below shows the numbers used to generate dataset 2 (D. melanogaster). The annotation includes roughly 23,000 genes with a single isoform and 3,000 genes with alternative splicing. Rounding up the ratio, we have decided to simulate 1/10 genes with alternative splicing and 9/10 genes without. We considered that each gene undergoing alternative splicing gave rise to three different isoforms with equal expression.

For simulation of PacBio reads, PBSIM parameters (read length, error probability by type, etc) were set to match those of dataset 5 containing reads of insert (see Supplement table 1).
For simulation of MinION ONT reads, PBSIM parameters (read length, error probability by type etc.) were set to match those for MinION reads from a R9 chemistry dataset obtained from the Loman lab website (http://lab.loman.net/2016/07/30/nanopore-r9-data-release). Only 2d reads statistics were used.
