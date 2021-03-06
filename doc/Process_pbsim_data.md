# Process_pbsim_data.py
Run Process_pbsim_data.py to evaluate mappings of simulated data generated by PBSIM simulator. Since PBSIM is a DNA read simulator, to generate RNA reads the simulation should be run on a transcriptome. To evalute mapping quality, the script will therefore require SAM file containing mapping results, annotations from which the transcriptome was generated and the simulation data generated by PBSIM.

Usage example:
 
    Process_pbsim_data.py process simulation_root_folder mappings.sam annotations.bed

## Evaluation method
Process_pbsim_data.py evaluates the quality of mapping by comparing the mapping file in SAM format to the simulation data generated by PBSIM. Simulation data is generated by applying PBSIM to a set of transcripts which is in turn generated from a genome and a set of annotations. To evaluate the mapping results, the script will therefore need both, reference genome and a set of annotations.

Aside from simulation parameters, PBSIM will take as input a fasta file containing one or more references. PBSIM will use all references to generate simulated reads. For each reference it will generate 3 files:
1. REF file containing the reference
2. FASTQ file containing simulated reads
3. MAF file containing the information from which part of the reference was each read generated

__IMPORTANT:__ The Process_pbsim_data.py script requires a default PBSIM naming convention to be used, generating filenames begining with "sd\_" and followed by a 4 digit number. While running PBSIM simulations DO NOT use option --prefix.

The evaluation process evaluates one read (alignment) at a time, and for each read (alignment) goes through the following steps:
1. From read header determine the files containing the information on read reference (transcript) and read origin (position on reference).
2. From PBSIM generated reference file determine the annotation used to simulate the read
3. From PBSIM generated MAF file, determine the position and length of the read on the reference
4. By mapping the read position from the reference (transcript) to the genome (using the correct annotation), determine the correct alignment of the read to the genome
5. Compare the calculated correct alignment to the alignment being tested

## Evaluating a complex simulated dataset
The script is written so that it can work for datasets constructed from multiple simulations simulations. This way each simulation can use different coverage and be run on a different set of transcripts, allowing more complex datasets. However, all simulations must be based on the same reference genome and the same total set of annotations (annotations with the same name used in multiple simulations must be defined identically).

Evaluating mapping quality for a complex dataset, constructed from multiple simulations, requires additional steps.
- Reads simulated within a single simulation need to have the same FASTA header prefix, unique for each simulation.
- All PBSIM output folders (for different simulations) must be within a single simulation_root_folder
- A dictionary simFolderDict (inside Process_pbsim_data.py script) must contain the connection between read header prefix belonging to a simulation and the corresponding pbsim output folder placed within the simulation_root_folder. An example of simFolderDict definition used for our RNA benchmark is given in [RNAseq_benchmark.py](RNAseq_benchmark.py).

## Evaluation example
Lets suppose that we wish to use three simulations to construct a synthetic dataset. We generate 3 sets of transcripts for the same genome and run a simulation on each of them. PBSIM output folders for each simulation are named group1, group2 and group3 and are placed within a folder named complex_simulation. Header prefixes for each simulation reads are SimG1, SimG2 and SimG3 respectively. 
Dictionary simFolderDict, within the Process_pbsim_data.py script, must be definied as:

    simFolderDict = {'SimG1' : 'group1'
                   , 'SimG2' : 'group2'
                   , 'SimG3' : 'group3'}

The evaluation is run as follows:

    Process_pbsim_data.py process complex_simulation  mappings.sam annotations.bed

A detailed description of the process used to prepare simulated datasets for our RNA benchmark is given in [RNAseq_benchmark/data_preparation.md](RNAseq_benchmark/data_preparation.md).

## Output
Process_pbsim_daty.py scripts generates a report containing various information.

It contains some general information on the input files:

    Original Samlines - Number of original SAM lines (in the input SAM file)
    Usable whole alignments (with valid CIGAR string) - number of alignments having a valid CIGAR string (in the input SAM file)
    Annotations - number of annotations (in the input GFF/GTF/BED file)
    Multiexon genes - number of annotations describing genes with more than one exon (in the input GFF/GTF/BED file)

It contains some general infomation gained by comparing alignments to annotations:

    Number of exon start hits - number of reads correctly aligned to a start of an exon
    Number of exon end hits - number of reads correctly aligned to an end of an exon
    Number of exon start and end hits - number of reads correctly aligned to both a start and an end of an exon
    Number of good whole alignments - number of alignments that overlap at least one exon

It contains some infomation obtained by comparing alignments from the input SAM file to the read origins from a set of MAF files withing the simulation folder:

    MAF: Correct alignment - the number of alignments that correctly map the read to its position of origin (within 5 bases)
    MAF: Hit all parts - the number of alignments that overlap all exons from the read origin, each overlap must be at least 5 bases
    MAF: Hit at least one part - the number of alignments that overlap at least one exon from the read origin, the overlap must be at least 5 bases
    MAF: Equals at least one part - the number alignments that correctly map the read to at least one exon from its position of origin (within 5 bases)
    MAF: Number of split reads - the number of reads that originate from more than one exon
    MAF: Correct alignment, SPLIT read - same as "MAF: Correct alignment", but only for "split reads"
    MAF: Hit all parts, SPLIT read - same as "MAF: Hit all parts", but only for "split reads"
    MAF: Hit at least one part, SPLIT read - same as "MAF: Hit at least one part", but only for "split reads"
    MAF: Equals at least one part, SPLIT read - same as "MAF: Equals at least one part", but only for "split reads"
