# RNAseqEval.py
Run RNAseqEval.py for general evaulation of mappings in sam file against reference and optionally annotations. This script is intended to evaluate real dataset mapping. Run RNAseqEval.py without any arguments to print options.

Usage:
     
    RNAseqEval.py eval-mapping <reference FASTA file> <input SAM file> options

## Usage modes
RNAseqEval.py script can be used in three differents modes, determined by the first argument. Each mode requires different parameters and allowes different options.

### eval-mapping
Used in eval-mapping mode, RNAseqEval.py script is used to evaluate RNAseq mappings against known FASTA reference and annotations. Annotations can be omitted, but in that case the script will provide only basic output.

Usage:

    RNAseqEval.py eval-mapping <reference FASTA file> <input SAM file> options
    
Allowed options:

    -a <file> : a reference annotation (GFF/GTF/BED) file
    -o (--output) <file> : output file to which the report will be written

### eval-annotations
Used in eval-annotations mode, RNAseqEval.py script will print out basic information on a annotations file.

Usage:

    RNAseqEval.py eval-annotations <annotations file> options

Allowed options:

    -o (--output) <file> : output file to which the report will be written

### eval-maplength
Used in eval-maplength mode, RNAseqEval script will return mapped percentage for each read

Usage:

    RNAseqEval.py eval-maplength <input SAM file> options

Options:

    -o (--output) <file> : output file to which the report will be written

## Output
Depending on the usage mode, RNAseqEval.py script will display various information about input files and the results of the analysis.

General information on FASTA reference and mapping SAM file:

     Reference length
     Number of chromosomes
     List of chromosomes
     Number of alignments in SAM file (total / unique)
     Alignments with / without CIGAR string
     Mapping quality without zeroes (avg / min / max)
     Alignments with mapping quality (>0 / =0)
     Number of matches / mismatches / inserts / deletes
     Percentage of matches / mismatches / inserts / deletes
