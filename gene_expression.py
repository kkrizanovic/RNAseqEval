#! /usr/bin/python

# This file contains mathods for manipulating gene expression data
# Mainly, it tries to use Biopython to connect gene expression files with annotations

import sys, os
from Bio import Entrez

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))

import utility_sam
from fastqparser import read_fastq


def verbose_usage_and_exit():
    sys.stderr.write('RNAseqEval - A tool for evaulating RNAseq results.\n')
    sys.stderr.write('This script is used to work with gene expression data.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode] [filename]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tsearch - Process S.Cerevisiae genome\n')
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) != 3):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'sc-genome'):
        genome_file = sys.argv[2]
        prepare_sc_genome(genome_file)

    elif (mode == 'sc-annotations'):
        annotations_file = sys.argv[2]
        prepare_sc_annotations(annotations_file)

    elif (mode == 'dm-genome'):
        genome_file = sys.argv[2]
        prepare_dm_genome(genome_file)

    elif (mode == 'dm-annotations'):
        annotations_file = sys.argv[2]
        prepare_dm_annotations(annotations_file)

    elif (mode == 'h-genome'):
        genome_file = sys.argv[2]
        prepare_human_genome(genome_file)

    elif (mode == 'h-annotations'):
        annotations_file = sys.argv[2]
        prepare_human_annotations(annotations_file)

    else:
        print 'Invalid mode: %s!' % mode
