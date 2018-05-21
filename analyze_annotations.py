#! /usr/bin/python

# This file contains mathods for manipulating fasta and annotation files
# To make them appropriate for testing RNAseq mapping

# All headers will be changed to chr[ID]

import sys, os
import random

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import utility_sam
import Annotation_formats

from fastqparser import read_fastq


def analyze(annotations_file):

    filename, file_extension = os.path.splitext(annotations_file)

    if file_extension.lower() in ['.gtf', '.gff']:
        filetype = 'GTF'
    elif file_extension.lower() in ['.bed']:
        filetype = 'BED'
    else:
        raise Exception('Invalid annotation file type: %s' % file_extension)

    # Reading annotation file
    # annotations = Annotation_formats.Load_Annotation_From_File(annotations_file, check_duplicates = True)
    annotations = Annotation_formats.Load_Annotation_From_File(annotations_file)

    # for annotation in annotations:
    #     if len(annotation.items) > 1 and annotation.genename[0] == 'Q':
    #         import pdb
    #         pdb.set_trace()

    # Analyzing annotations to discover alternate splicings
    # Grouping annotations which overlap and are on the same strand
    annotation_groups = {}
    group_found = True
    gene_start = gene_end = trcnt = iden = 0

    for annotation in annotations:
        group_found = False
        for idgroup, group in annotation_groups.iteritems():
            gene_start = group[0]
            gene_end = group[1]
            trcnt = group[2]
            iden = idgroup
            if annotation.overlapsGene(gene_start, gene_end):
                group_found = True
                break

        if group_found:
            if annotation.start < gene_start:
                gene_start = annotation.start
            if annotation.end > gene_end:
                gene_end = annotation.end
            trcnt += 1
            annotation_groups[iden] = (gene_start, gene_end, trcnt)
        else:
            iden = annotation.start
            annotation_groups[iden] = (annotation.start, annotation.end, 1)


    new = True
    new_groups = {}
    groupid = group_start = group_end = 0
    for iden, group in sorted(annotation_groups.iteritems(), key=lambda(k,v):v[0]):
        # group = annotation_groups[iden]

        if new:
            groupid = group[0]
            group_start = group[0]
            group_end = group[1]
            trcnt = group[2]
            new = False
        else:
            # If overlaps with the current group join it
            if  not (group[0] > group_end or group[1] < group_start):
                if group[0] < group_start:
                    group_start = group[0]
                if group[1] > group_end:
                    group_end = group[1]
                trcnt += group[2]

            # And if it doesnt overlap, add old group to the new group dictionary
            # And start a new group
            else:
                new_groups[groupid] = (group_start, group_end, trcnt)
                groupid = group[0]
                group_start = group[0]
                group_end = group[1]
                trcnt = group[2]

    # Add last group to thenew  dictionary
    new_groups[groupid] = (group_start, group_end, trcnt)


    sys.stderr.write("\nWritting annotation groups (%d)\n" % len(new_groups))
    sys.stdout.write("ID\tSTART\tEND\tTRCNT\n")
    for idgroup in sorted(new_groups.iterkeys()):
        group = new_groups[idgroup]
        sys.stdout.write("%d\t%d\t%d\t%d\n" % (idgroup, group[0], group[1], group[2]))

        

def verbose_usage_and_exit():
    sys.stderr.write('This script is used to analyze gene annotations.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode] [filename]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tanalyze - analyze annotations\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 3):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'analyze'):
        annotations_file = sys.argv[2]
        analyze(annotations_file)

    else:
        print 'Invalid mode: %s!' % mode
