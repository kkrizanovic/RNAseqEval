#! /usr/bin/python

import sys, os
import setup_RNAseqEval

from datetime import datetime

# To enable importing from samscripts submodulew
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import utility_sam
import Annotation_formats

from fastqparser import read_fastq
from report import EvalReport, ReportType

def cleanup():
    pass

DEFAULT_OPTIONS = {}


# TODO:
def eval_mapping_annotations(ref_file, sam_file, options):

    sys.stdout.write('\n')
    sys.stdout.write('\n(%s) Evaluating mapping with annotations:' % datetime.now().time().isoformat())
    # Reading annotation file
    sys.stdout.write('\n(%s) Loading annotations ... ' % datetime.now().time().isoformat())
    annotations = Annotation_formats.Load_Annotation_From_File(ref_file)

    sys.stdout.write('\n(%s) Sorting annotations according to position ... ' % datetime.now().time().isoformat())
    annotations.sort(reverse=False, key=lambda annotation: annotation.start)

    # Loading SAM file into hash
    # Keeping only SAM lines with regular CIGAR string, and sorting them according to position
    sys.stdout.write('\n(%s) Loading SAM file with mappings ... ' % datetime.now().time().isoformat())
    qnames_with_multiple_alignments = {}
    [sam_hash, sam_hash_num_lines, sam_hash_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_file, qnames_with_multiple_alignments)

    sys.stdout.write('\n(%s) Reorganizing SAM lines ... ' % datetime.now().time().isoformat())
    tmpsamlines = [sam_hash[qname][0] for qname in sam_hash.iterkeys()]
    samlines = [samline for samline in tmpsamlines if samline.cigar <> '*']

    sys.stdout.write('\n(%s) Sorting SAM lines ... ' % datetime.now().time().isoformat())
    samlines.sort(key = lambda samline: samline.pos)

    report = EvalReport(ReportType.ANNOTATION_REPORT)

    sys.stdout.write('\n(%s) Analyzing SAM file ... ' % datetime.now().time().isoformat())
    for lines in sam_hash.itervalues():
        if len(lines) > 1:
            report.num_multi_alignments += 1
    report.num_alignments = sam_hash_num_lines
    report.num_unique_alignments = sam_hash_num_unique_lines
    report.num_non_alignments = report.num_alignments - len(samlines)

    sys.stdout.write('\n(%s) Analyzing annotations ... ' % datetime.now().time().isoformat())
    # Analyzing annotations
    report.reflength = 0
    report.num_genes = len(annotations)
    report.num_exons = 0
    sumGeneLength = 0.0
    sumExonLength = 0.0
    for annotation in annotations:
        report.reflength += annotation.end - annotation.start
        report.num_exons += len(annotation.items)
        glength = annotation.getLength()
        if glength < report.min_gene_length or report.min_gene_length == 0:
            report.min_gene_length = glength
        if glength > report.max_gene_length or report.max_gene_length == 0:
            report.max_gene_length = glength
        sumGeneLength += glength
        for item in annotation.items:
            elength = item.getLength()
            if elength < report.min_exon_length or report.min_exon_length == 0:
                report.min_exon_length = elength
            if elength > report.max_exon_length or report.max_exon_length == 0:
                report.max_exon_length = elength
            sumExonLength += elength
    report.avg_gene_length = sumGeneLength / report.num_genes
    report.avg_exon_length = sumExonLength / report.num_exons

    numq = 0
    sumq = 0.0

    # import pdb
    # pdb.set_trace()

    sys.stdout.write('\n(%s) Analyzing mappings ... ' % datetime.now().time().isoformat())
    # Looking at SAM lines to estimate general mapping quality
    # TODO: This used to take a long time, but I managed to speed it up
    #       should be looked at a bit more to see if additional improvements could be made.

    # Looking at expressed genes, ones that overlap with at least one read in SAM file
    # Storing them in a dictionary together with a number of hits for each exon in the gene
    expressed_genes = {}
    report.num_near_miss_alignments = 0

    for samline in samlines:
        # KK: Ovo sada ne moguizracunati ovdje, trebam premjestiti
        # if len(samline) > 1:
        #    report.num_multi_alignments += 1
        quality = samline.chosen_quality
        if quality > 0:
            report.num_good_quality += 1
            if report.max_mapping_quality == 0 or report.max_mapping_quality < quality:
                report.max_mapping_quality = quality
            if report.min_mapping_quality == 0 or report.min_mapping_quality > quality:
                report.min_mapping_quality = quality
            numq += 1
            sumq += quality
        else:
            report.num_zero_quality += 1

        cigar = samline.cigar
        operations = samline.SplitCigar()
        split = False
        for op in operations[1:-1]:
            if op[1] == 'N':
                split = True
        if split:
            report.num_split_alignments += 1

        startpos = samline.pos
        length = samline.CalcReferenceLengthFromCigar()
        endpos = startpos + length

        hit = False
        exonHit = False
        for annotation in annotations:
            if annotation.overlapsGene(startpos, endpos):
                hit = True
                # Initializing a list of counters for a gene
                # Each gene has one global counted (index 0), and one counter for each exon
                expressed_genes[annotation.genename] = [0 for i in xrange(len(annotation.items) + 1)]
                if annotation.genename in expressed_genes.keys():
                    expressed_genes[annotation.genename][0] += 1
                else:
                    expressed_genes[annotation.genename][0] = 1
                if annotation.insideGene(startpos, endpos):
                    partial = False
                else:
                    partial = True
                exonhit = False
                item_idx = 0
                for item in annotation.items:
                    item_idx += 1
                    if item.overlapsItem(startpos, endpos):
                        expressed_genes[annotation.genename][item_idx] += 1
                        exonHit = True

        if hit and not partial:
            report.num_hit_alignments += 1
        elif hit and partial:
            report.num_partial_alignments += 1
        else:
            report.num_missed_alignments += 1

        if hit and not exonHit:
            report.num_near_miss_alignments += 1

    # How many genes were covered by alignments
    report.num_genes_covered = len(expressed_genes)
    report.num_exons_covered = 0
    for genecnt in expressed_genes.itervalues():
        for cnt in genecnt[1:]:
            if cnt > 0:
                report.num_exons_covered += 1

    # TODO: calculate coverage of partial alignments
    #       work with split alignments (ignore Ns)
    #       expand the same logic to exons instead of comlpete genes

    if numq > 0:
        report.avg_mapping_quality = sumq / numq

    sys.stdout.write('\n(%s) Done!' % datetime.now().time().isoformat())

    return report


def eval_mapping_fasta(ref_file, sam_file, options = DEFAULT_OPTIONS):

    # Reading FASTA reference
    [headers, seqs, quals] = read_fastq(ref_file)
    reference_seq = seqs[0]

    # Loading SAM file into hash
    qnames_with_multiple_alignments = {}
    [sam_hash, sam_hash_num_lines, sam_hash_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_file, qnames_with_multiple_alignments)

    report = EvalReport(ReportType.FASTA_REPORT)

    report.reflength = len(reference_seq)
    report.num_alingments = sam_hash_num_lines
    report.num_unique_alignments = sam_hash_num_unique_lines

    numq = 0
    sumq = 0.0

    # Looking at SAM lines to estimate general mapping quality
    for qname in sam_hash.iterkeys():
        samline = sam_hash[qname]
        if len(samline) > 1:
            report.num_multi_alignments += 1
        quality = samline[0].chosen_quality
        if quality > 0:
            report.num_good_quality += 1
            if report.max_mapping_quality == 0 or report.max_mapping_quality < quality:
                report.max_mapping_quality = quality
            if report.min_mapping_quality == 0 or report.min_mapping_quality > quality:
                report.min_mapping_quality = quality
            numq += 1
            sumq += quality
        else:
            report.num_zero_quality += 1

        cigar = samline[0].cigar
        operations = samline[0].SplitCigar()
        split = False
        for op in operations[1:-1]:
            if op[1] == 'N':
                split = True
        if split:
            report.num_split_alignments += 1

        cigar = samline[0].CalcExtendedCIGAR(reference_seq)
        pos = samline[0].pos
        quals = samline[0].qual

    if numq > 0:
        report.avg_mapping_quality = sumq / numq

    return report


def eval_mapping(input_type, ref_file, sam_file, options = DEFAULT_OPTIONS):

    if input_type == '-r':
        report = eval_mapping_fasta(ref_file, sam_file, options)
    elif input_type == '-a':
        report = eval_mapping_annotations(ref_file, sam_file, options)
    else:
        sys.stdout.write('\nInvalid reference type! Should be either -r for fasta or -a for annotation formats!\n')
        exit(1)

    return report


def verbose_usage_and_exit():
    sys.stderr.write('RNAseqEval - A tool for evaulating RNAseq results.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tsetup\n')
    sys.stderr.write('\t\tcleanup\n')
    sys.stderr.write('\t\teval-mapping\n')
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'setup'):
        if (len(sys.argv) != 2):
            sys.stderr.write('Setup the folder structures and install necessary tools.\n')
            sys.stderr.write('Requires no additional parameters to run.\n')
            sys.stderr.write('\n')
            exit(1)

        setup_RNAseqEval.setup_all()

    elif (mode == 'cleanup'):
        if (len(sys.argv) != 2):
            sys.stderr.write('Cleans up intermediate files.\n')
            sys.stderr.write('Requires no additional parameters to run.\n')
            sys.stderr.write('\n')
            exit(1)

        cleanup()

    elif (mode == 'eval-mapping'):
        params = {}
        if (len(sys.argv) < 5):
            sys.stderr.write('Evaluates RNAseq mapping from a SAM file.\n')
            sys.stderr.write('Can use either annotated or unannotaded reference.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s -r <reference FASTA file> <input SAM file> options\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('%s %s -a <reference annotation> (GFF/GTF/BED) file> <input SAM file> options\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('options:"\n')
            sys.stderr.write('-o <option1> - placeholder for option descriptions"\n')
            sys.stderr.write('\n')
            sys.stderr.write('\n')
            exit(1)

        input_type = sys.argv[2]
        ref_file = sys.argv[3]
        sam_file = sys.argv[4]

        # placeholder for options
        options = sys.argv[5:]

        report = eval_mapping(input_type, ref_file, sam_file)

        sys.stdout.write(report.toString())

    else:
        print 'Invalid option!'
        print 'Invalid option!'
