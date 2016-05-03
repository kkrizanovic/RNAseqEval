#! /usr/bin/python

import sys, os
import re
import setup_RNAseqEval, paramsparser

from datetime import datetime

# To enable importing from samscripts submodulew
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import utility_sam
import Annotation_formats

from fastqparser import read_fastq
from report import EvalReport, ReportType

# Parameter definitions for paramparser
paramdefs = {'-a' : 1,
             '--version' : 0,
             '-v' : 0}


def cleanup():
    pass


# A helper function that extracts a chromosome name from a fasta header
# Chromosome names should be either chromosome [designation] or chr[designation]
# In the case header represents a mitochondrion, 'chrM' is returned!
def getChromName(header):
    chromname = ''
    # regular expressions for searching for long and short chromosome names
    longre = r'(chromosome )(\w*)'
    shortre = r'(chr)(\w*)'

    if header.find('mitochondrion') > -1 or header.find('chrM') > -1:
        chromname = 'chrM'
    else:
        match1 = re.search(longre, header)
        match2 = re.search(shortre, header)
        if match1:
            designation = match1.group(2)
            chromname = 'chr%s' % designation
        elif match2:
            designation = match2.group(2)
            chromname = 'chr%s' % designation
        else:
            chromname = 'genome'        # In case I can't detect chromosome designation or mitochondrion
                                        # I decide its a genome

    return chromname



def load_and_process_reference(ref_file, paramdict, report):
    # Reading FASTA reference
    [headers, seqs, quals] = read_fastq(ref_file)

    # Analyzing FASTA file

    # Since annotation file and SAM file (mapper output) do not have to have unified sequence names
    # for individual chromosomes, I'm creating a translation dictionary that will help me access the right
    # sequence more quickly. This will be based on the function which will attempt to analyze a string
    # and determine whether it refers to a genome, a particular chromosome or a mitochndia
    # The chromname2seq dictionary will have a inferred name as key, and index of the corresponding
    # sequence and headeer as value
    chromname2seq = {}

    if len(headers) == 1:
        report.reflength = len(seqs[0])
        chromname = getChromName(headers[0])
        report.chromlengths = {chromname : report.reflength}
        chromname2seq[chromname] = 0
    else:
        for i in xrange(len(headers)):
            header, seq = headers[i], seqs[i]
            chromname = getChromName(header)
            if chromname in report.chromlengths:
                raise Exception('\nERROR: Duplicate chromosome name: %s' % chromname)
                # sys.stderr.write('\nERROR: Duplicate chromosome name: %s' % chromname)
                exit()
            else:
                report.chromlengths[chromname] = len(seq[i])
                report.reflength += len(seq[i])
                chromname2seq[chromname] = i

    return [chromname2seq, headers, seqs, quals]


def load_and_process_SAM(sam_file, paramdict, report):
    # Loading SAM file into hash
    # Keeping only SAM lines with regular CIGAR string, and sorting them according to position
    qnames_with_multiple_alignments = {}
    [sam_hash, sam_hash_num_lines, sam_hash_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_file, qnames_with_multiple_alignments)

    # Reorganizing SAM lines
    tmpsamlines = [sam_hash[qname][0] for qname in sam_hash.iterkeys()]
    samlines = [samline for samline in tmpsamlines if samline.cigar <> '*']

    # Sorting SAM lines
    samlines.sort(key = lambda samline: samline.pos)

    # Analyzing SAM file
    for lines in sam_hash.itervalues():
        if len(lines) > 1:
            report.num_multi_alignments += 1
    report.num_alignments = sam_hash_num_lines
    report.num_unique_alignments = sam_hash_num_unique_lines
    report.num_non_alignments = report.num_alignments - len(samlines)

    return samlines



def load_and_process_annotations(annotations_file, paramdict, report):
    # Reading annotation file
    annotations = Annotation_formats.Load_Annotation_From_File(annotations_file)

    # Sorting annotations according to position
    annotations.sort(reverse=False, key=lambda annotation: annotation.start)

    # Analyzing annotations

    # Looking at expressed genes, ones that overlap with at least one read in SAM file
    # Storing them in a dictionary together with a number of hits for each exon in the gene
    expressed_genes = {}

    report.totalGeneLength = 0
    report.num_genes = len(annotations)
    report.num_exons = 0
    sumGeneLength = 0.0
    sumExonLength = 0.0
    for annotation in annotations:
        # Initializing a list of counters for a gene
        # Each gene has one global counted (index 0), and one counter for each exon
        expressed_genes[annotation.genename] = [0 for i in xrange(len(annotation.items) + 1)]

        report.totalGeneLength += annotation.getLength()
        chromname = getChromName(annotation.seqname)
        if chromname in report.chromlengths:
            report.chromlengths[chromname] += annotation.getLength()
        else:
            report.chromlengths[chromname] = annotation.getLength()
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

    return annotations



# TODO: Refactor code, place some code in functions
#       Rewrite analyzing SAM file, detecting multi and split alignments
def eval_mapping_annotations(ref_file, sam_file, annotations_file, paramdict):

    sys.stdout.write('\n')
    sys.stdout.write('\n(%s) START: Evaluating mapping with annotations:' % datetime.now().time().isoformat())

    # Reading FASTA reference
    sys.stdout.write('\n(%s) Loading FASTA reference ... ' % datetime.now().time().isoformat())
    [headers, seqs, quals] = read_fastq(ref_file)

    # Reading annotation file
    sys.stdout.write('\n(%s) Loading annotations ... ' % datetime.now().time().isoformat())
    annotations = Annotation_formats.Load_Annotation_From_File(annotations_file)

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

    # Analyzing FASTA file
    sys.stdout.write('\n(%s) Analyzing FASTA file ... ' % datetime.now().time().isoformat())

    # Since annotation file and SAM file (mapper output) do not have to have unified sequence names
    # for individual chromosomes, I'm creating a translation dictionary that will help me acces the right
    # sequence more quickly. This will be based on the function which will attempt to analyze a string
    # and determine whether it refers to a genome, a particular chromosome or a mitochndia
    # The chromname2seq dictionary will have a inferred name as key, and index of the corresponding
    # sequence and headeer as value
    chromname2seq = {}

    if len(headers) == 1:
        report.reflength = len(seqs[0])
        chromname = getChromName(headers[0])
        report.chromlengths = {chromname : report.reflength}
        chromname2seq[chromname] = 0
    else:
        for i in xrange(len(headers)):
            header, seq = headers[i], seqs[i]
            chromname = getChromName(header)
            if chromname in report.chromlengths:
                sys.stderr.write('\nERROR: Duplicate chromosome name: %s' % chromname)
                exit()
            else:
                report.chromlengths[chromname] = len(seq[i])
                report.reflength += len(seq[i])
                chromname2seq[chromname] = i

    sys.stdout.write('\n(%s) Analyzing SAM file ... ' % datetime.now().time().isoformat())
    for lines in sam_hash.itervalues():
        if len(lines) > 1:
            report.num_multi_alignments += 1
    report.num_alignments = sam_hash_num_lines
    report.num_unique_alignments = sam_hash_num_unique_lines
    report.num_non_alignments = report.num_alignments - len(samlines)

    sys.stdout.write('\n(%s) Analyzing annotations ... ' % datetime.now().time().isoformat())
    # Analyzing annotations

    # Looking at expressed genes, ones that overlap with at least one read in SAM file
    # Storing them in a dictionary together with a number of hits for each exon in the gene
    expressed_genes = {}

    report.totalGeneLength = 0
    report.num_genes = len(annotations)
    report.num_exons = 0
    sumGeneLength = 0.0
    sumExonLength = 0.0
    for annotation in annotations:
        # Initializing a list of counters for a gene
        # Each gene has one global counted (index 0), and one counter for each exon
        expressed_genes[annotation.genename] = [0 for i in xrange(len(annotation.items) + 1)]

        report.totalGeneLength += annotation.getLength()
        chromname = getChromName(annotation.seqname)
        if chromname in report.chromlengths:
            report.chromlengths[chromname] += annotation.getLength()
        else:
            report.chromlengths[chromname] = annotation.getLength()
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

    sys.stdout.write('\n(%s) Analyzing mappings ... ' % datetime.now().time().isoformat())
    # Looking at SAM lines to estimate general mapping quality
    # TODO: This used to take a long time, but I managed to speed it up
    #       should be looked at a bit more to see if additional improvements could be made.

    numMatch = 0
    numMisMatch = 0
    numInsert = 0
    numDelete = 0

    for samline in samlines:
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

        chromname = getChromName(samline.rname)
        if chromname not in chromname2seq:
            raise Exception('\nERROR: Unknown choromosome name in SAM file! (chromname:"%s", samline.rname:"%s")' % (chromname, samline.rname))
        chromidx = chromname2seq[chromname]

        cigar = samline.CalcExtendedCIGAR(seqs[chromidx])
        pos = samline.pos
        quals = samline.qual

        # Using regular expressions to find repeating digit and skipping one character after that
        # Used to separate CIGAR string into individual operations
        pattern = '(\d+)(.)'
        operations = re.findall(pattern, cigar)

        # I'll assume that N's will not appear at the beginning or at the end of CIGAR string
        # so any N will mean that the alignment is split
        split = False
        for op in operations:
            if op[1] == 'N':
                split = True
            elif op[1] in ('M', '='):
                numMatch += int(op[0])
            elif op[1] == 'I':
                numInsert += int(op[0])
            elif op[1] == 'D':
                numDelete += int(op[0])
            elif op[1] =='X':
                numMisMatch += int(op[0])
            elif op[1] in ('S', 'H', 'P'):
                pass
            else:
                sys.stderr.write('\nERROR: Invalid CIGAR string operation (%s)' % op[1])

        if split:
            report.num_split_alignments += 1

        startpos = samline.pos
        length = samline.CalcReferenceLengthFromCigar()
        endpos = startpos + length

        hit = False
        exonHit = False
        exon_cnt = 0        # counting exons spanned by an alignement
        genescovered = []   # For testing, genes covered by an alignment
        for annotation in annotations:
            if annotation.overlapsGene(startpos, endpos):
                genescovered.append(annotation.genename)
                hit = True

                # Updating gene expression
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
                        exon_cnt += 1
                        expressed_genes[annotation.genename][item_idx] += 1
                        exonHit = True
                        if item.insideItem(startpos, endpos):
                            exonPartial = False
                        else:
                            exonPartial = True

                # For split alignments (alignments that span multiple exons)
                # I'm checking whether each part completeley falls inside an exon
                # this only makes sense if an alignment is completely inside a gene
                # NOTE: CIGAR string operations are already created
                badsplit = False
                if split and partial:
                    badsplit = True
                elif split and not partial:
                    t_startpos = startpos
                    for op in operations:
                        if op[1] == 'N':                    # Skipping Ns
                            t_startpos += int(op[0])
                        elif op[1] in ('M', '=', 'D', 'S', 'H', 'P', 'X'):      # NOTE: Inserts are ignored because they do not alter reference position!
                            t_endpos = t_startpos + int(op[0])
                            if not annotation.insideItems(t_startpos, t_endpos):
                                badsplit = True
                                break
                            t_startpos = t_endpos
                        else:
                            pass

        if exon_cnt > 1:
            report.num_multiple_exons += 1

        if len(genescovered) > 1:
            report.num_multi_gene_alignments += 1

        if badsplit:
            report.num_bad_split_alignments += 1

        if hit and not partial:
            report.num_hit_alignments += 1
        elif hit and partial:
            report.num_partial_alignments += 1
        else:
            report.num_missed_alignments += 1

        if exonHit and not exonPartial:
            report.num_exon_hit += 1
        elif exonHit and exonPartial:
            report.num_exon_partial += 1
        else:
            report.num_exon_miss += 1

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
    #       expand the same logic to exons instead of complete genes

    report.num_match = numMatch
    report.num_mismatch = numMisMatch
    report.num_insert = numInsert
    report.num_delete = numDelete

    total = numMatch + numMisMatch + numInsert + numDelete

    if total > 0:
        report.match_percentage = float(report.num_match)/total
        report.mismatch_percentage = float(report.num_mismatch)/total
        report.insert_percentage = float(report.num_insert)/total
        report.delete_percentage = float(report.num_delete)/total

    if numq > 0:
        report.avg_mapping_quality = sumq / numq

    sys.stdout.write('\n(%s) Done!' % datetime.now().time().isoformat())

    return report


def eval_mapping_fasta(ref_file, sam_file, paramdict):

    # Reading FASTA reference
    sys.stdout.write('\n(%s) Loading FASTA reference ... ' % datetime.now().time().isoformat())
    [headers, seqs, quals] = read_fastq(ref_file)

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

    report = EvalReport(ReportType.FASTA_REPORT)

    # Analyzing FASTA file
    sys.stdout.write('\n(%s) Analyzing FASTA file ... ' % datetime.now().time().isoformat())

    # Since annotation file and SAM file (mapper output) do not have to have unified sequence names
    # for individual chromosomes, I'm creating a translation dictionary that will help me acces the right
    # sequence more quickly. This will be based on the function which will attempt to analyze a string
    # and determine whether it refers to a genome, a particular chromosome or a mitochndia
    # The chromname2seq dictionary will have a inferred name as key, and index of the corresponding
    # sequence and headeer as value
    chromname2seq = {}

    if len(headers) == 1:
        report.reflength = len(seqs[0])
        chromname = getChromName(headers[0])
        report.chromlengths = {chromname : report.reflength}
        chromname2seq[chromname] = 0
    else:
        for i in xrange(len(headers)):
            header, seq = headers[i], seqs[i]
            chromname = getChromName(header)
            if chromname in report.chromlengths:
                sys.stderr.write('\nERROR: Duplicate chromosome name: %s' % chromname)
                exit()
            else:
                report.chromlengths[chromname] = len(seq[i])
                report.reflength += len(seq[i])
                chromname2seq[chromname] = i

    sys.stdout.write('\n(%s) Analyzing SAM file ... ' % datetime.now().time().isoformat())
    for lines in sam_hash.itervalues():
        if len(lines) > 1:
            report.num_multi_alignments += 1
    report.num_alignments = sam_hash_num_lines
    report.num_unique_alignments = sam_hash_num_unique_lines
    report.num_non_alignments = report.num_alignments - len(samlines)

    report.num_alingments = sam_hash_num_lines
    report.num_unique_alignments = sam_hash_num_unique_lines

    numq = 0
    sumq = 0.0

    # Analyzing mappings
    sys.stdout.write('\n(%s) Analyzing mappings ... ' % datetime.now().time().isoformat())

    numMatch = 0
    numMisMatch = 0
    numInsert = 0
    numDelete = 0

    # Looking at SAM lines to estimate general mapping quality
    for samline in samlines:
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

        chromname = getChromName(samline.rname)
        if chromname not in chromname2seq:
            raise Exception('\nERROR: Unknown choromosome name in SAM file! (chromname:"%s", samline.rname:"%s")' % (chromname, samline.rname))
        chromidx = chromname2seq[chromname]

        cigar = samline.CalcExtendedCIGAR(seqs[chromidx])
        pos = samline.pos
        quals = samline.qual

        # Using regular expressions to find repeating digit and skipping one character after that
        pattern = '(\d+)(.)'
        operations = re.findall(pattern, cigar)

        # I'll assume that N's will not appear et the beginning or at the end of CIGAR string
        # so any N will mean that the alignment is split
        split = False
        for op in operations:
            if op[1] == 'N':
                split = True
            elif op[1] in ('M', '='):
                numMatch += int(op[0])
            elif op[1] == 'I':
                numInsert += int(op[0])
            elif op[1] == 'D':
                numDelete += int(op[0])
            elif op[1] =='X':
                numMisMatch += int(op[0])
            elif op[1] in ('S', 'H', 'P'):
                pass
            else:
                sys.stderr.write('\nERROR: Invalid CIGAR string operation (%s)' % op[1])
        if split:
            report.num_split_alignments += 1

    report.num_match = numMatch
    report.num_mismatch = numMisMatch
    report.num_insert = numInsert
    report.num_delete = numDelete

    total = numMatch + numMisMatch + numInsert + numDelete

    if total > 0:
        report.match_percentage = float(report.num_match)/total
        report.mismatch_percentage = float(report.num_mismatch)/total
        report.insert_percentage = float(report.num_insert)/total
        report.delete_percentage = float(report.num_delete)/total

    if numq > 0:
        report.avg_mapping_quality = sumq / numq

    return report


def eval_mapping(ref_file, sam_file, paramdict):

    if '-a' in paramdict:
        annotations_file = paramdict['-a'][0]
        report = eval_mapping_annotations(ref_file, sam_file, annotations_file, paramdict)
    else:
        report = eval_mapping_fasta(ref_file, sam_file, paramdict)

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
        if (len(sys.argv) < 4):
            sys.stderr.write('Evaluates RNAseq mapping from a SAM file.\n')
            sys.stderr.write('Can use annotations if provided.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s <reference FASTA file> <input SAM file> options\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('options:"\n')
            sys.stderr.write('-a <file> : a reference annotation (GFF/GTF/BED) file"\n')
            sys.stderr.write('\n')
            sys.stderr.write('\n')
            exit(1)

        ref_file = sys.argv[2]
        sam_file = sys.argv[3]

        pparser = paramsparser.Parser(paramdefs)
        paramdict = pparser.parseCmdArgs(sys.argv[4:])

        ref_file = sys.argv[2]
        sam_file = sys.argv[3]

        report = eval_mapping(ref_file, sam_file, paramdict)

        sys.stdout.write(report.toString())

    else:
        print 'Invalid option!'
        print 'Invalid option!'
