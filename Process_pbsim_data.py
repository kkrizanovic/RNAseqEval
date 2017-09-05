#! /usr/bin/python

import sys, os
import paramsparser

from datetime import datetime

# To enable importing from samscripts submodulew
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import utility_sam
import Annotation_formats
import RNAseqEval
from report import EvalReport, ReportType
from RNAseq_benchmark import benchmark_params

from fastqparser import read_fastq

# Determines whether to check the strand whene analyzing data
# Due to complications in generating simulated RNA reads, this is False
P_CHECK_STRAND = False


# OLD: Predefined dictionaries for analyzing different datasets
# simFolderDict_d1 = {'SimG1' : 'group1'
#                   , 'SimG2' : 'group2'
#                   , 'SimG3' : 'group3'}
#
# simFolderDict_all = {'SimG1' : 'group1'
#                    , 'SimG2' : 'group2'
#                    , 'SimG3' : 'group3'
#                    , 'SimG1AS' : 'group1_AS'
#                    , 'SimG1SS' : 'group1_SS'
#                    , 'SimG2AS' : 'group2_AS'
#                    , 'SimG2SS' : 'group2_SS'
#                    , 'SimG3AS' : 'group3_AS'
#                    , 'SimG3SS' : 'group3_SS'}


# A dictionary connecting fasta/fastq header prefix with the folder with pbsim generated data
# Containing information for reads with each prefix
# This is used because data is simulated using several pbsim runs to get different
# coverages for different sets of references (in this case transcripts)
# NOTE: this should be changed for different simulations
simFolderDict = benchmark_params.simFolderDict


def interval_equals(interval1, interval2, allowed_inacc = Annotation_formats.DEFAULT_ALLOWED_INACCURACY):
    if interval1[0] < interval2[0] - allowed_inacc:
        return False
    if interval1[0] > interval2[0] + allowed_inacc:
        return False
    if interval1[1] < interval2[1] - allowed_inacc:
        return False
    if interval1[1] > interval2[1] + allowed_inacc:
        return False

    return True

def interval_overlaps(interval1, interval2, allowed_inacc = Annotation_formats.DEFAULT_ALLOWED_INACCURACY):

    if (interval1[1] <= interval2[0] + allowed_inacc) or (interval1[0] >= interval2[1] - allowed_inacc):
        return False
    else:
        return True


def processData(datafolder, resultfile, annotationfile):

    # Loading results SAM file
    report = EvalReport(ReportType.FASTA_REPORT)    # not needed
    paramdict = {}

    sys.stderr.write('\n(%s) Loading and processing SAM file with mappings ... ' % datetime.now().time().isoformat())
    all_sam_lines = RNAseqEval.load_and_process_SAM(resultfile, paramdict, report, BBMapFormat = True)


    # Reading annotation file
    annotations = Annotation_formats.Load_Annotation_From_File(annotationfile)

    s_num_multiexon_genes = 0

    # Hashing annotations according to name
    annotation_dict = {}
    for annotation in annotations:
        if annotation.genename in annotation_dict:
            sys.stderr.write('\nWARNING: anotation with name %s already in the dictionary!' % annotation.genename)
        else:
            annotation_dict[annotation.genename] = annotation
        if len(annotation.items) > 1:
            s_num_multiexon_genes += 1


    # Statistical information for evaluating the qualitiy of mapping
    s_gene_hits = 0
    s_gene_misses = 0
    s_whole_alignment_hits = 0
    s_whole_alignment_misses = 0
    s_partial_alignment_hits = 0
    s_partial_alignment_misses = 0
    s_num_start_hits = 0
    s_num_end_hits = 0
    s_num_start_end_hits = 0

    s_num_fw_strand = 0
    s_num_rv_strand = 0

    s_num_split_alignment = 0
    s_num_oversplit_alignment = 0       # Alignments that have more parts than exons

    s_num_good_alignments = 0

    s_num_badchrom_alignments = 0

    s_maf_suspicious_alignments = 0
    s_maf_bad_alignments = 0
    s_maf_good_alignments = 0

    s_maf_split_reads = 0
    s_maf_good_split_alignments = 0
    s_maf_bad_split_alignments = 0

    s_maf_hit_all_parts = 0
    s_maf_hit_one_part = 0
    s_maf_eq_one_part = 0
    s_maf_multihit_parts = 0

    s_maf_split_hit_all_parts = 0
    s_maf_split_hit_one_part = 0
    s_maf_split_eq_one_part = 0

    s_num_potential_bad_strand = 0

    # allowed_inacc = Annotation_formats.DEFAULT_ALLOWED_INACCURACY       # Allowing some shift in positions
    # Setting allowed inaccuracy
    allowed_inacc = 5

    # All samlines in a list should have the same query name
    for samline_list in all_sam_lines:
        qname = samline_list[0].qname

        isSplitAlignment = False
        if len(samline_list) > 1:
            s_num_split_alignment += 1
            isSplitAlignment = True

        # Checking the SAM file if all samlines in a list have the same qname
        for samline in samline_list[1:]:
            if samline.qname != qname:
                sys.stderr.write('\nWARNING: two samlines in the same list with different query names (%s/%s)' % (qname, samline.qname))

        # Look for the first underscore in query name
        # Everything before that is the simulation folder name
        # Everything after that is simulated query name
        pos = qname.find('_')
        if pos < 0:
            raise Exception('Invalid query name in results file (%s)!' % qname)

        simFolderKey = qname[:pos]
        if simFolderKey not in simFolderDict:
            # import pdb
            # pdb.set_trace()
            raise Exception('Bad simulation folder short name (%s)!' % simFolderKey)
        simFolder = simFolderDict[simFolderKey]
        simQName = qname[pos+1:]

        # Due to error in data preparation, have to make some extra processing
        if simQName[:6] == 'SimG2_':
            simQName = simQName[6:]


#        if simFolderKey == 'SimG1':
#            simFileSuffix = 'g1'
#        elif simFolderKey == 'SimG2':
#            simFileSuffix = 'g2'
#        elif simFolderKey == 'SimG3':
#            simFileSuffix = 'g3'
#        else:
#            simFileSuffix = 'sd'

        simFileSuffix = 'sd'


        pos = simQName.find('_')
        pos2 = simQName.find('_part')
        if pos < 0:
            raise Exception('Invalid simulated query name in results file (%s)!' % simQName)

        simQLetter = simQName[0]       # Should always be S

        # BBMap separates a query into smaller parts he can manage
        # Extends query with '_part_#', which has to be ignored
        if pos2 <> -1:
            simQName = simQName[:pos2]

        simRefNumber = int(simQName[1:pos])
        simQNumber = int(simQName[pos+1:])
        simFileName = simFileSuffix + '_%04d' % simRefNumber
        simRefFileName = simFileName + '.ref'
        simSeqFileName = simFileName + '.fastq'
        simMafFileName = simFileName + '.maf'

        simFilePath = os.path.join(datafolder, simFolder)
        simRefFilePath = os.path.join(simFilePath, simRefFileName)
        simSeqFilePath = os.path.join(simFilePath, simSeqFileName)
        simMafFilePath = os.path.join(simFilePath, simMafFileName)

        if not os.path.exists(simRefFilePath):
            # import pdb
            # pdb.set_trace()
            raise Exception('Reference file for simulated read %s does not exist!' % qname)
        if not os.path.exists(simSeqFilePath):
            # import pdb
            # pdb.set_trace()
            raise Exception('Sequence file for simulated read %s does not exist!' % qname)
        if not os.path.exists(simMafFilePath):
            # import pdb
            # pdb.set_trace()
            raise Exception('Sequence alignment (MAF) for simulated read %s does not exist!' % qname)

        # Reading reference file
        [headers, seqs, quals] = read_fastq(simRefFilePath)
        simGeneName = headers[0]
        annotation = annotation_dict[simGeneName]       # Getting the correct annotation

        if len(samline_list) > len(annotation.items):
            # sys.stderr.write('\nWARNING: A number of partial alignments exceeds the number of exons for query %s! (%d / %d)' % (qname, len(samline_list), len(annotation.items)))
            s_num_oversplit_alignment += 1

        # Reading MAF file to get original position and length of the simulated read
        # Query name should be a second item
        maf_startpos = maf_length = 0
        i = 0
        with open(simMafFilePath, 'rU') as maffile:
            i += 1
            for line in maffile:
                if line[0] == 's':
                    elements = line.split()
                    maf_qname = elements[1]
                    if maf_qname == 'ref':              # Have to remember data for the last reference before the actual read
                        maf_startpos = int(elements[2])
                        maf_length = int(elements[3])
                    if maf_qname == simQName:
                        # maf_startpos = int(elements[2])
                        # maf_length = int(elements[3])
                        break

        if maf_qname != simQName:
            # import pdb
            # pdb.set_trace()
            raise Exception('ERROR: could not find query %s in maf file %s' % (qname, simMafFileName))

        # Calculating expected partial alignmetns from MAF and annotations

        # Saving "maf_length" to be able to check it later
        t_maf_length = maf_length

        # 1. Calculating the index of the first exon
        # i - the index of exon currently being considered
        i = 0
        while annotation.items[i].getLength() < maf_startpos:
            maf_startpos -= annotation.items[i].getLength()
            i += 1

        # Calculating expected partial alignments by filling up exons using maf_length
        expected_partial_alignments = []
        while maf_length > 0:
            # try:
            #     start = annotation.items[i].start + maf_startpos
            #     end = annotation.items[i].end
            # except Exception:
            #     import pdb
            #     pdb.set_trace()
            #if not start < end:
            #    import pdb
            #    pdb.set_trace()
            start = annotation.items[i].start + maf_startpos
            end = annotation.items[i].end
            assert start <= end
            length = end-start+1
            if length <= maf_length:
                expected_partial_alignments.append((start, end))
                maf_length -= length
                i += 1
            else:
                expected_partial_alignments.append((start, start + maf_length))
                maf_length = 0
                i += 1

            # Start position should only be considered for the first exon
            maf_startpos = 0

        numparts = len(expected_partial_alignments)
        # For each part of expected partial alignments, these maps will count
        # how many real partial alignments overlap or equal it
        parthitmap = {(i+1):0 for i in xrange(numparts)}
        parteqmap = {(i+1):0 for i in xrange(numparts)}

        isSplitRead = False
        if len(expected_partial_alignments) > 1:
            s_maf_split_reads += 1
            isSplitRead = True

        if RNAseqEval.getChromName(samline_list[0].rname) != RNAseqEval.getChromName(annotation.seqname):
            # import pdb
            # pdb.set_trace()
            s_num_badchrom_alignments += 1
        else:
            if len(samline_list) != len(expected_partial_alignments):
            # sys.stderr.write('\nWARNING: suspicious number of alignments for query %s!' % qname)
                s_maf_suspicious_alignments += 1
            # import pdb
            # pdb.set_trace()

            good_alignment = True
            k = 0
            for samline in samline_list:
                # sl_startpos = samline.pos - 1   # SAM positions are 1-based
                sl_startpos = samline.pos
                reflength = samline.CalcReferenceLengthFromCigar()
                sl_endpos = sl_startpos + reflength

                # Comparing a samline to the corresponding expected partial alignment
                if k < len(expected_partial_alignments):
                    expected_alignement = expected_partial_alignments[k]
                    maf_startpos = expected_alignement[0]
                    maf_endpos = expected_alignement[1]
                    if abs(sl_startpos - maf_startpos) > allowed_inacc or abs(sl_endpos - maf_endpos) > allowed_inacc:
                        good_alignment = False
                else:
                    good_alignment = False
                k += 1

                # Comparing a samline to all expected partial alignments
                for i in xrange(len(expected_partial_alignments)):
                    expected_alignement = expected_partial_alignments[i]
                    maf_startpos = expected_alignement[0]
                    maf_endpos = expected_alignement[1]

                    if interval_equals((sl_startpos, sl_endpos), (maf_startpos, maf_endpos), allowed_inacc):
                        parteqmap[i+1] += 1
                    if interval_overlaps((sl_startpos, sl_endpos), (maf_startpos, maf_endpos), allowed_inacc):
                        parthitmap[i+1] += 1

            if good_alignment:
                s_maf_good_alignments += 1
                if isSplitRead:
                    s_maf_good_split_alignments += 1
            else:
                # import pdb
                # pdb.set_trace()
                s_maf_bad_alignments += 1
                if isSplitRead:
                    s_maf_bad_split_alignments += 1
                # TODO: check which alignments are bad and why
                # If the choromosome is different its obviously a bad alignment
                if RNAseqEval.getChromName(samline.rname) == RNAseqEval.getChromName(annotation.seqname):
                    # import pdb
                    # pdb.set_trace()
                    pass
                else:
                    s_num_badchrom_alignments += 1


        # Analyzing parthitmap and parteqmap
        oneHit = False
        allHits = True
        oneEq = False
        multiHit = False
        for i in xrange(numparts):
            if parthitmap[i+1] > 0:
                oneHit = True
            if parthitmap[i+1] == 0:
                allHits = False
            if parthitmap[i+1] > 1:
                multiHit = True
            if parteqmap[i+1] > 0:
                oneEq = True

        if oneHit:
            s_maf_hit_one_part += 1
            if isSplitRead:
                s_maf_split_hit_one_part += 1
        if allHits:
            s_maf_hit_all_parts += 1
            if isSplitRead:
                s_maf_split_hit_all_parts += 1
                #import pdb
                #pdb.set_trace()
        if oneEq:
            s_maf_eq_one_part += 1
            if isSplitRead:
                s_maf_split_eq_one_part += 1
        if multiHit:
            s_maf_multihit_parts += 1

        num_start_hits = 0
        num_end_hits = 0
        num_hits = 0

        num_partial_alignements = len(samline_list)
        whole_alignment_hit = False
        for samline in samline_list:
            startpos = samline.pos - 1
            reflength = samline.CalcReferenceLengthFromCigar()
            endpos = startpos + reflength

            if samline.flag & 16 == 0:
                readstrand = Annotation_formats.GFF_STRANDFW
                s_num_fw_strand += 1
            else:
                readstrand = Annotation_formats.GFF_STRANDRV
                s_num_rv_strand += 1

            chromname = RNAseqEval.getChromName(samline.rname)

            if chromname == RNAseqEval.getChromName(annotation.seqname) and readstrand != annotation.strand and annotation.overlapsGene(startpos, endpos):
                s_num_potential_bad_strand += 1

            if chromname == RNAseqEval.getChromName(annotation.seqname) and annotation.overlapsGene(startpos, endpos) and (not P_CHECK_STRAND or readstrand == annotation.strand):
                whole_alignment_hit = True
                s_partial_alignment_hits += 1
            else:
                s_partial_alignment_misses +=1

            # Checking how well partial alignments match exons
            startsItem = False
            endsItem = False
            for item in annotation.items:
                if item.overlapsItem(startpos, endpos):
                    num_hits += 1
                if item.startsItem(startpos, endpos):
                    num_start_hits += 1
                    startsItem = True
                if item.endsItem(startpos, endpos):
                    num_end_hits += 1
                    endsItem = True
                if startsItem and endsItem:
                    s_num_start_end_hits += 1

        s_num_start_hits += num_start_hits
        s_num_end_hits += num_end_hits

        # I'm allowing one start and one end not to match starts and ends of exons
        if (num_hits == num_partial_alignements) and (num_start_hits + num_end_hits >= 2*num_partial_alignements - 2) :
            s_num_good_alignments += 1
        # else:
        #     if num_hits > 0:
        #         import pdb
        #         pdb.set_trace()

        if whole_alignment_hit:
            s_whole_alignment_hits += 1
        else:
            s_whole_alignment_misses += 1


    # Printing out results : NEW
    # Variables names matching RNA benchmark paper
    sys.stdout.write('\n\nAnalysis results:')
    sys.stdout.write('\nOriginal Samlines: %d' % report.num_alignments)
    sys.stdout.write('\nUsable whole alignments (with valid CIGAR string): %d' % len(all_sam_lines))
    sys.stdout.write('\nAnnotations: %d' % len(annotation_dict))
    sys.stdout.write('\nMultiexon genes: %d' % s_num_multiexon_genes)

    sys.stdout.write('\nNumber of exon start hits: %d' % s_num_start_hits)
    sys.stdout.write('\nNumber of exon end hits: %d' % s_num_end_hits)
    sys.stdout.write('\nNumber of exon start and end hits: %d' % s_num_start_end_hits)
    sys.stdout.write('\nNumber of good whole alignments: %d' % s_num_good_alignments)

    sys.stdout.write('\nMAF: Hit both ends: %d' % s_maf_good_alignments)
    sys.stdout.write('\nMAF: Hit all parts: %d' % s_maf_hit_all_parts)
    sys.stdout.write('\nMAF: Hit at least one part: %d' % s_maf_hit_one_part)
    sys.stdout.write('\nMAF: Equals at least one part: %d' % s_maf_eq_one_part)

    sys.stdout.write('\nMAF: Number of split reads: %d' % s_maf_split_reads)
    sys.stdout.write('\nMAF: Hit both ends, SPLIT read: %d' % s_maf_good_split_alignments)
    sys.stdout.write('\nMAF: Hit all parts, SPLIT read: %d' % s_maf_split_hit_all_parts)
    sys.stdout.write('\nMAF: Hit at least one part, SPLIT read: %d' % s_maf_split_hit_one_part)
    sys.stdout.write('\nMAF: Equals at least one part, SPLIT read: %d' % s_maf_split_eq_one_part)

    sys.stdout.write('\nDone!\n')


    # # Printing out results
    # sys.stdout.write('\n\nAnalysis results:')
    # sys.stdout.write('\nOriginal Samlines: %d' % report.num_alignments)
    # sys.stdout.write('\nUsable whole alignments: %d' % len(all_sam_lines))
    # sys.stdout.write('\nSplit alignments: %d' % s_num_split_alignment)
    # sys.stdout.write('\nAnnotations: %d' % len(annotation_dict))
    # sys.stdout.write('\nMultiexon genes: %d' % s_num_multiexon_genes)
    # sys.stdout.write('\nPartial alignment hits: %d' % s_partial_alignment_hits)
    # sys.stdout.write('\nPartial alignment misses: %d' % s_partial_alignment_misses)
    # sys.stdout.write('\nWhole alignment hits: %d' % s_whole_alignment_hits)
    # sys.stdout.write('\nWhole alignment misses: %d' % s_whole_alignment_misses)
    # sys.stdout.write('\nNumber of oversplit alignments: %d' % s_num_oversplit_alignment)
    # sys.stdout.write('\nNumber of exon start hits: %d' % s_num_start_hits)
    # sys.stdout.write('\nNumber of exon end hits: %d' % s_num_end_hits)
    # sys.stdout.write('\nNumber of exon start and end hits: %d' % s_num_start_end_hits)
    # sys.stdout.write('\nNumber of good whole alignments: %d' % s_num_good_alignments)
    # sys.stdout.write('\nNumber of alignments mapped to an incorrect chromosome: %d' % s_num_badchrom_alignments)
    # sys.stdout.write('\nPartial alignments on strand (FW / RV): (%d / %d)' % (s_num_fw_strand, s_num_rv_strand))
    # sys.stdout.write('\nPotential bad strand alignments: %d' % s_num_potential_bad_strand)
    # sys.stdout.write('\nMAF: Suspicious alignments: %d' % s_maf_suspicious_alignments)
    # sys.stdout.write('\nMAF: Hit both ends: %d' % s_maf_good_alignments)
    # sys.stdout.write('\nMAF: Didn\'t hit both ends: %d' % s_maf_bad_alignments)
    # sys.stdout.write('\nMAF: Hit all parts: %d' % s_maf_hit_all_parts)
    # sys.stdout.write('\nMAF: Hit at least one part: %d' % s_maf_hit_one_part)
    # sys.stdout.write('\nMAF: Equals at least one part: %d' % s_maf_eq_one_part)
    # sys.stdout.write('\nMAF: Multihit parts (fragmented) alignments: %d' % s_maf_multihit_parts)
    # sys.stdout.write('\nMAF: Number of split reads: %d' % s_maf_split_reads)
    # sys.stdout.write('\nMAF: Hit both ends, SPLIT alignments: %d' % s_maf_good_split_alignments)
    # sys.stdout.write('\nMAF: Didn\'t hit both ends, SPLIT alignments: %d' % s_maf_bad_split_alignments)
    # sys.stdout.write('\nMAF: Hit all parts on split read: %d' % s_maf_split_hit_all_parts)
    # sys.stdout.write('\nMAF: Hit at least one part on split read: %d' % s_maf_split_hit_one_part)
    # sys.stdout.write('\nMAF: Equals at least one part on split read: %d' % s_maf_split_eq_one_part)
    # sys.stdout.write('\nDone!\n')


def verbose_usage_and_exit():
    sys.stderr.write('Process pbsim data - A tool for processing data generated by pbsim.\n')
    sys.stderr.write('                   - Collects data generated for multiple references.\n')
    sys.stderr.write('                   - And adjusts headers to reflect a reference of origin.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tprocess\n')
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'process'):
        if (len(sys.argv) != 5):
            sys.stderr.write('Processes a folder containing data generated by pbsim.\n')
            sys.stderr.write('Joins all generated reads into a single FASTQ file.\n')
            sys.stderr.write('Expands existing headers with the name of originating reference.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s <pbsim data folder> <results file> <annotations file>\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        datafolder = sys.argv[2]
        resultfile = sys.argv[3]
        annotationfile = sys.argv[4]

        processData(datafolder, resultfile, annotationfile)

    else:
        print 'Invalid mode!'
