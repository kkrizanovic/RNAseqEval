#! /usr/bin/python

import sys, os
import re

# For copying SAM lines
import copy
from datetime import datetime

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import utility_sam
import Annotation_formats

from RNAseqEval import load_and_process_reference, getChromName
from fastqparser import read_fastq
from report import EvalReport, ReportType

def test_cigars(samfile, fastaref):

    paramdict = {}
    report = EvalReport(ReportType.TEMP_REPORT)

    sys.stderr.write('\n(%s) Loading and processing FASTA reference ... ' % datetime.now().time().isoformat())
    [chromname2seq, headers, seqs, quals] = load_and_process_reference(fastaref, paramdict, report)

    sys.stderr.write('\n(%s) Loading and processing SAM file with mappings ... ' % datetime.now().time().isoformat())
    qnames_with_multiple_alignments = {}
    [sam_hash, sam_hash_num_lines, sam_hash_num_unique_lines] = utility_sam.HashSAMWithFilter(samfile, qnames_with_multiple_alignments)

    sys.stdout.write('\nTYPE\tQNAME\tMAX numMatch\tLENGTH\tFLAG\n')
    for (samline_key, samline_list) in sam_hash.iteritems():
        if samline_list[0].cigar <> '*' and samline_list[0].cigar <> '':
            for samline in samline_list:
                chromname = getChromName(samline.rname)
                if chromname not in chromname2seq:
                    # import pdb
                    # pdb.set_trace()
                    raise Exception('\nERROR: Unknown chromosome name in SAM file! (chromname:"%s", samline.rname:"%s")' % (chromname, samline.rname))
                chromidx = chromname2seq[chromname]
                cigar = samline.cigar
                length = samline.CalcReadLengthFromCigar()
                numMatch = numMatch1 = numMatch2 = 0
                flag = -1
                try:
                    # Using regular expressions to find repeating digit and skipping one character after that
                    # Used to separate CIGAR string into individual operations
                    pattern = '(\d+)(.)'
                    pos = samline.pos
                    flag = samline.flag

                    # Calculating regular matches
                    extcigar = samline.CalcExtendedCIGAR(seqs[chromidx])
                    operations = re.findall(pattern, extcigar)
                    for op in operations:
                        if op[1] in ('M', '='):
                            numMatch += int(op[0])
                        elif op[1] in ('I', 'D', 'X', 'N', 'S', 'H', 'P'):
                            pass
                        else:
                            sys.stderr.write('\nERROR: Invalid CIGAR string operation (%s)' % op[1])

                    # Calculating for pos + 1 
                    samline.pos = pos + 1
                    extcigar = samline.CalcExtendedCIGAR(seqs[chromidx])
                    operations = re.findall(pattern, extcigar)
                    for op in operations:
                        if op[1] in ('M', '='):
                            numMatch1 += int(op[0])
                        elif op[1] in ('I', 'D', 'X', 'N', 'S', 'H', 'P'):
                            pass
                        else:
                            sys.stderr.write('\nERROR: Invalid CIGAR string operation (%s)' % op[1])

                    # Calculating for pos - 1
                    samline.pos = pos - 1
                    extcigar = samline.CalcExtendedCIGAR(seqs[chromidx])
                    operations = re.findall(pattern, extcigar)
                    for op in operations:
                        if op[1] in ('M', '='):
                            numMatch2 += int(op[0])
                        elif op[1] in ('I', 'D', 'X', 'N', 'S', 'H', 'P'):
                            pass
                        else:
                            sys.stderr.write('\nERROR: Invalid CIGAR string operation (%s)' % op[1])

                except Exception, Argument:
                    # import pdb
                    # pdb.set_trace()
                    sys.stderr.write('ERROR: querry/ref/pos/message = %s/%s/%d/%s \n' % (samline.qname, samline.rname, samline.pos, Argument))
                    pass

                if (numMatch > numMatch1 and numMatch > numMatch2):
                    sys.stdout.write('REGULAR\t%s\t%d\t%d\t%d\n' % (samline.qname, numMatch, length, flag))
                elif (numMatch1 > numMatch and numMatch1 > numMatch2):
                    sys.stdout.write('PLUS ONE\t%s\t%d\t%d\t%d\n' % (samline.qname, numMatch1, length, flag))
                elif (numMatch2 > numMatch and numMatch2 > numMatch1):
                    sys.stdout.write('MINUS ONE\t%s\t%d\t%d\t%d\n' % (samline.qname, numMatch2, length, flag))
                else:
                    sys.stdout.write('NONE\t%s\t%d\t%d\t%d\n' % (samline.qname, numMatch, length, flag))



def verbose_usage_and_exit():
    sys.stderr.write('test_cigar - a script for testing cigar strings.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    # sys.stderr.write('\t\tsetup\n')
    # sys.stderr.write('\t\tcleanup\n')
    sys.stderr.write('\t\ttest-pos\n')
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'test-pos'):
        if (len(sys.argv) != 4):
            sys.stderr.write('Test CIGAR strings. Check if changing position results in a "better" CIGAR\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s test-pos [SAM file] [FASTA reference]\n' % sys.argv[0])
            sys.stderr.write('\n')
            exit(1)

        samfile = sys.argv[2]
        fastaref = sys.argv[3]
        test_cigars(samfile, fastaref)

