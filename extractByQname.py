#! /usr/bin/python

import sys, os

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import utility_sam
from fastqparser import read_fastq

def extractFromSAM(sam_fname, qnames_fname):
    sys.stderr.write('\nLoading qnames file!')
    qnames = []
    qnames_dict = {}
    with open(qnames_fname, 'rU') as qnames_f:
        qnames = qnames_f.readlines()
        qnames_f.close()
    # Creating a dictionary for faster search
    # Also removing '\n' from the end
    for qname in qnames:
        qnames_dict[qname[:-1]] = 1

    sys.stderr.write('\nLoading SAM file!')
    # Loading SAM file into hash
    # Keeping only SAM lines with regular CIGAR string, and sorting them according to position
    qnames_with_multiple_alignments = {}
    [sam_hash, sam_hash_num_lines, sam_hash_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_fname, qnames_with_multiple_alignments)

    sys.stderr.write('\nExtracting ...')
    # Keys in the dictionary (hash) correspond to qnames
    for (samline_key, samline_list) in sam_hash.iteritems():
        if samline_key in qnames_dict:
            for samline in samline_list:
                sys.stdout.write(samline.original_line + '\n')
        else:
            # import pdb
            # pdb.set_trace()
            pass

    sys.stderr.write('\nFinished!')



def extractFromFAST(fast_fname, qnames_fname):

    sys.stderr.write('\nLoading qnames file!')
    qnames = []
    qnames_dict = {}
    with open(qnames_fname, 'rU') as qnames_f:
        qnames = qnames_f.readlines()
        qnames_f.close()
    # Creating a dictionary for faster search
    # Also removing '\n' from the end
    for qname in qnames:
        qnames_dict[qname[:-1]] = 1

    sys.stderr.write('\nLoading FASTA/FASTQ file!')
    [headers, seqs, quals] = read_fastq(fast_fname)

    sys.stderr.write('\nExtracting ...')
    for i in xrange(len(headers)):
        header = headers[i]
        seq = seq[i]
        qual = quals[i]
        qname = header[i]
        if qname in qnames_dict:
            if fext.lower() in ['.fa', '.fna', 'faa', '.fasta']:
                sys.stdout.write('>' + header + '\n')
                sys.stdout.write(seq + '\n')
            elif fext.lower() in ['.fq', '.fastq']:
                sys.stdout.write('@' + header + '\n')
                sys.stdout.write(seq + '\n')
                sys.stdout.write('+' + header + '\n')
                sys.stdout.write(qual + '\n')

    sys.stderr.write('\nFinished!')



def verbose_usage_and_exit():
    sys.stderr.write('extractByQname - extract lines from FASTA/FASTQ or SAM files\n')
    sys.stderr.write('                 relating to given qnames and print them out.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s SAM_or_FASTQ_file qnames_file\n' % sys.argv[0])
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) != 3):
        verbose_usage_and_exit()

    lines_fname = sys.argv[1]
    qnames_fname = sys.argv[2]

    fname, fext = os.path.splitext(lines_fname)
    if fext.lower() in ['.fa', '.fna', 'faa', '.fasta'] or fext.lower() in ['.fq', '.fastq']:
        sys.stderr.write('\nExtracting from FASTA/FASTQ')
        extractFromFAST(lines_fname, qnames_fname)
    elif fext.lower() in ('.sam'):
        sys.stderr.write('\nExtracting from SAM')
        extractFromSAM(lines_fname, qnames_fname)

    # END
