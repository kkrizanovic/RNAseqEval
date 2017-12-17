#! /usr/bin/python

import sys, os
import re

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
from fastqparser import read_fastq
import utility_sam


allowed_inacc = 0
test_flag = 2048

def analyze_chimeric_PAF(filename):

    fname, fext = os.path.splitext(filename)
    count = 0
    chimeric_reads = {}

    if fext != '.PAF' and fext != '.paf':
        raise Exception('File format need to be PAF!: %s' % fext)

    paffile = open(filename, 'rU')
    # sys.stdout.write('\nQNAME\tSTRAND\tQLEN:TLEN\tQSTART:TSTART\tQEND:TEND\n')
    for line in paffile:
        elements = line.split('\t')
        if len(elements) < 12:
            raise Exception('Invalid number of elements in a PAF line: %s' % line)

        qname = elements[0]
        qlen = int(elements[1])
        qstart = int(elements[2])
        qend = int(elements[3])
        strand = elements[4]
        tname = elements[5]
        tlen = int(elements[6])
        tstart = int(elements[7])
        tend = int(elements[8])

        # If read is considered chimeric, print its information
        if (qname == tname):
            assert qlen == tlen
            # if (qstart == tstart and qend == tend):
            if (abs(qstart - tstart) <= allowed_inacc and abs(qend - tend) <= allowed_inacc):
                count += 1
                chimeric_reads[qname] = 1
                # sys.stdout.write('%s\t%s\t%d:%d\t%d:%d\t%d:%d\n' % (qname, strand, qlen, tlen, qstart, tstart, qend, tend))

    paffile.close()

    # KK: printing out only the names of chimeric reads
    for name in chimeric_reads.iterkeys():
        sys.stdout.write('%s\n' % name)

    sys.stderr.write('\n Count occurences: %d' % count)
    sys.stderr.write('\n Count reads: %d' % len(chimeric_reads))


def analyze_chimeric_SAM(filename):

    fname, fext = os.path.splitext(filename)
    count = 0
    chimeric_reads = {}

    if fext != '.SAM' and fext != '.sam':
        raise Exception('File format need to be SAM!: %s' % fext)

    # Loading SAM file into hash
    # Keeping only SAM lines with regular CIGAR string, and sorting them according to position
    qnames_with_multiple_alignments = {}
    [sam_hash, sam_hash_num_lines, sam_hash_num_unique_lines] = utility_sam.HashSAMWithFilter(filename, qnames_with_multiple_alignments)

    for (samline_key, samline_list) in sam_hash.iteritems():
        for samline in samline_list:
            # import pdb
            # pdb.set_trace()
            flag = samline.flag
            # Detecting a chimeric alignment from a SAM file
            if flag & test_flag > 0:
                chimeric_reads[samline.qname] = 1
                count += 1

    # KK: printing out only the names of chimeric reads
    for name in chimeric_reads.iterkeys():
        sys.stdout.write('%s\n' % name)

    sys.stderr.write('\n Count occurences: %d' % count)
    sys.stderr.write('\n Count reads: %d' % len(chimeric_reads))



def split(readsfile, namesfile):
    fname, fext = os.path.splitext(readsfile)

    readsfile1 = fname + '1' + fext
    readsfile2 = fname + '2' + fext

    file1 = open(readsfile1, 'w')
    file2 = open(readsfile2, 'w')

    [headers, seqs, quals] = read_fastq(readsfile)
    names = []
    nfile = open(namesfile, 'rU')
    for line in nfile:
        names.append(line[:-1])

    i = 0
    count1 = count2 = 0
    for i in range(len(headers)):
        header = headers[i]
        # Removing everything after the first space
        pos = header.find(' ')
        header = header[:pos]
        seq = seqs[i]
        qual = quals[i]
        if header in names:
            if fext.lower() in ['.fa', '.fna', 'faa', '.fasta']:
                file1.write('>' + header + '\n')
                file1.write(seq + '\n')
            elif fext.lower() in ['.fq', '.fastq']:
                file1.write('@' + header + '\n')
                file1.write(seq + '\n')
                file1.write('+' + header + '\n')
                file1.write(qual + '\n')
            else:
                raise Exception('Invalid extension for reads file: %s' % fext)
            count1 += 1
        else:
            if fext.lower() in ['.fa', '.fna', 'faa', '.fasta']:
                file2.write('>' + header + '\n')
                file2.write(seq + '\n')
            elif fext.lower() in ['.fq', '.fastq']:
                file2.write('@' + header + '\n')
                file2.write(seq + '\n')
                file2.write('+' + header + '\n')
                file2.write(qual + '\n')
            count2 += 1

        i += 1

    file1.close()
    file2.close()

    sys.stderr.write('\n%d reads in file1; %d reads n file2\n' % (count1, count2))


def verbose_usage_and_exit():
    sys.stderr.write('get_chimer_from_paf - prints out the name of chimeric reads, from a paf file with overlaps.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tanalyze-PAF\n')
    sys.stderr.write('\t\tanalyze-SAM\n')
    sys.stderr.write('\t\tsplit\n')
    sys.stderr.write('\t\tsplit_all\n')
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'analyze-PAF'):
        if (len(sys.argv) != 3):
            sys.stderr.write('Analyzes a PAF file with overlaps to detect chimeric reads.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s <PAF file with overlaps>\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        filename = sys.argv[2]

        analyze_chimeric_PAF(filename)

    if (mode == 'analyze-SAM'):
        if (len(sys.argv) != 3):
            sys.stderr.write('Analyzes a SAM file with mappings to detect chimeric reads.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s <SAM file with mappings>\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        filename = sys.argv[2]

        analyze_chimeric_SAM(filename)

    elif (mode == 'split'):
        if (len(sys.argv) != 4):
            sys.stderr.write('Splits reads in a fasta file into two files according to a set of read names.\n')
            sys.stderr.write('One file contains reads from a given set of names, and the other file contains other reads.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s <FASTA file with reads> <file with read names>\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        readsfile = sys.argv[2]
        namesfile = sys.argv[3]

        split(readsfile, namesfile)

    else:
        print 'Invalid mode!'
