#! /usr/bin/python

import sys, os


# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import utility_sam

from fastqparser import read_fastq


# A function that copies the original file, but replaces the sequences that have a consensus_file
# with those from a consensus file
# Consensus has no quals so they are not important
# Output will be a FASTA file
# NOTE: headers in a consensus file start with "Consensus_"
def fixAfterRacon(consensus_file, original_file, output_file = sys.stdout):
    [cheaders, cseqs, cquals] = read_fastq(consensus_file)
    [oheaders, oseqs, oquals] = read_fastq(original_file)

    clen = len(cheaders)
    olen = len(oheaders)

    for oidx in xrange(olen):
        csame = 0
        oheader = oheaders[oidx]
        oseq = oseqs[oidx]
        for cidx in xrange(clen):
            cheader = cheaders[cidx]
            cseq = cseqs[cidx]
            if oheader == cheader[10:]:
                csame += 1
                # Write consensus sequence to output
                sys.stdout.write('>%s\n' % cheader)
                sys.stdout.write('%s\n' % cseq)

        if csame == 0:
            # Write original sequence to output
            sys.stdout.write('>%s\n' % oheader)
            sys.stdout.write('%s\n' % oseq)

        if csame > 1:
            sys.stderr.write('\nFound an original with %d corresponding consensuses' % csame)
            sys.stderr.write('\n%s' % oheader)


    sys.stdout.write('\n');
    sys.stdout.write('\nNumber of sequences in original file: %d' % olen);
    sys.stdout.write('\nNumber of sequences in consensus file: %d' % clen);
    sys.stdout.write('\n');

    pass


def verbose_usage_and_exit():
    sys.stderr.write('FixAfterRacon - A tool that adds left out sequences to a consensus after aplying Racon.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s consensus_file original_file\n' % sys.argv[0])
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) != 3):
        verbose_usage_and_exit()

    consensus_file = sys.argv[1]
    original_file = sys.argv[2]

    fixAfterRacon(consensus_file, original_file)
