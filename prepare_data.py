#! /usr/bin/python

# This file contains mathods for manipulating fasta and annotation files
# To make them appropriate for testing RNAseq mapping

# All headers will be changed to chr[ID]

import sys, os

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import utility_sam

from fastqparser import read_fastq


# String that will disqualify a fasta or annotation line
# All lines containing a bad string will be dropped
bad_strings = []


# Prepare genome reference for drosophila melanogaster
def prepare_dm_genome(genome_file):
    pass



# Prepare genome annotations for drosophila melanogaster
def prepare_dm_annotations(annotations_file):
    pass



# Prepare genome reference for saccharomyces cerevisiae
# Processed file will be added '_P' before extension
def prepare_sc_genome(genome_file):
    filename, file_extension = os.path.splitext(genome_file)
    processed_genome_file = filename + '_P' + file_extension
    [headers, seqs, quals] = read_fastq(genome_file)

    with open(processed_genome_file, 'w') as pgfile:
        for i in range(len(headers)):
            header = headers[i]
            new_header = 'ERROR!'       # In case it somehow slips through
            seq = seqs[i]
            qual = quals[i]
            goodLine = True

            # Check if line contains any disqualifying enteries
            for badstring in bad_strings:
                if header.find(badstring) > -1:
                    goodLine = False
                    break

            # If the line is still good, transform header name to desired form
            if goodLine:
                pos = header.find('chromosome')
                if pos > -1:
                    pos2 = header[pos:].find(',')
                    new_header = 'chr' + header[pos+11:pos+pos2]
                else:
                    pos = header.find('mitochondrion')
                    if pos > -1:
                        new_header = 'chrM'
                    else:
                        # This shouldn't happens
                        raise Exception('Invalid SC genome header: %s!') % header

            if goodLine:
                if file_extension.lower() in ['.fa', '.fna', 'faa', '.fasta']:
                    pgfile.write('>' + new_header + '\n')
                    pgfile.write(seq + '\n')
                elif file_extension.lower() in ['.fq', '.fastq']:
                    pgfile.write('@' + new_header + '\n')
                    pgfile.write(seq + '\n')
                    pgfile.write('+' + new_header + '\n')
                    pgfile.write(qual + '\n')
                else:
                    raise Exception('Invalid file extension: %s' % file_extension)


# Prepare genome annotations for saccharomyces cerevisiae
def prepare_sc_annotations(annotations_file):
    # ATM these annotations seem to be OK
    sys.stderr.write('\nNothing to do with SC annotations!')
    pass


# Prepare genome reference for homo sapiens
def prepare_human_genome(genome_file):
    pass



# Prepare genome annotations for homo sapiens
def prepare_human_annotations(annotations_file):
    pass



def verbose_usage_and_exit():
    sys.stderr.write('RNAseqEval - A tool for evaulating RNAseq results.\n')
    sys.stderr.write('This script is used to prepare genome reference and annotation files.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode] [filename]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tsc-genome - Process S.Cerevisiae genome\n')
    sys.stderr.write('\t\tsc-annotations - Process S.Cerevisiae annotations\n')
    sys.stderr.write('\t\tdm-genome - Process D-Melanogaster genome\n')
    sys.stderr.write('\t\tdm-annotations - Process D-Melanogaster annotations\n')
    sys.stderr.write('\t\th-genome - Process human genome\n')
    sys.stderr.write('\t\th-annotations - Process human annotations\n')
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
