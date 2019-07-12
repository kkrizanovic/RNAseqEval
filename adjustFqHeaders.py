#! /usr/bin/python

import sys, os

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import utility_sam
import Annotation_formats

from fastqparser import read_fastq

def adjustFqHeaders(fastqfile, findStr, replaceStr):
    # Reading fastq file
    [headers, seqs, quals] = read_fastq(fastqfile)
    filename, file_extension = os.path.splitext(fastqfile)

    totalSeqs = len(headers)
    findLen = len(findStr)
    replaceLen = len(replaceStr)
    replaced = 0
    notreplaced = 0

    for i in xrange(totalSeqs):
        header = headers[i]
        seq = seqs[i]               # Not really needed
        qual = quals[i]             # Not really needed

        if header[:findLen] == findStr:
            newheader = replaceStr + header[findLen:]
            headers[i] = newheader
            replaced += 1
        else:
            notreplaced += 1

    with open(fastqfile, 'w') as ffile:
        for i in xrange(totalSeqs):
            header = headers[i]
            seq = seqs[i]
            qual = quals[i]

            if file_extension.lower() in ['.fa', '.fna', 'faa', '.fasta']:
                ffile.write('>' + header + '\n')
                ffile.write(seq + '\n')
            elif file_extension.lower() in ['.fq', '.fastq']:
                ffile.write('@' + header + '\n')
                ffile.write(seq + '\n')
                ffile.write('+' + header + '\n')
                ffile.write(qual + '\n')
            else:
                ffile.write(r'@ERROR occured. File is NOT COMPLETE!')
                raise Exception('Invalid file extension: %s' % file_extension)

    return replaced, notreplaced



def verbose_usage_and_exit():
    sys.stderr.write('adjustFqHeaders - Replace a string at the beginning of each fastq header\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [Fastq/Fasta file] [find string] [replace string]\n' % sys.argv[0])
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) != 4):
        verbose_usage_and_exit()

    fastqfile = sys.argv[1]
    findStr = sys.argv[2]
    replaceStr = sys.argv[3]

    replaced, notreplaced = adjustFqHeaders(fastqfile, findStr, replaceStr)

    print('\nStatistics:\n')
    print('Adjusted headers: %d\n' % replaced)
    print('Non adjusted headers: %d\n' % notreplaced)
