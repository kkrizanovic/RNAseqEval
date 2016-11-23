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


# String that will disqualify a fasta or annotation line
# All lines containing a bad string will be dropped
bad_strings = ['NW_']

bad_strings_annnotations = ['random', 'chrUn']

bad_strings_genomes = ['scaffold', 'patch']

# Predefined values for spliting transcriptomes for used organisms
split_sc = {1: 4000, 2: 1000, 3: 1000}
limits_sc = [4000, 5000, 6000]

split_dm_AS = {1: 1500, 2: 750, 3: 750}
limits_dm_AS = [1500, 2250, 3000]

split_dm_SS = {1: 4500, 2: 1750, 3: 1750}
limits_dm_SS = [4500, 6250, 8000]

split_hChr19_AS = {1: 800, 2: 150, 3: 50}
limits_hChr19_AS = [800, 950, 1000]

split_hChr19_SS = {1: 400, 2: 100, 3: 20}
limits_hChr19_SS = [400, 500, 520]


# A maximum number of alternate splicings to keep per gene
# When preparing transcriptomes for benchmark
# If set to 0, all alternate splcings will be kept
ALTERNATE_SPLICINGS_TO_KEEP = 3

# In some annotations file there are duplicate annotations (different annotations with the same name)
# This parameter decides whether to keep all of them or only the first one
KEEP_DUPLICATES = False



# Prepare genome reference for drosophila melanogaster
def prepare_dm_genome(genome_file):
    filename, file_extension = os.path.ext(genome_file)
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
            for badstring in bad_strings_genomes:
                if header.find(badstring) > -1:
                    goodLine = False
                    break

            if goodLine:
                pos = header.find('chromosome')
                if pos > -1:
                    pos2 = header[pos:].find(' ')
                    pos3 = header[pos+pos2+1:].find(' ')    # Looking for second space
                    if pos3 == -1:
                        new_header = 'chr' + header[pos+pos2+1:]
                    else:
                        new_header = 'chr' + header[pos+pos2+1:pos+pos2+1+pos3]
                elif header.find('chr') > -1:
                    # If we can find chr and not chromosome, assume that this header is as it should be
                    new_header = header
                else:
                    pos = header.find('mitochondrion')
                    if pos > -1:
                        new_header = 'chrM'
                    else:
                        # This shouldn't happens
                        import pdb
                        pdb.set_trace()
                        raise Exception('Invalid DM genome header: %s!') % header

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
                    pgfile.write(r'@ERROR occured. File is NOT COMPLETE!')
                    raise Exception('Invalid file extension: %s' % file_extension)



# Prepare genome annotations for drosophila melanogaster
def prepare_dm_annotations(annotations_file):
    filename, file_extension = os.path.splitext(annotations_file)
    processed_annotations_file = filename + '_P' + file_extension

    with open(processed_annotations_file, 'w') as pafile, open(annotations_file, 'r') as afile:
        for line in afile:
            goodLine = True

            # Check if line contains any disqualifying enteries
            for badstring in bad_strings_annnotations:
                if line.find(badstring) > -1:
                    goodLine = False
                    break

            if goodLine:
                # Lines already contain new line character
                # pafile.write(line + '\n')

                pafile.write(line)



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
                elif header.find('chr') > -1:
                    # If we can find chr and not chromosome, assume that this header is as it should be
                    new_header = header
                else:
                    pos = header.find('mitochondrion')
                    if pos > -1:
                        new_header = 'chrM'
                    else:
                        # This shouldn't happens
                        raise Exception('Invalid SC genome header: %s!' % header)

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
                    pgfile.write(r'@ERROR occured. File is NOT COMPLETE!')
                    raise Exception('Invalid file extension: %s' % file_extension)


# Prepare genome annotations for saccharomyces cerevisiae
def prepare_sc_annotations(annotations_file):
    # ATM these annotations seem to be OK
    sys.stdout.write('\nNothing to do with SC annotations!')
    pass


# Split a transcriptome into 3 parts, to simulate each with different coverage
def split_transcriptome(transcriptome_file):
    # split = {1: 4000, 2: 1000, 3: 1000}     # Split ratio
    # limits = [4000, 5000, 6000]
    split = split_sc
    limits = limits_sc

    filename, file_extension = os.path.splitext(transcriptome_file)
    g1_filename = filename + '_G1' + file_extension
    g2_filename = filename + '_G2' + file_extension
    g3_filename = filename + '_G3' + file_extension
    [headers, seqs, quals] = read_fastq(transcriptome_file)

    total = sum(split.values())
    if len(headers) > total:
        total = len(headers)

    random.seed()

    with open(g1_filename, 'w') as g1file, open(g2_filename, 'w') as g2file, open(g3_filename, 'w') as g3file:
        for i in xrange(len(headers)):
            header = headers[i]
            seq = seqs[i]
            qual = quals[i]
            rnum = random.randint(0, total)     # Generate random number
            gfile = None
            if rnum < limits[0]:
                gfile = g1file
            elif rnum < limits[1]:
                gfile = g2file
            elif rnum < limits[2]:
                gfile = g3file
            else:
                continue       # Skip this sequence

            if file_extension.lower() in ['.fa', '.fna', 'faa', '.fasta']:
                gfile.write('>' + header + '\n')
                gfile.write(seq + '\n')
            elif file_extension.lower() in ['.fq', '.fastq']:
                gfile.write('@' + header + '\n')
                gfile.write(seq + '\n')
                gfile.write('+' + header + '\n')
                gfile.write(qual + '\n')


# Prepare genome reference for homo sapiens
# Using only the chromosome 19 Primary assembly
def prepare_human_genome(genome_file):
    filename, file_extension = os.path.splitext(genome_file)
    processed_genome_file = filename + '_P' + file_extension
    [headers, seqs, quals] = read_fastq(genome_file)

    with open(processed_genome_file, 'w') as pgfile:
        for i in range(len(headers)):
            header = headers[i]
            new_header = 'chr19'
            seq = seqs[i]
            qual = quals[i]

            if header.find('chromosome 19') > -1 and header.find('Primary Assembly') > -1:

                if file_extension.lower() in ['.fa', '.fna', 'faa', '.fasta']:
                    pgfile.write('>' + new_header + '\n')
                    pgfile.write(seq + '\n')
                elif file_extension.lower() in ['.fq', '.fastq']:
                    pgfile.write('@' + new_header + '\n')
                    pgfile.write(seq + '\n')
                    pgfile.write('+' + new_header + '\n')
                    pgfile.write(qual + '\n')
                else:
                    pgfile.write(r'@ERROR occured. File is NOT COMPLETE!')
                    raise Exception('Invalid file extension: %s' % file_extension)

                break


# Prepare genome annotations for homo sapiens
# Using only annotations for chromosome 19 Primary assembly
def prepare_human_annotations(annotations_file):
    filename, file_extension = os.path.splitext(annotations_file)
    processed_annotations_file = filename + '_P' + file_extension

    with open(processed_annotations_file, 'w') as pafile, open(annotations_file, 'r') as afile:
        for line in afile:
            goodLine = True

            if goodLine and line.find('chr19') > -1:
                # Lines already contain new line character
                # pafile.write(line + '\n')

                pafile.write(line)



def split_alternate(annotations_file):

    filename, file_extension = os.path.splitext(annotations_file)
    processed_annotations_file_AS = filename + '_AS' + file_extension
    processed_annotations_file_SS = filename + '_SS' + file_extension

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
    # Groupign annotations which overlap and are on the same strand
    start_new_group = True
    grouped_annotations = []

    for new_annotation in annotations:
        if start_new_group:
            annotation_group = []
            annotation_group.append(new_annotation)
            group_start = new_annotation.start
            group_end = new_annotation.end
            group_strand = new_annotation.strand
            group_chrom = new_annotation.seqname
            start_new_group = False
        else:
            if new_annotation.overlapsGene(group_start, group_end) and group_strand == new_annotation.strand and group_chrom == new_annotation.seqname:
                # Add annotation to current group
                annotation_group.append(new_annotation)
                # Adjust group start and end
                if new_annotation.start < group_start:
                    group_start = new_annotation.start
                if new_annotation.end > group_end:
                    group_end = new_annotation.end
            else:
                # Save the current group and start the new one
                grouped_annotations.append(annotation_group)
                annotation_group = []
                annotation_group.append(new_annotation)
                group_start = new_annotation.start
                group_end = new_annotation.end
                group_strand = new_annotation.strand
                group_chrom = new_annotation.seqname

    # At the end, add last group if it exists
    if len(annotation_group) > 0:
        grouped_annotations.append(annotation_group)


    # Annotations with alternate splicing
    as_annotations = []

    # Annotation with single splicing
    ss_annotations = []

    # Separate annotations into those for genes with alternate splicing and genes with single splicing
    # For genes with alternate splicing, keep only ALTERNATE_SPLICINGS_TO_KEEP annotations
    # Have to watch out for duplicate annotation names. Annotation is considered only the first times
    # It enters a list. If it already exists in any list, it is ignored.
    duplicate_genename = False
    for annotation_group in grouped_annotations:
        if len(annotation_group) > 1:
            i = 0
            tr = 0
            for annotation in annotation_group:
                if ALTERNATE_SPLICINGS_TO_KEEP > 0 and ALTERNATE_SPLICINGS_TO_KEEP <= tr:
                    break
                if annotation_group[i].genename in ss_annotations or annotation_group[i].genename in as_annotations:
                    duplicate_genename = True
                else:
                    as_annotations.append(annotation_group[i].genename)
                    tr += 1
                i += 1
        else:
            if annotation_group[0].genename in ss_annotations or annotation_group[0].genename in as_annotations:
                duplicate_genename = True
            else:
                ss_annotations.append(annotation_group[0].genename)

    if duplicate_genename:
        sys.stderr.write('\nWARNING: there were duplicate annotations!\n')


    # Reading original annotations file and writing lines in separate files for
    # single-splicing and alternate-splicing
    # Variable old_genename is used to detect genename change in gtf files, in case duplicate annotations need to be skipped.
    # It is assumed that duplicate genename enteries do not come one after the other (there are other enteries inbetween)
    old_genename = ''
    with open(processed_annotations_file_AS, 'w') as pafile_AS, open(processed_annotations_file_SS, 'w') as pafile_SS, open(annotations_file) as afile:
        for line in afile:
            is_AS = False
            is_SS = False

            count = 0       # Used for sanity check
            # extracting genename from annotation line
            if filetype == 'BED':
                if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                    pass
                else:
                    elements = line.split()
                    genename = elements[3]
            elif filetype == 'GTF':
                genename = 'Unknown'
                elements = line.split('\t')
                att_line = elements[8]
                att_list = att_line.split(';')          # Separating attribute definitions
                for i in xrange(len(att_list)):
                    elements = att_list[i].split()      # Separating key and value for each attribute
                    if len(elements) > 1 and elements[0] == 'transcript_id':
                        genename = elements[1][1:-1]


            # Checking if the line is for an alternate spliced gene
            if genename in as_annotations:
                is_AS = True
                pafile_AS.write(line)
                if not KEEP_DUPLICATES:
                    if filetype == 'BED':
                        as_annotations.remove(genename)

            # Checking if the line is for a single spliced gene
            if genename in ss_annotations:
                is_SS = True
                pafile_SS.write(line)
                if not KEEP_DUPLICATES:
                    if filetype == 'BED':
                        ss_annotations.remove(genename)

            if not KEEP_DUPLICATES and filetype == 'GTF' and old_genename != '' and old_genename != genename:
                if old_genename in as_annotations:
                    as_annotations.remove(old_genename)
                if old_genename in ss_annotations:
                    ss_annotations.remove(old_genename)

            old_genename = genename

            # For testing purposes
            # if (not is_AS) and (not is_SS):
            #     import pdb
            #     pdb.set_trace()

            if is_AS and is_SS:
                sys.stderr.write('\nERROR: genename found in both lists (single splices and alternate spliced)\n')
                sys.stderr.write(line)


# Expand each header in a given fasta/fastq file with a given string.
# String is added at the beginning of a header.
def expandHeader(fastfile, sstring):
    filename, file_extension = os.path.splitext(fastfile)
    [headers, seqs, quals] = read_fastq(fastfile)

    with open(fastfile, 'w') as ffile:
        for i in range(len(headers)):
            header = headers[i]
            new_header = sstring + header
            seq = seqs[i]
            qual = quals[i]

            if file_extension.lower() in ['.fa', '.fna', 'faa', '.fasta']:
                ffile.write('>' + new_header + '\n')
                pgfffileile.write(seq + '\n')
            elif file_extension.lower() in ['.fq', '.fastq']:
                ffile.write('@' + new_header + '\n')
                ffile.write(seq + '\n')
                ffile.write('+' + new_header + '\n')
                ffile.write(qual + '\n')
            else:
                ffile.write(r'@ERROR occured. File is NOT COMPLETE!')
                raise Exception('Invalid file extension: %s' % file_extension)



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
    sys.stderr.write('\t\ttrans-split - Split a transcriptome into 3 groups\n')
    sys.stderr.write('\t\tdm-genome - Process D-Melanogaster genome\n')
    sys.stderr.write('\t\tdm-annotations - Process D-Melanogaster annotations\n')
    sys.stderr.write('\t\th-genome - Process human genome\n')
    sys.stderr.write('\t\th-annotations - Process human annotations\n')
    sys.stderr.write('\t\tsplit-alternate - Split transcriptome into transcripts with only one splicing\n')
    sys.stderr.write('\t\t                  and transcripts with alternate splicings.\n')
    sys.stderr.write('\t\t                  Keep only a given number of alternate splicings.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Alternate usage:\n')
    sys.stderr.write('\t%s expand-headers [filename] [string]\n' % sys.argv[0])
    sys.stderr.write('\t\t - Add a string at the beginning of each header in a fasta/fastq file\n')
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 3):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'sc-genome'):
        genome_file = sys.argv[2]
        prepare_sc_genome(genome_file)

    elif (mode == 'sc-annotations'):
        annotations_file = sys.argv[2]
        prepare_sc_annotations(annotations_file)

    elif (mode == 'trans-split'):
        transcriptome_file = sys.argv[2]
        split_transcriptome(transcriptome_file)

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

    elif (mode == 'split-alternate'):
        annotations_file = sys.argv[2]
        split_alternate(annotations_file)

    elif (mode == 'expand-headers'):
        fastfile = sys.argv[2]
        sstring = sys.argv[3]
        expandHeader(fastfile, sstring)

    else:
        print 'Invalid mode: %s!' % mode
