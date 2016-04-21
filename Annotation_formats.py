#! /usr/bin/python

import sys, os

GFF_STRANDFW = '+'
GFF_STRANDRV = '-'
GFF_FRAME = [0, 1, 2]


class GeneItem:
    def __init__(self):
        self.itemName = ''
        self.start = 0
        self.end = 0
        self.frame = 0

    def getLength(self):
        return self.end - self.start

    def insideItem(self, startpos, endpos):
        if startpos >= self.start and endpos <= self.end:
            return True
        else:
            return False

    # this implementation of overlap, will include inside
    def overlapsItem(self, startpos, endpos):
        if endpos <= self.start or startpos >= self.end:
            return False
        else:
            return True

    def basesInside(self, startpos, endpos):
        count = 0
        if self.insideItem(startpos, endpos):
            count = startpos - endpos
        elif self.overlapsItem(startpos, endpos):
            if endpos < self.end:
                count = endpos - self.start
            if startpos > self.start:
                count = self.end - startpos

        return count


class GeneDescription:
    def __init__(self):
        self.seqname = ''
        self.source = ''
        self.genename = ''
        self.transcriptname = ''
        self.strand = GFF_STRANDFW
        self.start = -1
        self.end = -1
        self.score = 0.0
        self.items = []

    def getLength(self):
        return self.end - self.start

    def insideGene(self, startpos, endpos):
        if startpos >= self.start and endpos <= self.end:
            return True
        else:
            return False

    # this implementation of overlap, will include inside
    def overlapsGene(self, startpos, endpos):
        if endpos <= self.start or startpos >= self.end:
            return False
        else:
            return True

    def basesInsideGene(self, startpos, endpos):
        count = 0
        if self.insideGene(startpos, endpos):
            count = startpos - endpos
        elif self.overlapsGene(startpos, endpos):
            if endpos < self.end:
                count = endpos - self.start
            if startpos > self.start:
                count = self.end - startpos

        return count

    def basesInsideItems(self, startpos, endpos):
        for item in self.items:
            if item.insideItem(startpos, endpos):
                return True

        return False

    def overlapsItems(self, startpos, endpos):
        for item in self.items:
            if item.overlapsItem(startpos, endpos):
                return True

    def basesInsideItems(self, startpos, endpos):
        count = 0
        for item in self.items:
            bases += item.basesInside(startpos, endpos)

        return count



class GFFLine:
    def __init__(self):
        self.seqname = ''
        self.source = ''
        self.feature = ''
        self.start = 0
        self.end = 0
        self.score = 0.0
        self.strand = GFF_STRANDFW
        self.frame = 0
        self.attribute = {}

class BEDLine:
    def __init__(self):
        self.chrom = ''
        self.chromStart = -1
        self.chromEnd = -1
        self.name = ''
        self.score = -1
        self.strand = GFF_STRANDFW
        self.thickStart = -1
        self.thickEnd = -1
        self.itemRGB = ''
        self.blockCount = 0
        self.blockSizes = []
        self.blockStarts = []


# Opposed to BED data, annotation data here will have absolute positions
def Annotation_From_BED(bedline):
    genedscp = GeneDescription()
    genedscp.seqname = bedline.chrom
    genedscp.genename = bedline.name
    genedscp.score = bedline.score
    genedscp.strand = bedline.strand
    genedscp.start = bedline.chromStart
    genedscp.end = bedline.chromEnd
    for i in range(bedline.blockCount):
        geneitem = GeneItem()
        geneitem.start = genedscp.start + bedline.blockStarts[i]
        geneitem.end = geneitem.start + bedline.blockSizes[i]
        genedscp.items.append(geneitem)

    return genedscp


# TODO: the implementation is currently Faulty
# It assumes one exon per gene, which is true  for bacteria but not for eucaryota
def Annotation_From_GFF(gffline):
    genedscp = GeneDescription()
    genedscp.seqname = gffline.seqname
    genedscp.source = gffline.source
    genedscp.start = gffline.start
    genedscp.end = gffline.end
    genedscp.score = gffline.score
    genedscp.strand = gffline.strand

    # Extracting from GFF attributes
    genedcsp.genename = gffline.attribute['gene_id']
    genedcsp.transcriptname = gffline.attribute['transcript_id']

    # constructing a single gene item (exon)
    geneitem = GeneItem()
    geneitem.frame = gffline.frame
    geneitem.start = gffline.start
    geneitem.end = gffline.end
    genedscp.items.append(geneitem)

    return genedscp


def Load_Annotation_From_File(filename):

    fname, fext = os.path.splitext(filename)
    if fext == '.gff':
        type = 'GFF'
    elif fext == '.gtf':
        type = 'GTF'
    elif fext == '.bed':
        type = 'BED'
    else:
        raise Exception('Invalid annotation file type: %s' % fext)

    annotations = []

    if type == 'GFF' or type == 'GTF':
        gff_lines = Load_GFF_From_File(filename)
        for gffline in gff_lines:
            annt = Annotation_From_GFF(gffline)
            annotations.append(annt)

    elif type == 'BED':
        bed_lines = Load_BED_From_File(filename)
        for bedline in bed_lines:
            annt = Annotation_From_BED(bedline)
            annotations.append(annt)

    return annotations


def Load_GFF_From_File(filename):
    gff_lines = []
    if not (filename.endswith('.gff') or filename.endswith('.gtf')):
        sys.stderr.write('\nWARNING: file %s does not have GFF/GTF extension!\n' % filename)

    fname, fext = os.path.splitext(filename)
    if fext == '.gff':
        type = 'GFF'
    elif fext == '.gtf':
        type = 'GTF'
    else:
        raise Exception('Invalid annotation file type: %s' % fext)

    file = open(filename, 'rU')
    for line in file:
        elements = line.split('\t')
        gffline = GFFLine()

        if elements[0] == '.':
            gffline.seqname = ''
        else:
            gffline.seqname = elements[0]
        if elements[1] == '.':
            gffline.source = ''
        else:
            gffline.source = elements[1]
        if elements[2] == '.':
            gffline.feature = ''
        else:
            gffline.feature = elements[2]
        if elements[3] == '.':
            gffline.start = 0
        else:
            gffline.start =  int(elements[3])
        if elements[4] == '.':
            gffline.end = 0
        else:
            gffline.end = int(elements[4])
        if elements[5] == '.':
            gffline.score = 0.0
        else:
            gffline.score = float(elements[5])
        if elements[6] not in [GFF_STRANDFW, GFF_STRANDRV]:
            gffline.strand = GFF_STRANDFW
        else:
            gffline.strand = elements[6]
        if elements[7] not in GFF_FRAME:
            gffline.frame = 0
        else:
            gffline.frame = int(elements[7])

        if elements[8] == '.':
            gffline.attribute = {}
        else:
            att_line = elements[8]
            att_list = att_line.split(';')
            for i in xrange(len(att_list)/2):
                gffline.attribute[att_list[2*i]] = att_list[2*i+1]

        # TODO: GFF and GTF contain start and stop codons, CDSs and exons
        # currently using only exons (maybe CDS would be a better choice)
        if gffline.feature == 'exon':
            gff_lines.append(gffline)

    return gff_lines


def Load_BED_From_File(filename):
    bed_lines = []
    if not (filename.endswith('.bed')):
        sys.stderr.write('\nWARNING: file %s does not have BED extension!\n' % filename)

    # Copied from GFF, might be useful in the future
    fname, fext = os.path.splitext(filename)
    if fext == '.bed':
        type = 'BED'
    else:
        raise Exception('Invalid annotation file type: %s' % fext)

    file = open(filename, 'rU')
    for line in file:
        # Ignoring header lines
        if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
            pass
        else:
            elements = line.split()    # splitting with default delimitters
            attcount = len(elements)
            bedline = BEDLine()
            bedline.chrom = elements[0]
            bedline.chromStart = int(elements[1])
            bedline.chromEnd = int(elements[2])
            if attcount >= 4:
                bedline.name = elements[3]
            if attcount >= 5:
                bedline.score = int(elements[4])
            if attcount >= 6:
                bedline.strand = elements[5]
            if attcount >= 7:
                bedline.thickStart = int(elements[6])
            if attcount >= 8:
                bedline.thickEnd = int(elements[7])
            if attcount >= 9:
                bedline.itemRGB = elements[8]
            if attcount >= 10:
                bedline.blockCount = int(elements[9])
            if attcount >= 11:
                if elements[10].endswith(','):
                    elements[10] = elements[10][:-1]
                bedline.blockSizes = [int(el) for el in elements[10].split(',')]
            if attcount >= 12:
                if elements[11].endswith(','):
                    elements[11] = elements[11][:-1]
                bedline.blockStarts = [int(el) for el in elements[11].split(',')]

            bed_lines.append(bedline)

    return bed_lines



if __name__ == "__main__":
	pass;
