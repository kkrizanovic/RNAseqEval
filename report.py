#! /usr/bin/python


class ReportType:
    FASTA_REPORT = 0
    ANNOTATION_REPORT = 1

class EvalReport:


    def __init__(self, rtype = ReportType.FASTA_REPORT):
        # Basic information
        self.rtype = rtype

        # Reference information
        self.chromlengths = {}        # All chromosome lengths for the reference
        self.reflength = 0
        self.totalGeneLength = 0

        # SAM file statistics
        self.num_alignments = 0
        self.num_unique_alignments = 0
        self.num_multi_alignments = 0
        self.num_split_alignments = 0
        self.num_non_alignments = 0         # Alignments that do not have a CIGAR string

        self.avg_mapping_quality = 0.0      # Not counting quality = 0
        self.min_mapping_quality = 0
        self.max_mapping_quality = 0
        self.num_good_quality = 0           # Number of alignments with quality > 0
        self.num_zero_quality = 0

        # Maybe additional information should be considered
        # i.e. avg indel length
        self.num_match = -1
        self.num_mismatch = -1
        self.num_insert = -1
        self.num_delete = -1
        self.match_percentage = 0.0
        self.mismatch_percentage = 0.0
        self.insert_percentage = 0.0
        self.delete_percentage = 0.0

        # Information gained from annotation file
        self.num_genes = 0
        self.num_exons = 0
        self.max_exons_per_gene = 0         # ATM not used

        self.min_gene_length = 0
        self.max_gene_length = 0
        self.avg_gene_length = 0.0

        self.min_exon_length = 0
        self.max_exon_length = 0
        self.avg_exon_length = 0.0

        # Advanced mapping information
        self.num_genes_covered = 0      # How many genes from the annotation file were
                                        # covered (intersected) by one or more alignments
        self.num_missed_alignments = 0     # Alignements that missed genes, aligned outside genes
        self.num_partial_alignments = 0    # partially overlapping genes
        self.num_hit_alignments = 0        # completely inside genes

        # Information on exons analogous to above information on genes
        self.num_exons_covered = 0
        self.num_exon_hit = 0
        self.num_exon_partial = 0
        self.num_exon_miss = 0
        self.num_multiple_exons = 0         # alignments spanning multiple exons

        self.num_near_miss_alignments = 0  # Alignments that are inside genes, but missed all exons

        self.num_multi_gene_alignments = 0  # Alignments that overlap more then one gene
        self.num_bad_split_alignments = 0   # A number of alignments that are reported as split
                                            # but one or more parts fall outside exons

        # Statistical information on partial alignments
        # TODO: have to see what would be interesting to put in here
        self.avg_alignment_hit_perc = 0.0
        self.avg_exon_hit_perc = 0.0

    def chromosomes(self):
        output = ''
        for chrom, chromlen in self.chromlengths.iteritems():
            output += "%s: %dbp\n" % (chrom, chromlen)

        return output


    def toString(self):
        if self.rtype == ReportType.FASTA_REPORT:
            report = """\n
            Reference format: FASTA
            General information:
            Reference length = %d bp
            Number of chromosomes = %d
            Chromosomes:
            %s
            Number of alignments (total / unique / multi / split) = %d / %d / %d / %d
            Alignments with / without CIGAR string = %d / %d
            Mapping quality without zeroes (avg / min / max) = %.2f / %d / %d
            Alignments with mapping quality (>0 / =0) = %d / %d
            Number of matches / mismatches / inserts / deletes = %d / %d / %d / %d
            Percentage of matches / mismatches / inserts / deletes = %.2f / %.2f / %.2f / %.2f
            """ % (self.reflength, len(self.chromlengths), self.chromosomes(), \
                   self.num_alignments, self.num_unique_alignments, self.num_multi_alignments, self.num_split_alignments, \
                   self.num_alignments - self.num_non_alignments, self.num_non_alignments, \
                   self.avg_mapping_quality, self.min_mapping_quality, self.max_mapping_quality, self.num_good_quality, self.num_zero_quality, \
                   self.num_match, self.num_mismatch, self.num_insert, self.num_delete, \
                   self.match_percentage, self.mismatch_percentage, self.insert_percentage, self.delete_percentage)
            return report + '\n'
        elif self.rtype == ReportType.ANNOTATION_REPORT:
            report = """\n
            Reference format: ANNOTATION
            General information:
            Reference length = %d bp
            Number of chromosomes = %d
            Chromosomes:
            %s
            Number of alignments (total / unique / multi / split) = %d / %d / %d / %d
            Alignments with / without CIGAR string = %d / %d
            Mapping quality without zeroes (avg / min / max) = %.2f / %d / %d
            Alignments with mapping quality (>0 / =0) = %d / %d
            Number of matches / mismatches / inserts / deletes = %d / %d / %d / %d
            Percentage of matches / mismatches / inserts / deletes = %.2f / %.2f / %.2f / %.2f
            """ % (self.reflength, len(self.chromlengths), self.chromosomes(), \
                   self.num_alignments, self.num_unique_alignments, self.num_multi_alignments, self.num_split_alignments, \
                   self.num_alignments - self.num_non_alignments, self.num_non_alignments, \
                   self.avg_mapping_quality, self.min_mapping_quality, self.max_mapping_quality, self.num_good_quality, self.num_zero_quality, \
                   self.num_match, self.num_mismatch, self.num_insert, self.num_delete, \
                   self.match_percentage, self.mismatch_percentage, self.insert_percentage, self.delete_percentage)

            report += """\n
            Annotation statistics:
            Total gene length = %d
            Number of Genes / Exons = %d / %d
            Gene size (Min / Max / Avg) = %d / %d / %.2f
            Exon size (Min / Max / Avg) = %d / %d / %.2f
            """ % (self.totalGeneLength, self.num_genes, self.num_exons, self.min_gene_length, self.max_gene_length, self.avg_gene_length, \
                   self.min_exon_length, self.max_exon_length, self.avg_exon_length)

            report += """\n
            Mapping quality information:
            Genes covered / missed / total = %d / %d / %d
            Exons covered / missed / total = %d / %d / %d
            Alignments on genes hit / partial / missed = %d / %d / %d
            Alignments on exons hit / partial / missed = %d / %d / %d
            Alignments spanning multiple genes = %d
            Alignments spanning multiple exons = %d
            Alignments with hit on gene and miss on exons (near miss) = %d
            Split alignments with exon miss (or partial hit) = %d
            Avg. alignment hit percentage = %.2f
            Avg. exon hit percentage = %.2f
            """ % (self.num_genes_covered, self.num_genes - self.num_genes_covered, self.num_genes, \
                   self.num_exons_covered, self.num_exons - self.num_exons_covered, self.num_exons, \
                   self.num_hit_alignments, self.num_partial_alignments, self.num_missed_alignments, \
                   self.num_exon_hit, self.num_exon_partial, self.num_exon_miss, \
                   self.num_multi_gene_alignments, self.num_multiple_exons, self.num_near_miss_alignments, self.num_bad_split_alignments, \
                   self.avg_exon_hit_perc, self.avg_alignment_hit_perc)

            return report + '\n'
        else:
            return "\nERROR: Report not initialized!\n"


if __name__ == "__main__":
	pass;
