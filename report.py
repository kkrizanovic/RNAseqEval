#! /usr/bin/python


class ReportType:
    FASTA_REPORT = 0
    ANNOTATION_REPORT = 1

class EvalReport:


    def __init__(self, rtype = ReportType.FASTA_REPORT):
        # Basic information
        self.rtype = rtype
        self.reflength = 0

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
        self.num_genes_covered = -1     # How many genes from the annotation file were
                                        # covered (intersected) by one or more alignments

        self.num_missed_alignments = -1     # Alignements that missed genes, aligned outside genes
        self.num_partial_alignments = -1    # partially overlapping genes
        self.num_hit_alignments = -1        # completely inside genes

        self.avg_alignment_hit_perc = 0.0

        self.num_exons_covered = -1
        self.num_near_miss_alignments = -1  # Alignments that are inside genes, but missed all exons

    def toString(self):
        if self.rtype == ReportType.FASTA_REPORT:
            return """\n
            Reference format: FASTA
            General information:
            Reference length = %d bp
            Number of alignments (total / unique / multi / split / non) = %d / %d / %d / %d / %d
            Mapping quality (avg / min / max) = %.2f / %d / %d
            Alignments with mapping quality (>0 / =0) = %d / %d
            """ % (self.reflength, self.num_alignments, self.num_unique_alignments, self.num_multi_alignments, self.num_split_alignments, self.num_non_alignments, \
                   self.avg_mapping_quality, self.min_mapping_quality, self.max_mapping_quality, self.num_good_quality, self.num_zero_quality)
        elif self.rtype == ReportType.ANNOTATION_REPORT:
            report = """\n
            Reference format: ANNOTATION
            General information:
            Reference length = %d bp
            Number of alignments (total / unique / multi / split / non) = %d / %d / %d / %d / %d
            Mapping quality (avg / min / max) = %.2f / %d / %d
            Alignments with mapping quality (>0 / =0) = %d / %d
            """ % (self.reflength, self.num_alignments, self.num_unique_alignments, self.num_multi_alignments, self.num_split_alignments, self.num_non_alignments, \
                   self.avg_mapping_quality, self.min_mapping_quality, self.max_mapping_quality, self.num_good_quality, self.num_zero_quality)

            report += """\n
            Annotation statistics:
            Number of Genes / Exons = %d / %d
            Gene size (Min / Max / Avg) = %d / %d / %.2f
            Exon size (Min / Max / Avg) = %d / %d / %.2f
            """ % (self.num_genes, self.num_exons, self.min_gene_length, self.max_gene_length, self.avg_gene_length, \
                   self.min_exon_length, self.max_exon_length, self.avg_exon_length)

            report += """\n
            Mapping quality information:
            Genes covered / missed / total = %d / %d / %d
            Exons covered / missed / total = %d / %d / %d
            Alignments hit / partial / missed / near missed = %d / %d / %d / %d
            Avg. alignment hit percentage = %.2f
            """ % (self.num_genes_covered, self.num_genes - self.num_genes_covered, self.num_genes, \
                   self.num_exons_covered, self.num_exons - self.num_exons_covered, self.num_exons, \
                   self.num_hit_alignments, self.num_partial_alignments, self.num_missed_alignments, self.num_near_miss_alignments, \
                   self.avg_alignment_hit_perc)

            return report
        else:
            return "\nERROR: Report not initialized!\n"


if __name__ == "__main__":
	pass;
