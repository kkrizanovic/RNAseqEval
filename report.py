#! /usr/bin/python


class ReportType:
    FASTA_REPORT = 0
    MAPPING_REPORT = 1
    ANNOTATION_REPORT = 2
    TEMP_REPORT = 10        # Report used to temporarily store some data

class EvalReport:


    def __init__(self, rtype = ReportType.FASTA_REPORT):
        self.commandline = ''       # Command which generated the report

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
        self.num_possibly_split_alignements = 0     # Multiple alignments that are candidate split alignments
                                                    # Some aligners, such as BWA, do not generate split alignments
                                                    # using multiple alignments with clipping instead
        self.num_split_alignments = 0
        self.num_non_alignments = 0         # Alignments that do not have a CIGAR string
        self.num_real_alignments = 0        # A final number of actual alignments with a correct CIGAR string
                                            # To get this, alignments without a CIGAR string are dropped, several alignments of the
                                            # same query are grouped into a split alignement and only the best alignment for each
                                            # query is taken
        self.num_real_split_alignments = 0

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
        self.num_multiexon_genes = 0
        self.max_exons_per_gene = 0         # ATM not used

        self.min_gene_length = 0
        self.max_gene_length = 0
        self.avg_gene_length = 0.0

        self.min_exon_length = 0
        self.max_exon_length = 0
        self.avg_exon_length = 0.0

        # Information on groupped annotations / alternate splicings
        self.num_annotation_groups = 0
        self.num_alternate_spliced_genes = 0
        self.max_spliced_alignments = 0
        self.min_spliced_alignments = 0
        self.max_spliced_exons = 0
        self.min_spliced_exons = 0

        # Advanced mapping information
        self.allowed_inacc = 0
        self.sum_read_length = 0
        self.sum_bases_aligned = 0
        self.percentage_bases_aligned = 0.0

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
        self.num_halfbases_hit = 0          # The number of reads whose alignment falling within an annotation
                                            # covers at least half of the bases of the read

        self.num_lowmatchcnt = 0            # The number of reads with a very low match count from extended CIGAR
                                            # the number of matches is lower then the number of matches, inserts and deletes 
        self.num_partial_exon_miss = 0      # Then number of alignments that have a disjoined part completely missing an exon
                                            # This disjoined part (partial alignment) needs to be greater than allowed_inaccuracy 
                                            # To be considered
        self.num_almost_good = 0            # The number of alignments that would be contiguous, but have at least one partial alignment
                                            # that missses all exons


        self.num_multi_exon_alignments = 0      # alignments spanning multiple exons / this can be correct
        self.num_multi_gene_alignments = 0      # Alignments that overlap more then one gene / this is definitely incorrect

        self.num_inside_miss_alignments = 0     # Alignments that are inside genes, but missed all exons

        self.num_bad_split_alignments = 0   # A number of alignments that are reported as split, but do not cover exons
                                            # as expected
        self.num_good_split_alignments = 0


        # Statistical information on partial alignments
        # TODO: have to see what would be interesting to put in here
        self.avg_alignment_hit_perc = 0.0
        self.avg_exon_hit_perc = 0.0

        # Statistical information on split alignments
        # NOTE: This variable already exists
        # self.num_split_alignments
        self.num_cover_all_exons = 0
        self.num_cover_some_exons = 0
        self.num_cover_no_exons = 0
        self.num_multicover_exons = 0       # Exons covered by more than one part of a split alignment
        self.num_equal_exons = 0            # Exons compeltely covered by a single alignment
        self.num_partial_exons = 0          # Alignments partially covering an exon
        self.num_undercover_alignments = 0
        self.num_overcover_alignments = 0
        self.num_possible_spliced_alignment = 0     # A number of alignments that are possibly spliced
                                                    # I.e. skipping an one or more exons between two covered exons
        self.num_hit_all = 0				# A number of alignemnts that hit a series of exons without skipping any
        									# Each exon does not need to be covered correctly, only "hit"
        									# Should be fairly similar to Hit_all statistics for simulated reads
        self.hit_all_percent = 0.0

        # TODO: Some statistics about a precision of start and end exon points
        self.num_good_starts = 0
        self.num_good_ends = 0

        # Final information that sums up the statistices
        self.num_good_alignment = 0
        self.num_bad_alignment = 0
        self.good_alignment_percent = 0.0
        self.bad_alignment_percent = 0.0

        # Gene expression
        # Looking at expressed genes, ones that overlap with at least one read in SAM file
        # Storing them in a dictionary together with a number of hits for each exon in the gene
        self.expressed_genes = {}
        self.output_gene_expression = False

        # Calculating gene coverage, how many bases of each gene and exon are covered by ready
        # Bases covered multiple times are taken into account multiple times
        # The structure of this dictionary is similar to expressed_genes above
        # Each gene has one global counter (index 0), and one counter for each exon
        self.gene_coverage = {}

        # Alternate splicing
        # Information on genes that have alternate splicing
        self.alternate_splicing = {}
        self.output_alternate_splicing = True

        # Lists for saving qnames
        # - qnames that overlap one exon
        # - qnames that have contiguous alignment
        # - qames that have incorrect alignment
        self.hitone_names = []
        self.hithalfbases_names = []
        self.contig_names = []
        self.incorr_names = []
        self.unmapped_names = []

        # Number of alignments after preprocessing - the number of alignments which are evaluated
        self.num_evaluated_alignments = 0


    def chromosomes(self):
        output = '\t\t'
        for chrom in sorted(self.chromlengths.keys()):
            chromlen = self.chromlengths[chrom]

            output += "\t\t%s: %dbp\n" % (chrom, chromlen)

        return output

    def get_hitone_names(self):
        report = ''
        for name in self.hitone_names:
            report += name + '\n'
        return report

    def get_hithalfbases_names(self):
        report = ''
        for name in self.hithalfbases_names:
            report += name + '\n'
        return report    

    def get_contig_names(self):
        report = ''
        for name in self.contig_names:
            report += name + '\n'
        return report

    def get_incorr_names(self):
        report = ''
        for name in self.incorr_names:
            report += name + '\n'
        return report

    def get_unmapped_names(self):
        report = ''
        for name in self.unmapped_names:
            report += name + '\n'
        return report

    # New toString function for printing out reports
    # Values reported are consistent with the benchmark paper!
    def toString(self):
        if self.rtype == ReportType.FASTA_REPORT:
            report = """\n
            Reference format: FASTA
            General information:
            Reference length = %d bp
            Number of chromosomes = %d
            Chromosomes:\n%s
            Number of alignments (total / unique) = %d / %d
            Alignments with / without CIGAR string = %d / %d
            Mapping quality without zeroes (avg / min / max) = %.2f / %d / %d
            Alignments with mapping quality (>0 / =0) = %d / %d
            Number of matches / mismatches / inserts / deletes = %d / %d / %d / %d
            Percentage of matches / mismatches / inserts / deletes = %.2f / %.2f / %.2f / %.2f
            """ % (self.reflength, len(self.chromlengths), self.chromosomes(), \
                   self.num_alignments, self.num_unique_alignments, \
                   self.num_alignments - self.num_non_alignments, self.num_non_alignments, \
                   self.avg_mapping_quality, self.min_mapping_quality, self.max_mapping_quality, self.num_good_quality, self.num_zero_quality, \
                   self.num_match, self.num_mismatch, self.num_insert, self.num_delete, \
                   self.match_percentage, self.mismatch_percentage, self.insert_percentage, self.delete_percentage)
            return report + '\n'
        elif self.rtype == ReportType.MAPPING_REPORT:
            report = """\n
            Command Line: %s
            Reference format: ANNOTATION
            General information:
            Reference length = %d bp
            Number of chromosomes = %d
            Chromosomes:
            %s
            Number of alignments in SAM file (total / unique) = %d / %d
            Alignments with / without CIGAR string = %d / %d
            Mapping quality without zeroes (avg / min / max) = %.2f / %d / %d
            Alignments with mapping quality (>0 / =0) = %d / %d
            Number of matches / mismatches / inserts / deletes = %d / %d / %d / %d
            Percentage of matches / mismatches / inserts / deletes = %.2f / %.2f / %.2f / %.2f
            """ % (self.commandline, self.reflength, len(self.chromlengths), self.chromosomes(), \
                   self.num_alignments, self.num_unique_alignments, \
                   self.num_alignments - self.num_non_alignments, self.num_non_alignments, \
                   self.avg_mapping_quality, self.min_mapping_quality, self.max_mapping_quality, self.num_good_quality, self.num_zero_quality, \
                   self.num_match, self.num_mismatch, self.num_insert, self.num_delete, \
                   self.match_percentage, self.mismatch_percentage, self.insert_percentage, self.delete_percentage)

            report += """\n
            Annotation statistics:
            Total gene length = %d
            Number of Transcripts / Exons (Multiexon transcripts) = %d / %d (%d)
            Maximum number of exons in a gene = %d
            Gene size (Min / Max / Avg) = %d / %d / %.2f
            Exon size (Min / Max / Avg) = %d / %d / %.2f
            """ % (self.totalGeneLength, self.num_genes, self.num_exons, self.num_multiexon_genes, self.max_exons_per_gene, \
                   self.min_gene_length, self.max_gene_length, self.avg_gene_length, \
                   self.min_exon_length, self.max_exon_length, self.avg_exon_length)

            # report += """\n
            # Grouped annotation (alternate splicing) statistics:
            # Number of annotation groups (genes) = %d
            # Number of genes with alternate splicing = %d
            # Maximum / minimum number of alternate spliced alignments for a gene = %d / %d
            # Maximum / minimum number of exons in spliced alignments = %d / %d
            # """ % (self.num_annotation_groups, self.num_alternate_spliced_genes, \
            #        self.max_spliced_alignments, self.min_spliced_alignments, \
            #        self.max_spliced_exons, self.min_spliced_exons)

            report += """\n
            Mapping quality information:
            Bases in reads (aligned / total) (percent) = (%d / %d) (%.2f%%)
            Transcripts covered / missed / total = %d / %d / %d
            Exons covered / missed / total = %d / %d / %d
            Total number of evaluated alignments = %d
            Alignments on transcript hit / missed = %d / %d
            Alignments on exons hit / missed = %d / %d
            Alignments for reads where more than 50%% bases falls within an annotation = %d
            Alignments with low match count (lower than mismatch, insert and delete combined) = %d
            Alignments hitting an exon (start / end / both) = %d / %d / %d
            Alignments with a partial miss (not good / "almost "good) = %d / %d
            Contiguous / non contiguous alignments: %d (%.2f%%) / %d (%.2f%%)
            Hit all for real reads: %d (%.2f%%)
            """ % (self.sum_bases_aligned, self.sum_read_length, self.percentage_bases_aligned, \
                   self.num_genes_covered, self.num_genes - self.num_genes_covered, self.num_genes, \
                   self.num_exons_covered, self.num_exons - self.num_exons_covered, self.num_exons, \
                   self.num_evaluated_alignments, \
                   self.num_hit_alignments, self.num_missed_alignments, \
                   self.num_exon_hit, self.num_exon_miss, \
                   self.num_halfbases_hit, \
                   self.num_lowmatchcnt, \
                   self.num_good_starts, self.num_good_ends, self.num_equal_exons, \
                   self.num_partial_exon_miss, self.num_almost_good, \
                   self.num_good_alignment, self.good_alignment_percent, self.num_bad_alignment, self.bad_alignment_percent, \
            	   self.num_hit_all, self.hit_all_percent)

            # Counting the number of expressed genes
            # This has already been written in the report, but this was the value can be double checked
            exp_gn_cnt = 0
            for expression in self.expressed_genes.values():
                if expression[0] > 0:
                    exp_gn_cnt += 1

            if self.output_gene_expression:
                report += """\n
            Transcript/exon expression and coverage information:
            Number of expressed Transcripts = %d
            Expressed transcripts:
            genename  number_of_exons  gene_hits / gene_covered_bases  [(exon_hits / exon_covered_bases)]...
                """ % exp_gn_cnt

                if len(self.expressed_genes) != len(self.gene_coverage):
                    raise Exception('ERROR: Gene expression and gene coverage dictionaries do not match in length! (%d <> %d)' % (len(expressed_genes), len(gene_coverage)))

                for genename in self.expressed_genes.keys():
                    numexons = len(self.expressed_genes[genename]) - 1
                    genehits = self.expressed_genes[genename][0]
                    gene_cov_bs = self.gene_coverage[genename][0]
                    if genehits > 0:
                        reportline = '%s  %d  %d  %d' % (genename, numexons, genehits, gene_cov_bs)
                        for i in range(1, numexons+1):
                            exonhits = self.expressed_genes[genename][i]
                            exon_cov_bs = self.gene_coverage[genename][i]
                            reportline += '  (%d / %d)' % (exonhits, exon_cov_bs)

                        report += reportline + '\n'

            return report + '\n'
        elif self.rtype == ReportType.ANNOTATION_REPORT:
            report = """\n
            Command Line: %s
            Reference format: ANNOTATION
            General information:
            Reference length = %d bp
            Number of chromosomes = %d
            Chromosomes:
            %s
            """ % (self.commandline, self.reflength, len(self.chromlengths), self.chromosomes())

            report += """\n
            Annotation statistics:
            Total gene length = %d
            Number of Genes / Exons (Multiexon genes) = %d / %d (%d)
            Maximum number of exons in a gene = %d
            Gene size (Min / Max / Avg) = %d / %d / %.2f
            Exon size (Min / Max / Avg) = %d / %d / %.2f
            """ % (self.totalGeneLength, self.num_genes, self.num_exons, self.num_multiexon_genes, self.max_exons_per_gene, \
                   self.min_gene_length, self.max_gene_length, self.avg_gene_length, \
                   self.min_exon_length, self.max_exon_length, self.avg_exon_length)

            # report += """\n
            # Grouped annotation (alternate splicing) statistics:
            # Number of annotation groups (genes) = %d
            # Number of genes with alternate splicing = %d
            # Maximum / minimum number of alternate spliced alignments for a gene = %d / %d
            # """ % (self.num_annotation_groups, self.num_alternate_spliced_genes, \
            #        self.max_spliced_alignments, self.min_spliced_alignments)

            if self.output_alternate_splicing:
                report += """\n
            Information on genes with alternate splicing:
            Number of genes with alternate splicing = %d
            Alternate splicing genes:
            genename: [transcript name (number of exons)]...
                """ % len(self.alternate_splicing)
                for genename, alternate_splicing_info in self.alternate_splicing.items():
                    report += '%s: %s\n' % (genename, alternate_splicing_info)

            return report + '\n'

        elif self.rtype == ReportType.TEMP_REPORT:
            return "\nTEMP report! Nothing to report!\n"
        else:
            return "\nERROR: Report not initialized!\n"



    # DEPRECATED
    # Old toString function, not used any more!
    # def toString_OLD(self):
    #     if self.rtype == ReportType.FASTA_REPORT:
    #         report = """\n
    #         Reference format: FASTA
    #         General information:
    #         Reference length = %d bp
    #         Number of chromosomes = %d
    #         Chromosomes:\n%s
    #         Number of alignments (total / unique / non / real) = %d / %d / %d / %d
    #         Number of alignments (multi / split /possibly split) = %d / %d / %d
    #         Alignments with / without CIGAR string = %d / %d
    #         Mapping quality without zeroes (avg / min / max) = %.2f / %d / %d
    #         Alignments with mapping quality (>0 / =0) = %d / %d
    #         Number of matches / mismatches / inserts / deletes = %d / %d / %d / %d
    #         Percentage of matches / mismatches / inserts / deletes = %.2f / %.2f / %.2f / %.2f
    #         """ % (self.reflength, len(self.chromlengths), self.chromosomes(), \
    #                self.num_alignments, self.num_unique_alignments, self.num_non_alignments, self.num_real_alignments, \
    #                self.num_multi_alignments, self.num_split_alignments, self.num_possibly_split_alignements, \
    #                self.num_alignments - self.num_non_alignments, self.num_non_alignments, \
    #                self.avg_mapping_quality, self.min_mapping_quality, self.max_mapping_quality, self.num_good_quality, self.num_zero_quality, \
    #                self.num_match, self.num_mismatch, self.num_insert, self.num_delete, \
    #                self.match_percentage, self.mismatch_percentage, self.insert_percentage, self.delete_percentage)
    #         return report + '\n'
    #     elif self.rtype == ReportType.MAPPING_REPORT:
    #         report = """\n
    #         Command Line: %s
    #         Reference format: ANNOTATION
    #         General information:
    #         Reference length = %d bp
    #         Number of chromosomes = %d
    #         Chromosomes:
    #         %s
    #         Number of alignments in SAM file (total / unique) = %d / %d
    #         Number of alignments (multi / split /possibly split) = %d / %d / %d
    #         Number of alignments (real / non / real split) = %d / %d / %d
    #         Alignments with / without CIGAR string = %d / %d
    #         Mapping quality without zeroes (avg / min / max) = %.2f / %d / %d
    #         Alignments with mapping quality (>0 / =0) = %d / %d
    #         Number of matches / mismatches / inserts / deletes = %d / %d / %d / %d
    #         Percentage of matches / mismatches / inserts / deletes = %.2f / %.2f / %.2f / %.2f
    #         """ % (self.commandline, self.reflength, len(self.chromlengths), self.chromosomes(), \
    #                self.num_alignments, self.num_unique_alignments, \
    #                self.num_multi_alignments, self.num_split_alignments, self.num_possibly_split_alignements, \
    #                self.num_real_alignments, self.num_non_alignments, self.num_real_split_alignments, \
    #                self.num_alignments - self.num_non_alignments, self.num_non_alignments, \
    #                self.avg_mapping_quality, self.min_mapping_quality, self.max_mapping_quality, self.num_good_quality, self.num_zero_quality, \
    #                self.num_match, self.num_mismatch, self.num_insert, self.num_delete, \
    #                self.match_percentage, self.mismatch_percentage, self.insert_percentage, self.delete_percentage)
    #
    #         report += """\n
    #         Annotation statistics:
    #         Total gene length = %d
    #         Number of Transcripts / Exons (Multiexon transcripts) = %d / %d (%d)
    #         Maximum number of exons in a gene = %d
    #         Gene size (Min / Max / Avg) = %d / %d / %.2f
    #         Exon size (Min / Max / Avg) = %d / %d / %.2f
    #         """ % (self.totalGeneLength, self.num_genes, self.num_exons, self.num_multiexon_genes, self.max_exons_per_gene, \
    #                self.min_gene_length, self.max_gene_length, self.avg_gene_length, \
    #                self.min_exon_length, self.max_exon_length, self.avg_exon_length)
    #
    #         report += """\n
    #         Grouped annotation (alternate splicing) statistics:
    #         Number of annotation groups (genes) = %d
    #         Number of genes with alternate splicing = %d
    #         Maximum / minimum number of alternate spliced alignments for a gene = %d / %d
    #         Maximum / minimum number of exons in spliced alignments = %d / %d
    #         """ % (self.num_annotation_groups, self.num_alternate_spliced_genes, \
    #                self.max_spliced_alignments, self.min_spliced_alignments, \
    #                self.max_spliced_exons, self.min_spliced_exons)
    #
    #         report += """\n
    #         Mapping quality information:
    #         Bases in reads (aligned / total) (percent) = (%d / %d) (%.2f%%)
    #         Genes covered / missed / total = %d / %d / %d
    #         Exons covered / missed / total = %d / %d / %d
    #         Alignments on genes hit / partial / missed = %d / %d / %d
    #         Alignments on exons hit / partial / missed = %d / %d / %d
    #         Alignments spanning multiple genes = %d
    #         Alignments spanning multiple exons = %d
    #         Alignments with hit on gene and miss on exons (inside miss) = %d
    #         Split alignments with exon miss (or partial hit) = %d
    #         Alignments covering exons (all / some / none / multi) = %d / %d / %d / %d
    #         Alignments covering an exon (completely / partially) = %d / %d
    #         Alignments hitting an exon (start / end / both) = %d / %d / %d
    #         Candidate spliced alignments = %d
    #         Avg. alignment hit percentage = %.2f
    #         Avg. exon hit percentage = %.2f
    #         SUMMARY: Good / bad alignments: %d (%.2f%%) / %d (%.2f%%)
    #         """ % (self.sum_bases_aligned, self.sum_read_length, self.percentage_bases_aligned, \
    #                self.num_genes_covered, self.num_genes - self.num_genes_covered, self.num_genes, \
    #                self.num_exons_covered, self.num_exons - self.num_exons_covered, self.num_exons, \
    #                self.num_hit_alignments, self.num_partial_alignments, self.num_missed_alignments, \
    #                self.num_exon_hit, self.num_exon_partial, self.num_exon_miss, \
    #                self.num_multi_gene_alignments, self.num_multi_exon_alignments, self.num_inside_miss_alignments, self.num_bad_split_alignments, \
    #                self.num_cover_all_exons, self.num_cover_some_exons, self.num_cover_no_exons, self.num_multicover_exons, \
    #                self.num_equal_exons, self.num_partial_exons, \
    #                self.num_good_starts, self.num_good_ends, self.num_equal_exons, \
    #                self.num_possible_spliced_alignment, \
    #                self.avg_exon_hit_perc, self.avg_alignment_hit_perc, \
    #                self.num_good_alignment, self.good_alignment_percent, self.num_bad_alignment, self.bad_alignment_percent)
    #
    #         # Counting the number of expressed genes
    #         # This has already been written in the report, but this was the value can be double checked
    #         exp_gn_cnt = 0
    #         for expression in self.expressed_genes.values():
    #             if expression[0] > 0:
    #                 exp_gn_cnt += 1
    #
    #         if self.output_gene_expression:
    #             report += """\n
    #         Gene/exon expression and coverage information:
    #         Number of expressed genes = %d
    #         Expressed genes:
    #         genename  number_of_exons  gene_hits / gene_covered_bases  [(exon_hits / exon_covered_bases)]...
    #             """ % exp_gn_cnt
    #
    #             if len(self.expressed_genes) != len(self.gene_coverage):
    #                 raise Exception('ERROR: Gene expression and gene coverage dictionaries do not match in length! (%d <> %d)' % (len(expressed_genes), len(gene_coverage)))
    #
    #             for genename in self.expressed_genes.keys():
    #                 numexons = len(self.expressed_genes[genename]) - 1
    #                 genehits = self.expressed_genes[genename][0]
    #                 gene_cov_bs = self.gene_coverage[genename][0]
    #                 if genehits > 0:
    #                     reportline = '%s  %d  %d  %d' % (genename, numexons, genehits, gene_cov_bs)
    #                     for i in range(1, numexons+1):
    #                         exonhits = self.expressed_genes[genename][i]
    #                         exon_cov_bs = self.gene_coverage[genename][i]
    #                         reportline += '  (%d / %d)' % (exonhits, exon_cov_bs)
    #
    #                     report += reportline + '\n'
    #
    #         return report + '\n'
    #     elif self.rtype == ReportType.ANNOTATION_REPORT:
    #         report = """\n
    #         Command Line: %s
    #         Reference format: ANNOTATION
    #         General information:
    #         Reference length = %d bp
    #         Number of chromosomes = %d
    #         Chromosomes:
    #         %s
    #         """ % (self.commandline, self.reflength, len(self.chromlengths), self.chromosomes())
    #
    #         report += """\n
    #         Annotation statistics:
    #         Total gene length = %d
    #         Number of Genes / Exons (Multiexon genes) = %d / %d (%d)
    #         Maximum number of exons in a gene = %d
    #         Gene size (Min / Max / Avg) = %d / %d / %.2f
    #         Exon size (Min / Max / Avg) = %d / %d / %.2f
    #         """ % (self.totalGeneLength, self.num_genes, self.num_exons, self.num_multiexon_genes, self.max_exons_per_gene, \
    #                self.min_gene_length, self.max_gene_length, self.avg_gene_length, \
    #                self.min_exon_length, self.max_exon_length, self.avg_exon_length)
    #
    #         report += """\n
    #         Grouped annotation (alternate splicing) statistics:
    #         Number of annotation groups (genes) = %d
    #         Number of genes with alternate splicing = %d
    #         Maximum / minimum number of alternate spliced alignments for a gene = %d / %d
    #         Maximum / minimum number of exons in spliced alignments = %d / %d
    #         """ % (self.num_annotation_groups, self.num_alternate_spliced_genes, \
    #                self.max_spliced_alignments, self.min_spliced_alignments, \
    #                self.max_spliced_exons, self.min_spliced_exons)
    #
    #         if self.output_alternate_splicing:
    #             report += """\n
    #         Information on genes with alternate splicing:
    #         Number of genes with alternate splicing = %d
    #         Alternate splicing genes:
    #         genename: [transcript name (number of exons)]...
    #             """ % len(self.alternate_splicing)
    #
    #             for genename, alternate_splicing_info in self.alternate_splicing.items():
    #                 report += '%s: %s\n' % (genename, alternate_splicing_info)
    #
    #         return report + '\n'
    #     elif self.rtype == ReportType.TEMP_REPORT:
    #         return "\nTEMP report! Nothing to report!\n"
    #     else:
    #         return "\nERROR: Report not initialized!\n"


if __name__ == "__main__":
	pass;
