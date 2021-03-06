# Benchmark of RNAseq tools on third generation sequencing data
Scripts in RNAseqEval reposiroty were used for evaluating the performance of RNA mapping tools on third generation sequencing data. The tools were tested on several real and several synthetic datasets. Synthetic datasets are easier to evaluate since the exact origin of each read is known, and the alignment generated by a mapping tool can be compared to it. However, synthetic datasets fail to mimic every aspect of real data potentially biasing the evaluation. Therefore, our benchmark included both types of data, and readers can judge the results for themselves.

All datasets, references, gene annotations, gene expressions data and simulation data are available on FigShare.
- Reference, annotation and expression files - https://figshare.com/articles/RNA_benchmark_references_zip/5360950
- Test datasets - https://figshare.com/articles/RNA_benchmark_datasets/5360998
- Simulation data generated by PBSIM - https://figshare.com/articles/Simulation_data_for_RNA_benchmark/5361013

## Evaluating real datasets
Real datasets were evaluated using the RNAseqEval.py script with the option eval-maping and annotations passed to the script as a parameter. Detailed description on the working of RNAseqEval.py script can be found at [doc/RNAseqEval.md](doc/RNAseqEval.md).

Usage example:

    RNAseqEval.py eval-mapping dmelanogaster_genome.fasta real_dataset_1.sam -a dmlanogaster_annotations.gtf

Additionally, match rate is calculated using errorrates.py script from https://github.com/isovic/samscripts, by comparing alignments from the SAM file to the reference genome.

Usage example:

    errorrates.py base dmelanogaster_genome.fasta real_dataset_1.sam

## Evaluating synthetic datasets
Synthetic datasets were evaluated using Process_pbsim_data.py script which compares the SAM file with alignments to the simulation data generated by PBSIM. Since datasets ued in testing can be constructed from several simulations, data generated by PBSIM for each simulation must be organized in a specific way. Required PBSIM data organization and detailed description of the way Process_pbsim_data.py works is given in [doc/Process_pbsim_data.md](doc/Process_pbsim_data.md).

Usage example:

    Process_pbsim_data.py process dataset1_simulation_folder dataset1.sam scerevisiae_annotations.gtf

As with real datasets match rate is alculated using errorrates.py script from https://github.com/isovic/samscripts, by comparing alignments from the SAM file to the reference genome.

## Resource usage measurement
Resource consumption is measured by a fork of cgmmetime tool avaiable at https://github.com/isovic/cgmemtime.git.

