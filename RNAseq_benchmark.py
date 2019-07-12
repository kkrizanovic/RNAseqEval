#! /usr/bin/python

# A dictionary used for RNAseq benchmark
SIM_FOLDER_DICT_ALL = {'SimG1' : 'group1'
                     , 'SimG2' : 'group2'
                     , 'SimG3' : 'group3'
                     , 'SimG1AS' : 'group1_AS'
                     , 'SimG1SS' : 'group1_SS'
                     , 'SimG2AS' : 'group2_AS'
                     , 'SimG2SS' : 'group2_SS'
                     , 'SimG3AS' : 'group3_AS'
                     , 'SimG3SS' : 'group3_SS'}


class benchmark_params:
    # A dictionary connecting fasta/fastq header prefix with the folder with PBSIM generated data
    # Containing information for reads with each prefix
    # This is used because data is simulated using several pbsim runs to get different
    # coverages for different sets of references (in this case transcripts)
    # NOTE: this should be changed for different simulations
    simFolderDict = SIM_FOLDER_DICT_ALL
