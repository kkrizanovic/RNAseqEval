#! /usr/bin/python

import sys, os

# from file_compare import compare2files, compare3files


datasets = ['d1', 'd2', 'd3', 'd4']
measures = ['correct', 'hitall', 'hitone', 'incorrect', 'unmapped']
aligners = ['gm', 'gmap', 'mm2']


def compare2files(filename1, filename2):
    list1 = []
    list2 = []
    dict1 = {}
    dict2 = {}

    list1only = []
    list2only = []
    dictboth = {}       # Need dictionary to check for errors

    # Comparing two files!
    # Loading first file
    with open(filename1, 'rU') as file1:
        list1 = file1.readlines()
        file1.close()
    # Loading second file
    with open(filename2, 'rU') as file2:
        list2 = file2.readlines()
        file2.close()

    count1 = len(list1)
    count2 = len(list2)
    countboth = 0

    # Generating dictionary 1
    for line in list1:
        dict1[line] = 1
    # Generating dictionary 2
    for line in list2:
        dict2[line] = 1

    # Checking names from file 1
    for line in list1:
        if line in dict2:
            countboth += 1
            dictboth[line] = 1
        else:
            list1only.append(line)
    # Checking names from file 2
    for line in list2:
        if line in dict1:
            if line not in dictboth:
                sys.stderr.write('\nERROR with line: %s!' % line)
        else:
            list2only.append(line)

    return (list1only, list2only, dictboth, count1, count2, countboth)

def compare3files(filename1, filename2, filename3):
    list1 = []
    list2 = []
    list3 = []
    dict1 = {}
    dict2 = {}
    dict3 = {}

    list1only = []
    list2only = []
    list3only = []
    dictall = {}       # Need dictionary to check for errors

    # Comparing three files
    # Loading first file
    with open(filename1, 'rU') as file1:
        list1 = file1.readlines()
        file1.close()
    # Loading second file
    with open(filename2, 'rU') as file2:
        list2 = file2.readlines()
        file2.close()
    # Loading third file
    with open(filename3, 'rU') as file3:
        list3 = file3.readlines()
        file3.close()

    count1 = len(list1)
    count2 = len(list2)
    count3 = len(list3)
    countall = 0

    # Generating dictionary 1
    for line in list1:
        dict1[line] = 1
    # Generating dictionary 2
    for line in list2:
        dict2[line] = 1
    # Generating dictionary 3
    for line in list3:
        dict3[line] = 1

    # Checking names from file 1
    for line in list1:
        if line in dict2 and line in dict3:
            countall += 1
            dictall[line] = 1
        if line not in dict2 and line not in dict3:
            list1only.append(line)
    # Checking names from file 2
    for line in list2:
        if line in dict1 and line in dict3:
            if line not in dictall:
                sys.stderr.write('\nERROR with line: %s' % line)
        if line not in dict1 and line not in dict3:
            list2only.append(line)
    # Checking names from file 3
    for line in list3:
        if line in dict1 and line in dict2:
            if line not in dictall:
                sys.stderr.write('\nERROR with line: %s' % line)
        if line not in dict1 and line not in dict2:
            list3only.append(line)

    return (list1only, list2only, list3only, dictall, count1, count2, count3, countall)


def verbose_usage_and_exit():
    sys.stderr.write('qname_analysis - Analyze qnames for multiple aligners and measures.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s name_folder results_folder\n' % sys.argv[0])
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) != 3):
        verbose_usage_and_exit()

    namesfolder = sys.argv[1]
    resultsfolder = sys.argv[2]

    gm_qnames_filename = ''
    gmap_qnames_filename = ''
    mm2_qnames_filename = ''

    statsAll_filename = os.path.join(resultsfolder, 'statsAll.tsv')
    statsGMvsMM2_filename = os.path.join(resultsfolder, 'statsGMvsMM2.tsv')
    statsGMvsGMAP_filename = os.path.join(resultsfolder, 'statsGMvsGMAP.tsv')

    statsAllList = []
    statsGMvsMM2List = []
    statsGMvsGMAPList = []

    if not os.path.isdir(os.path.join(os.getcwd(), resultsfolder)):
        sys.stderr.write('\nCreating results folder and subfolders: %s' % os.path.join(os.getcwd(), resultsfolder))
        os.mkdir(os.path.join(os.getcwd(), resultsfolder))
    for dataset in datasets:
        if not os.path.isdir(os.path.join(os.getcwd(), resultsfolder, dataset)):
            os.mkdir(os.path.join(os.getcwd(), resultsfolder, dataset))
        if not os.path.isdir(os.path.join(os.getcwd(), resultsfolder, dataset, 'all')):
            os.mkdir(os.path.join(os.getcwd(), resultsfolder, dataset, 'all'))
        if not os.path.isdir(os.path.join(os.getcwd(), resultsfolder, dataset, 'gmvsmm2')):
            os.mkdir(os.path.join(os.getcwd(), resultsfolder, dataset, 'gmvsmm2'))
        if not os.path.isdir(os.path.join(os.getcwd(), resultsfolder, dataset, 'gmvsgmap')):
            os.mkdir(os.path.join(os.getcwd(), resultsfolder, dataset, 'gmvsgmap'))

    filesInOrder = True
    sys.stderr.write('\n\nStarting qname analysis.....')
    for dataset in datasets:
        sys.stderr.write('\n\nWorking with dataset: %s' % dataset)
        for measure in measures:
            gm_qnames_filename = os.path.join(namesfolder, 'gm', 'gm_%s_test_%s.names' % (dataset, measure))
            gmap_qnames_filename = os.path.join(namesfolder, 'gmap', 'gmap_%s_test_%s.names' % (dataset, measure))
            mm2_qnames_filename = os.path.join(namesfolder, 'mm2', 'mm2_%s_test_%s.names' % (dataset, measure))

            if not os.path.exists(gm_qnames_filename):
                sys.stderr.write('\nFile does not exist: %s' % gm_qnames_filename)
                filesInOrder = False
            if not os.path.exists(gmap_qnames_filename):
                sys.stderr.write('\nFile does not exist: %s' % gmap_qnames_filename)
                filesInOrder = False
            if not os.path.exists(mm2_qnames_filename):
                sys.stderr.write('\nFile does not exist: %s' % mm2_qnames_filename)
                filesInOrder = False

            # comparing all three aligners
            (gmonly, gmaponly, mm2only, dictall, gmcount, gmapcount, mm2count, countall) = compare3files(gm_qnames_filename, gmap_qnames_filename, mm2_qnames_filename)
            statsAllList.append((dataset, measure, gmcount, gmapcount, mm2count, countall))
            gmname = 'gm_%s_%s_only.names' % (dataset, measure)
            filename = os.path.join(os.getcwd(), resultsfolder, dataset, 'all', gmname)
            with open(filename,'w') as f:
                for line in gmonly:
                    f.write(line)
                f.close()
            gmapname = 'gmap_%s_%s_only.names' % (dataset, measure)
            filename = os.path.join(os.getcwd(), resultsfolder, dataset, 'all', gmapname)
            with open(filename,'w') as f:
                for line in gmaponly:
                    f.write(line)
                f.close()
            mm2name = 'mm2_%s_%s_only.names' % (dataset, measure)
            filename = os.path.join(os.getcwd(), resultsfolder, dataset, 'all', mm2name)
            with open(filename,'w') as f:
                for line in mm2only:
                    f.write(line)
                f.close()
            name = '%s_%s_all.names' % (dataset, measure)
            filename = os.path.join(os.getcwd(), resultsfolder, dataset, 'all', name)
            with open(filename,'w') as f:
                for line in dictall.iterkeys():
                    f.write(line)
                f.close()

            # comparing GraphMap and Minimap2
            (gmonly, mm2only, dictboth, gmcount, mm2count, countboth) = compare2files(gm_qnames_filename, mm2_qnames_filename)
            statsGMvsMM2List.append((dataset, measure, gmcount, mm2count, countboth))
            gmname = 'gm_%s_%s_only.names' % (dataset, measure)
            filename = os.path.join(os.getcwd(), resultsfolder, dataset, 'gmvsmm2', gmname)
            with open(filename,'w') as f:
                for line in gmonly:
                    f.write(line)
                f.close()
            mm2name = 'mm2_%s_%s_only.names' % (dataset, measure)
            filename = os.path.join(os.getcwd(), resultsfolder, dataset, 'gmvsmm2', mm2name)
            with open(filename,'w') as f:
                for line in mm2only:
                    f.write(line)
                f.close()
            name = 'gm_mm2_%s_%s_both.names' % (dataset, measure)
            filename = os.path.join(os.getcwd(), resultsfolder, dataset, 'gmvsmm2', name)
            with open(filename,'w') as f:
                for line in dictboth.iterkeys():
                    f.write(line)
                f.close()

            # comparing GraphMap and GMAP
            (gmonly, gmaponly, dictboth, gmcount, gmapcount, countboth) = compare2files(gm_qnames_filename, gmap_qnames_filename)
            statsGMvsGMAPList.append((dataset, measure, gmcount, gmapcount, countboth))
            gmname = 'gm_%s_%s_only.names' % (dataset, measure)
            filename = os.path.join(os.getcwd(), resultsfolder, dataset, 'gmvsgmap', gmname)
            with open(filename,'w') as f:
                for line in gmonly:
                    f.write(line)
                f.close()
            gmapname = 'gmap_%s_%s_only.names' % (dataset, measure)
            filename = os.path.join(os.getcwd(), resultsfolder, dataset, 'gmvsgmap', gmapname)
            with open(filename,'w') as f:
                for line in gmaponly:
                    f.write(line)
                f.close()
            name = 'gm_gmap_%s_%s_both.names' % (dataset, measure)
            filename = os.path.join(os.getcwd(), resultsfolder, dataset, 'gmvsgmap', name)
            with open(filename,'w') as f:
                for line in dictboth.iterkeys():
                    f.write(line)
                f.close()

    if filesInOrder:
        sys.stderr.write('\n\nAll files in order!')

    # Writting statistics
    with open(statsAll_filename,'w') as f:
        f.write('Dataset\tMeasure\tGraphmap\tGMAP\tMinimap2\tAll\n')
        for line in statsAllList:
            f.write('%s\t%s\t%d\t%d\t%d\t%d\t\n' % (line[0], line[1], line[2], line[3], line[4], line[5]))
        f.close()
    with open(statsGMvsMM2_filename,'w') as f:
        f.write('Dataset\tMeasure\tGraphmap\tMinimap2\tBoth\n')
        for line in statsGMvsMM2List:
            f.write('%s\t%s\t%d\t%d\t%d\t\n' % (line[0], line[1], line[2], line[3], line[4]))
        f.close()
    with open(statsGMvsGMAP_filename,'w') as f:
        f.write('Dataset\tMeasure\tGraphmap\GMAP\tBoth\n')
        for line in statsGMvsGMAPList:
            f.write('%s\t%s\t%d\t%d\t%d\t\n' % (line[0], line[1], line[2], line[3], line[4]))
        f.close()
