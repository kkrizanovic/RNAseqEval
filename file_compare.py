#! /usr/bin/python

import sys, os

def compare2files(filename1, filename2):
    list1 = []
    list2 = []
    dict1 = {}
    dict2 = {}

    list1only = []
    list2only = []
    dictboth = {}       # Need dictionary to check for errors

    sys.stdout.write('\nComparing two files!')

    sys.stdout.write('\nLoading first file!')
    with open(filename1, 'rU') as file1:
        list1 = file1.readlines()
        file1.close()

    sys.stdout.write('\nLoading second file!')
    with open(filename2, 'rU') as file2:
        list2 = file2.readlines()
        file2.close()

    count1 = len(list1)
    count2 = len(list2)
    countboth = 0

    sys.stderr.write('\nGenerating dictionary 1!')
    for line in list1:
        dict1[line] = 1

    sys.stderr.write('\nGenerating dictionary 2!')
    for line in list2:
        dict2[line] = 1

    sys.stderr.write('\nSanity check: list1 / dic1 || list2 / dict2: %d/%d || %d / %d' % (len(list1), len(dict1), len(list2), len(dict2)))

    sys.stderr.write('\nChecking names from file 1!')
    for line in list1:
        if line in dict2:
            countboth += 1
            dictboth[line] = 1
        else:
            list1only.append(line)

    sys.stderr.write('\nChecking names from file 2!')
    for line in list2:
        if line in dict1:
            if line not in dictboth:
                sys.stderr.write('\nERROR with line: %s!' % line)
        else:
            list2only.append(line)

    sys.stderr.write('\nFile1 / file2 / both: %d / %d / %d\n' % (count1, count2, countboth))

    sys.stderr.write('\nWritting names for file 1 only')
    with open('file1_only.names', 'w+') as file1:
        for line in list1only:
            file1.write(line)
        file1.close()

    sys.stderr.write('\nWritting names for file 2 only')
    with open('file2_only.names', 'w+') as file2:
        for line in list2only:
            file2.write(line)
        file2.close()

    sys.stderr.write('\nWritting names for both files\n')
    with open('file12_both.names', 'w+') as file_both:
        for line in dictboth.iterkeys():
            file_both.write(line)
        file_both.close()

    # Names that occur in only one file are stored in lists and can be easily retreived andsaved to files

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

    sys.stdout.write('\nComparing three files!')

    sys.stdout.write('\nLoading first file!')
    with open(filename1, 'rU') as file1:
        list1 = file1.readlines()
        file1.close()

    sys.stdout.write('\nLoading second file!')
    with open(filename2, 'rU') as file2:
        list2 = file2.readlines()
        file2.close()

    sys.stdout.write('\nLoading third file!')
    with open(filename3, 'rU') as file3:
        list3 = file3.readlines()
        file3.close()

    count1 = len(list1)
    count2 = len(list2)
    count3 = len(list3)
    countall = 0

    sys.stderr.write('\nGenerating dictionary 1!')
    for line in list1:
        dict1[line] = 1

    sys.stderr.write('\nGenerating dictionary 2!')
    for line in list2:
        dict2[line] = 1

    sys.stderr.write('\nGenerating dictionary 3!')
    for line in list3:
        dict3[line] = 1

    sys.stderr.write('\nSanity check: list1 / dict1 || list2 / dict2 || list3 / dict3: %d/%d || %d/%d || %d/%d' \
                      % (len(list1), len(dict1), len(list2), len(dict2), len(list3), len(dict3)))

    sys.stderr.write('\nChecking names from file 1!')
    for line in list1:
        if line in dict2 and line in dict3:
            countall += 1
            dictall[line] = 1
        if line not in dict2 and line not in dict3:
            list1only.append(line)

    sys.stderr.write('\nChecking names from file 2!')
    for line in list2:
        if line in dict1 and line in dict3:
            if line not in dictall:
                sys.stderr.write('\nERROR with line: %s' % line)
        if line not in dict1 and line not in dict3:
            list2only.append(line)

    sys.stderr.write('\nChecking names from file 3!')
    for line in list3:
        if line in dict1 and line in dict2:
            if line not in dictall:
                sys.stderr.write('\nERROR with line: %s' % line)
        if line not in dict1 and line not in dict2:
            list3only.append(line)

    sys.stderr.write('\nFile1 / file2 / file3 / all: %d / %d / %d / %d\n' % (count1, count2, count3, countall))
    sys.stderr.write('\nFile1only / file2only / file3only: %d / %d / %d\n' % (len(list1only), len(list2only), len(list3only)))

    # NOTE: At the moment not interested in printing lines out
    # for line in dictall.iterkeys():
    #     sys.stdout.write(line)

    # Names that occur in only one file are stored in lists and can be easily retreived andsaved to files

def verbose_usage_and_exit():
    sys.stderr.write('File compare - compare lines from two or three files.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s file1 file2 [file3]\n' % sys.argv[0])
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 3 or len(sys.argv) > 4):
        verbose_usage_and_exit()

    filename1 = sys.argv[1]
    filename2 = sys.argv[2]

    if len(sys.argv) == 3:
        compare2files(filename1, filename2)
    else:
        filename3 = sys.argv[3]
        compare3files(filename1, filename2, filename3)



    # END
