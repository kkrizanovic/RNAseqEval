#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
import sys
sys.path.append(SCRIPT_PATH + '')
# sys.path.append(SCRIPT_PATH + '/wrappers')
import subprocess

import basicdefines

def create_folders():
    sys.stderr.write('Generating folders...\n')
    sys.stderr.write('\n')

    if not os.path.exists(basicdefines.TOOLS_ROOT_ABS):
        sys.stderr.write('Creating folder "%s".\n' % basicdefines.TOOLS_ROOT_ABS)
        os.makedirs(basicdefines.TOOLS_ROOT_ABS)

    if not os.path.exists(basicdefines.INTERMEDIATE_PATH_ROOT_ABS):
        sys.stderr.write('Creating folder "%s".\n' % basicdefines.INTERMEDIATE_PATH_ROOT_ABS)
        os.makedirs(basicdefines.INTERMEDIATE_PATH_ROOT_ABS)


def setup_tools():
    sys.stderr.write('Setting various tools up...\n')
    pass

def setup_all():
    create_folders()
    setup_tools()


def verbose_usage_and_exit():
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s\n' % sys.argv[0])
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    setup_all()
