#!/usr/bin/env python

import os
import sys

from argparse import ArgumentParser
from omim import MIM

__author__ = 'Tal Friedman (talf301@gmail.com)'

def script(path, out):
    # Load omim
    mim = MIM(os.path.join(data_path, 'phenotype_annotation.tab'))
    omim = filter(lambda d:d.db == 'OMIM', mim.diseases)
    omim_dict = {dis.id:dis for dis in omim}

    contents = os.listdir(path)
    results = filter(lambda f: f.endswith('.results'), contents)

    for file in results:
        info = list(open(os.path.join(path, file)))
            

def parse_args(args):
    parser = ArgumentParser(description='Find the suggested phenotypes from the phenomizer diff')
    parser.add_argument('path', metavar='DIR')
    parser.add_argument('out', metavar='OUT')
    return parser.parse_args()

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
