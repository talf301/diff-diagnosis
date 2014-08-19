#!/usr/bin/env python


import os
import sys
import logging

from argparse import ArgumentParser
from orpha import Orphanet

def script(path):
    logging.basicConfig(filename=os.path.join(path,'score.log'), level=logging.INFO)
    orphanet_lookup = '/dupa-filer/talf/matchingsim/patients/orphanet_lookup.xml'
    orphanet_inher = '/dupa-filer/talf/matchingsim/patients/orphanet_inher.xml'
    orphanet_geno_pheno = '/dupa-filer/talf/matchingsim/patients/orphanet_geno_pheno.xml'
    orph = Orphanet(orphanet_lookup, orphanet_inher, orphanet_geno_pheno)
    lookup = orph.lookup

    contents = os.listdir(path)
    results = filter(lambda f: f.endswith('.results'), contents)
   
    #pvalues = []
    counter = 0
    topcounter = 0
    top5counter = 0
    for file in results:
        info = list(open(os.path.join(path, file)))
        orphanum = file.split('_')[-3]
        omim = lookup[orphanum].pheno[0]
        if info:
            counter += any(x.split('\t')[2].split(':')[1].strip() == omim for x in info)
            topcounter += any(x.split('\t')[2].split(':')[1].strip() == omim for x in info[:1])
            top5counter += any(x.split('\t')[2].split(':')[1].strip() == omim for x in info[:5])
        #print file
        #pvalues.append(float(info[0].split('\t')[0]))


    #pvalues.sort()
    logging.info("Total contained in results: %d" % counter)
    logging.info("Total top hits correct: %d" % topcounter)
    logging.info("Total top 5 correct: %d" % top5counter)
    #logging.info("Highest pvalue: %f" %pvalues[-1])
    #logging.info("Lowest pvalue: %f" % pvalues[0])
    #logging.info("Average pvalue: %f" % sum(pvalues)/float(len(pvalues)))

def parse_args(args):
    parser = ArgumentParser(description='Score how well a phenomizer results directory is doing')
    parser.add_argument('path', metavar='DIR') 
    return parser.parse_args()

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
