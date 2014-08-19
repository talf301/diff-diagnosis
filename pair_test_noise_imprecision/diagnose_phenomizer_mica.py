import os
import sys
import hpo_lib
import get_hp_ic
import omim
from argparse import ArgumentParser
import time

__author__ = 'Tal Friedman (talf301@gmail.com)'

def script(path, out):
    omim_dict = get_hp_ic.parse_anno(open('/dupa-filer/talf/diff-diagnosis/phenotype_annotation.tab'))
    diseases = list(omim_dict.keys())
   
    ic_dict = get_hp_ic.load_ic(open('ic/ic_parent.txt'))
    contents = os.listdir(path)
    patients = [f for f in contents if f.endswith('_hpo.txt')]
    t1 = time.time()
    for i, p in enumerate(patients):
        print(i, time.time() - t1)
        try:
            phenotypes = open(os.path.join(path, p)).readlines()[0].split(',')
        except IndexError:
            continue
        omim_scores = []
        for o, phe in omim_dict.items():
            dist = get_hp_ic.calc_dis_sim_no_norm(omim.Disease(None, None, None, phe), omim.Disease(None, None, None, phenotypes), ic_dict)
            omim_scores.append((o, dist))

        omim_scores.sort(key=lambda x: x[1],reverse=True)
        with open(os.path.join(out, p + '.results'), 'w') as f:
            for line in omim_scores[:20]:
                f.write('\t'.join(['1.000', str(line[1]), 'OMIM:' + line[0]]) + '\n')

def parse_args(args):
    parser = ArgumentParser(description='Get diagnosis based on simgic score')
    parser.add_argument('path', metavar='DIR')
    parser.add_argument('out', metavar='OUT')
    return parser.parse_args()

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())

