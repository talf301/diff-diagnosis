import os
import sys
import hpo_lib
import get_hp_ic
import omim
from argparse import ArgumentParser

__author__ = 'Tal Friedman (talf301@gmail.com)'

def script(path,ic, out):
    omim_dict = get_hp_ic.parse_anno(open('/dupa-filer/talf/diff-diagnosis/phenotype_annotation.tab'))
    #diseases = list(omim_dict.keys())
   
    ic_dict = get_hp_ic.load_ic(open(ic))
    contents = os.listdir(path)
    patients = [f for f in contents if f.endswith('.txt')]

    fail_counter = 0
    failed_keys = set()
    for i, p in enumerate(patients):
        if i % 10 == 0: print(i)
        try:
            phenotypes = open(os.path.join(path, p)).readlines()[0].split(',')
        except IndexError:
            continue
        omim_scores = []
        for o, phe in omim_dict.items():
            try:
                sim = get_hp_ic.calc_simgic(omim.Disease(None, None, None, phe), omim.Disease(None, None, None, phenotypes), ic_dict)
            except KeyError as e:
                fail_counter += 1
                failed_keys.add(e.args[0])

            omim_scores.append((o, sim))

        omim_scores.sort(key=lambda x: x[1],reverse=True)
        with open(os.path.join(out, p + '.results'), 'w') as f:
            for line in omim_scores[:20]:
                f.write('\t'.join(['1.000', str(line[1]), line[0]]) + '\n')

    print("Number of occurences of missing phenotypes: " + str(fail_counter) + '\n')
    print("Number of unique phenos failing: " + str(len(failed_keys)) + '\n')
    print(failed_keys)

def parse_args(args):
    parser = ArgumentParser(description='Get diagnosis based on simgic score')
    parser.add_argument('path', metavar='DIR')
    parser.add_argument('ic', metavar = 'IC')
    parser.add_argument('out', metavar='OUT')
    return parser.parse_args()

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())

