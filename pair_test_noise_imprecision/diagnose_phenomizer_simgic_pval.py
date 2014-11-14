import os
import sys
import hpo_lib
import get_hp_ic
import omim
import time
from argparse import ArgumentParser
import random
import numpy

__author__ = 'Tal Friedman (talf301@gmail.com)'

def script(path,ic, out, pvals):
    omim_dict = get_hp_ic.parse_anno(open('/dupa-filer/talf/diff-diagnosis/phenotype_annotation.tab'))
    #diseases = list(omim_dict.keys())
   
    ic_dict = get_hp_ic.load_ic(open(ic))
    contents = os.listdir(path)
    patients = [f for f in contents if f.endswith('.txt')]

    (pval_table, pval_dict) = load_pvals(pvals)

    fail_counter = 0
    failed_keys = set()
    t1 = time.time()
    for i, p in enumerate(patients):
        if i % 10 == 0: print(i, time.time() - t1)
        try:
            phenotypes = open(os.path.join(path, p)).readlines()[0].split(',')
        except IndexError:
            continue
        omim_scores = []
        for o, phe in omim_dict.items():
            if not o in pval_dict: continue
            if len(phenotypes) == 0:
                (pval, sim, final_score) = (1.0, 0, 0)
            else:
                pval_list = pval_table[pval_dict[o]][min(len(phenotypes) - 1, 9)]
                try:
                    sim = get_hp_ic.calc_simgic(omim.Disease(None, None, None, phe), omim.Disease(None, None, None, phenotypes), ic_dict)
                    pval = get_hp_ic.get_p_value(sim, pval_list) * 0.1
                    final_score = (1 - pval) * 1000000000000 + sim
                except KeyError as e:
                    fail_counter += 1
                    failed_keys.add(e.args[0])

            omim_scores.append((o, final_score, pval, sim))
            print(pval)
        omim_scores.sort(key=lambda x: x[1], reverse=True)
        with open(os.path.join(out, p + '.results'), 'w') as f:
            for line in omim_scores[:20]:
                f.write('\t'.join([str(line[2]), str(line[3]), line[0]]) + '\n')

    print("Number of occurences of missing phenotypes: " + str(fail_counter) + '\n')
    print("Number of unique phenos failing: " + str(len(failed_keys)) + '\n')
    print(failed_keys)

def load_pvals_old(topdir):
    t1 = time.time()
    pval_dict = {}
    disdirs = os.listdir(topdir)
    for j, disdir in enumerate(disdirs):
        if j % 10 == 0: print('loaded {} diseases'.format(j), time.time() - t1)
        #pval_dict[disdir] = {i:None for i in range(1, 11)}
        pval_dict[disdir] = {}
        for i in range(1, 11):
            pvals_file = open('{}/{}/{}.txt'.format(topdir, disdir, str(i)))
            pvals = [float(num.strip()) for num in pvals_file.readlines()]
            #print(sys.getsizeof(pvals), sys.getsizeof(pval_dict[disdir]))
            #pvals = [random.random() for num in range(100000)]
            pval_dict[disdir][i] = pvals
            pvals_file.close()
    return pval_dict

def load_pvals(topdir):
    cutoff = 100000
    t1 = time.time()
    pval_dict = {}
    disdirs = os.listdir(topdir)
    pval_table = numpy.empty((len(disdirs), 10, cutoff))
    for j, disdir in enumerate(disdirs):
        if j % 10 == 0: print('loaded {} diseases'.format(j), time.time() - t1)
        pval_dict[disdir] = j
        for i in range(1, 11):
            l = 0
            f = open('{}/{}/{}.txt'.format(topdir, disdir, str(i)))
            for line in f:
                pval_table[j][i - 1][l] = float(line.strip())
                l += 1
                if(l >= cutoff): break
            f.close()
    return (pval_table, pval_dict)
                

def parse_args(args):
    parser = ArgumentParser(description='Get diagnosis based on simgic score')
    parser.add_argument('path', metavar='DIR')
    parser.add_argument('ic', metavar = 'IC')
    parser.add_argument('out', metavar='OUT')
    parser.add_argument('pvals', metavar='PVAL_DIR')
    return parser.parse_args()

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())

