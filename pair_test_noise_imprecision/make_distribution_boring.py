import sys
import get_hp_ic
import random
import omim
import time
import hpo_lib

if __name__ == '__main__':
    dis_code = sys.argv[1]
    outdir = sys.argv[2]
    ic_file = sys.argv[3]
    dis_dict = get_hp_ic.parse_anno(open('/dupa-filer/talf/diff-diagnosis/phenotype_annotation.tab', encoding='utf-8'))
    pheno_dist = []
    for d in dis_dict:
        valid = filter(lambda hp: hpo_lib.search_code(hp[3:]) != None, dis_dict[d])
        pheno_dist += valid
    dis_phenos = dis_dict[dis_code]
    dis = omim.Disease(None, None, None, dis_phenos)
    score_fn = get_hp_ic.calc_simgic
    random_fn = random.choice
    omim_cons = omim.Disease
    ic_dict = get_hp_ic.load_ic(open(ic_file))
    t1 = time.time()
    for query_size in range(1, 11):
        scores = []
        for trial in range(100000):
            if trial % 1000 == 0: print(trial, time.time() - t1)
            query = [random_fn(pheno_dist) for i in range(query_size)]
            score = score_fn(omim_cons(None, None, None, query), dis, ic_dict, ancestors=True, norm=False)
            scores.append(score)
        scores.sort()
        fout = open('{}/{}.txt'.format(outdir, str(query_size)), 'w')
        for score in scores:
            fout.write(str(score) + '\n')
        fout.close()
        
