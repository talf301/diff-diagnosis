import sys
import omim
import hpo_lib
import math
import os
import numpy
import sqlite3
import scipy.spatial.distance

from collections import defaultdict

#def calc_pheno_sim(p1, p2, ic_dict):
#    #p1 and p2 are HP codes - i.e. 0004322, ic_dict is dict from HP code to information content
#    return calc_mica(p1, p2, ic_dict)

def calc_dis_ic(dis, ic_dict): # Gets total ic for Disease dis, based on ic in ic_dict. ic_dict is a dict of (HP code, HP code): float.
    pheno = dis.phenotype_freqs
    ic = 0
    for p in pheno:
        if p:
            if not hpo_lib.search_code(p[3:]): continue
            ic += ic_dict[p[3:]]
    return ic 

def calc_mica(p1, p2, ic_dict): # Return ic of most informative common ancestor of p1, p2 (HP codes), based on ic in ic_dict
    p1anc = add_phenos(['HP:' + p1])
    p2anc = add_phenos(['HP:' + p2])
    common = p1anc.intersection(p2anc)
    if len(common) > 0:
        return max([ic_dict[p] for p in common])
    return 0

def calc_pheno_sim_db_cached(d1, d2, mica_dir): # Return similarity between Diseases d1, d2, by mica. mica values are cached in sqlite databases in mica_dir.
    #d1, d2 are disease objects
    pheno1 = d1.phenotype_freqs
    pheno2 = d2.phenotype_freqs
    ttl = 0
    if len(pheno2) > len(pheno1): # more efficient if query is shorter
        q = pheno1
        d = pheno2
    else:
        q = pheno2
        d = pheno1
    for p1 in q:
        dbname = '{}/{}.db'.format(mica_dir, p1[3:])
        if not os.path.isfile(dbname): continue

        conn = sqlite3.connect(dbname)
        c = conn.cursor()
        max_ic = 0
        for num in c.execute('SELECT ic FROM mica WHERE ' + ' OR '.join(['phe == "{}"'.format(x[3:]) for x in d])):
            if num[0] > max_ic: max_ic = num[0]

        ttl += max_ic
        conn.close()
    try:
        return (float(ttl)) 
    except ZeroDivisionError as e:
        return 0

def calc_pheno_sim_db_cached_symmetric(d1, d2, mica_dir):
    sim1 = calc_dis_sim_db_cached(d1, d2, mica_dir)
    sim2 = calc_dis_sim_db_cached(d2, d1, mica_dir)
    return 0.5 * sim1 + 0.5 * sim2

def calc_pheno_sim_cached(d1, d2, mica_dir): # Return similarity of Diseases d1, d2, by mica. mica is cached in txt files in mica_dir.
    #d1, d2 are disease objects
    pheno1 = d1.phenotype_freqs
    pheno2 = d2.phenotype_freqs
    ttl = 0
    norm = 0 #to normalize - perfect match will get 1
    all_hps = os.listdir(mica_dir)
    for p1 in pheno1:
        if not p1[3:] + '.txt' in all_hps: continue
        p1_file = open('{}/{}.txt'.format(mica_dir, p1[3:]))
        p1_micas = {info[1]: float(info[2]) for info in [line.split() for line in p1_file.readlines()]}
        p1_file.close()
        max = 0
        for p2 in pheno2:
            if p2[3:] in p1_micas:
                sim = p1_micas[p2[3:]]
                if sim > max:
                    max = sim
        ttl += max
        norm += p1_micas[p1[3:]]
    try:
        return (float(ttl) / norm) 
    except ZeroDivisionError as e:
        return 0

def calc_pheno_sim_cached_symmetric(d1, d2, mica_dir):
    sim1 = calc_dis_sim_cached(d1, d2, mica_dir)
    sim2 = calc_dis_sim_cached(d2, d1, mica_dir)
    return 0.5 * sim1 + 0.5 * sim2
    
    

def calc_dis_sim(d1, d2, ic_dict): # Return similarity of Diseases d1, d2 by mica. No caching.
    #d1, d2 are disease objects
    pheno1 = d1.phenotype_freqs
    pheno2 = d2.phenotype_freqs
    ttl = 0
    norm = 0 #to normalize - perfect match will get 1
    for p1 in pheno1:
        hp = hpo_lib.search_code(p1[3:])
        if not hp:
            continue;
        max = 0
        for p2 in pheno2:
            hp = hpo_lib.search_code(p2[3:])
            if not hp:
                continue;
            sim = calc_mica(p1[3:], p2[3:], ic_dict)
            if sim > max:
                max = sim
        ttl += max
        norm += ic_dict[p1[3:]]
    try:
        return (float(ttl) / norm) 
    except ZeroDivisionError as e:
        return 0



def calc_dis_sim_no_norm(d1, d2, ic_dict): # Same as calc_dis_sim but not normalized.
    #d1, d2 are disease objects
    pheno1 = d1.phenotype_freqs
    pheno2 = d2.phenotype_freqs
    ttl = 0
    norm = 0 #to normalize - perfect match will get 1
    mica_fxn = calc_mica
    search_fn = hpo_lib.search_code
    filtpheno1 = list(filter(lambda p1: search_fn(p1[3:]) != None, pheno1))
    filtpheno2 = list(filter(lambda p2: search_fn(p2[3:]) != None, pheno2))
    for p1 in filtpheno1:
        maxi = 0
        for p2 in filtpheno2:
            sim = mica_fxn(p1[3:], p2[3:], ic_dict)
            maxi = max(maxi, sim)
        ttl += maxi
        #norm += ic_dict[p1[3:]]
    try:
        return float(ttl)
    except ZeroDivisionError as e:
        return 0

def calc_dis_sim_symmetric(d1, d2, ic_dict):
    sim1 = calc_dis_sim(d1, d2, ic_dict)
    sim2 = calc_dis_sim(d2, d1, ic_dict)
    return 0.5 * sim1 + 0.5 * sim2

def calc_dis_sim_symmetric_no_norm(d1, d2, ic_dict):
    sim1 = calc_dis_sim_no_norm(d1, d2, ic_dict)
    sim2 = calc_dis_sim_no_norm(d2, d1, ic_dict)
    return 0.5 * sim1 + 0.5 * sim2

def calc_simgic(d1, d2, ic_dict, ancestors=True, norm=True): # Return similarity between Diseases d1, d2 with simgic metric by ic in ic_dict.
    pheno1 = d1.phenotype_freqs
    pheno2 = d2.phenotype_freqs
    allp1 = add_phenos(pheno1, ancestors)
    allp2 = add_phenos(pheno2, ancestors)
    common = allp1.intersection(allp2)
    all = allp1.union(allp2)
    #print([ic_dict[p] for p in common])

    common_ic = sum([ic_dict[p] for p in common])
    all_ic = sum([ic_dict[p] for p in all])
    if norm == True:
        try:
            return common_ic / all_ic
        except ZeroDivisionError as e:
            return 0
    elif norm == 'half':
        try:
            return common_ic / sum([ic_dict[p] for p in allp1])
        except ZeroDivisionError as e:
            return 0
    elif norm == False:
        try:
            #norm_sum =   sum([ic_dict[p] for p in allp1])
            #print(allp1, pheno1)
            #print(common_ic)
            return common_ic #/ norm_sum
        except ZeroDivisionError as e:
            return 0
    else:
        print('Invalid option: {}'.format(norm))
        sys.exit()

def calc_con_vec(d1, d2, vec_dict): # Return similarity between Diseases d1, d2 using context vectors in vec_dict - a dict of HP Code: HPO-length list of ints.
    pheno1 = d1.phenotype_freqs
    pheno2 = d2.phenotype_freqs
    allp1 = (add_phenos(pheno1)) # maybe we won't add ancestors - it would be redundant with what is already in the context vectors.
    allp2 = (add_phenos(pheno2))
    if len(allp1) == 0 or len(allp2) == 0: # if either patient has no symptoms - similarity will be 0
        return 0
    ttl1 = numpy.sum([vec_dict[v] for v in allp1], axis=0) # add up all context vectors of all HP terms.
    ttl2 = numpy.sum([vec_dict[v] for v in allp2], axis=0)
    try:
        sim = 1 - scipy.spatial.distance.cosine(ttl1, ttl2) # return cosine of the two diseases
        return sim
    except TypeError as e:
        print('ttl1: {}, ttl2: {}'.format(str(type(ttl1)), str(type(ttl2))))
        print(ttl1, ttl2)
        raise
    except ZeroDivisionError as e:
        return 0 
    

def load_vectors(f): # Load context vectors from file into dict. File is tab separated, one HP term per line. Each line has one number for each term in HPO.
    vec_dict = {}
    i = 0
    for line in f:
        if i % 1000 == 0 : print('{} vectors loaded'.format(i))
        i += 1
        info = line.split()
        vec_dict[info[0]] = [float(x) for x in info[1:]]
    return vec_dict

def parse_anno(f, names=False): # Parses annotation file. Returns dict of disease code: list of HP phenotypes.
    omim = {}
    omim_names = {}
    for line in f:
        info = line.strip().split('\t')
        dis_code = '{}:{}'.format(info[0], info[1])
        if not dis_code in omim:
            omim[dis_code] = []
        omim[dis_code].append(info[4])
        if names:
            if not dis_code in omim_names:
                omim_names[dis_code] = info[2]
    if names:
        return (omim, omim_names)
    else:
        return omim

def load_ic(f): # Return dict of HP code to info content (float)
    ic_dict = {}
    for line in f:
        info = line.strip().split()
        ic_dict[info[0]] = float(info[1])
    return ic_dict

def add_phenos(dis, ancestors=True): # dis is a list of phenotypes, return a set of all ancestors (including self) of all terms in dis.
    annotated = set()
    for phe in dis:
        phe = phe[3:]
        hp = hpo_lib.search_code(phe)
        if hp:
            annotated.add(phe)
            if not ancestors: continue  
            ancestors = hpo_lib.get_ancestors(hp)
            for anc in ancestors:
                annotated.add(anc.id[3:])           
    return annotated    

def get_p_value(score, score_list): # Return p-value of score, given score_list, an increasing list of floats, representing the null distribution.
    null_dist = score_list
    rank = binary_search(score, null_dist)
    return (len(null_dist) - rank + 1) / (len(null_dist) + 1)

def binary_search(item, L):
    if len(L) == 0:
        return 0
    if len(L) == 1:
        return int(item > float(L[0]))
    comp = float(L[len(L) // 2])
    if item <= comp:
        return binary_search(item, L[:len(L) // 2])
    else:
        return len(L) // 2 + binary_search(item, L[len(L) // 2:])

def get_dis_ic(dis, ic_dict):
    phenos = add_phenos(dis)
    return sum(ic_dict[hp] for hp in phenos)

if __name__ == '__main__':
    hpo_lib.build_decoder()
    print(len(hpo_lib.get_decoder()), len(hpo_lib.decoder))
    omim = parse_anno(open('/dupa-filer/talf/diff-diagnosis/phenotype_annotation.tab'))
    ttl_dis = len(omim)
    fout = open(sys.argv[1], 'w')
    freqs = defaultdict(lambda: 1) 
    for dis in omim:
        annotated = add_phenos(omim[dis])
        for a in annotated:
                freqs[a] += 1
    for phe in hpo_lib.decoder.keys():
        if hpo_lib.decoder[phe] != None:
            ic = -math.log(freqs[phe] / ttl_dis)
            fout.write('{}\t{}\n'.format(phe, str(ic)))
    print(len(hpo_lib.get_decoder()), len(hpo_lib.decoder))
    
    

