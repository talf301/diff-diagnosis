import hpo
import sys
import hpo_lib

HPO = hpo.script('hp.obo')
depths = {}
HPO_DICT = {}
decoder = {}
built = False
LCA_DICT = {}
parents = {}
grandparents = {}
ancestor_dict = {}

def build_decoder():
    get_codes(HPO.root)

def get_decoder():
    return decoder

def get_depth(code):
    hp = search_code(code)
    if hp:
        if hp in depths:
            return depths[hp]
        else:
            if not hp.parents and not hp._parent_hps:
                depth = 0
            else:
                p_depths = []
                for p in hp.parents:
                    p_depths.append(get_depth(p.id[3:]))
                depth = min(p_depths) + 1
            depths[hp] = depth
            return depth


def get_codes(root):
    if not root.id[3:] in decoder:
        decoder[root.id[3:]] = root
        for alt in root.alts:
            if not alt[3:] in decoder:
                decoder[alt[3:]] = root
        for c in root.children:
            get_codes(c)

def search_code(num, just_once=False):
    global built
    if not built and not just_once:
        build_decoder()
        built = True
    if num in decoder:
        return decoder[num]
    else:
        hp = search(num, HPO.root)
        decoder[num] = hp
        return hp

def get_parents_count(root):
    parents[root] = len(root.parents)
    for child in root.children:
        get_parents_count(child)

def get_grandparents_count(root):
    p = 0
    for par in root.parents:
        p += parents[par]
    grandparents[root] = p
    for child in root.children:
        get_grandparents_count(child)


def search(num, root):
    if root.id[3:] == num or num in [a[3:] for a in root.alts]:
        return root
    for c in root.children:
        s = search(num, c)
        if s:
            return s

def lca_distance(num1, num2, want_dist=True): 
    cutoff = 10000
    try:
        if (num1, num2) in LCA_DICT and want_dist:
            return LCA_DICT[(num1, num2)]
        else:
            hpa = search_code(num1)
            above_a = get_lineage(hpa) #a list of lists: all paths from hpa to root
            hpb = search_code(num2)
            above_b = get_lineage(hpb) #a list of lists: all paths from hpa to root
            (lca, dist) = find_lca(above_a, above_b, cutoff)
            LCA_DICT[(num1, num2)] = dist
            LCA_DICT[(num2, num1)] = dist
    except TypeError as e:
        print('ERROR on numbers:')
        print(num1, num2)
        sys.exit()
    if want_dist:
        return dist
    else:
        return lca

def lca(num1, num2):
    return lca_distance(num1, num2, False) 
    
def find_lca(a, b, cutoff=100000):#a and b are lineages
    global_min = 100000
    lca = None
    b_dists = {}
    a_dists = {}
    for path in b:
        for i in range(len(path)):
            if not i in b_dists:
                b_dists[i] = set()
            b_dists[i].add(path[i])
    for path in a:
        for i in range(len(path)):
            if not i in a_dists:
                a_dists[i] = set()
            a_dists[i].add(path[i])
    dist = 0
    while dist < max(list(a_dists.keys())) + max(list(b_dists.keys())) + 1:
        for i in range(0, dist + 1):
            j = dist - i
            if i in a_dists and j in b_dists:
                common = a_dists[i] & b_dists[j]
                if len(common) > 0:
                    return (common.pop(), dist)
        dist += 1


def get_ancestors(hp):
    if hp.id[3:] in ancestor_dict:
        return list(ancestor_dict[hp.id[3:]])
    ancestors = set()
    for p in hp.parents:
        ancestors = ancestors.union(get_ancestors(p))
    ancestors.add(hp)
    ancestor_dict[hp.id[3:]] = ancestors
    return list(ancestors)
    
def get_lineage(hp):
    if hp.id[3:] in HPO_DICT:
        return HPO_DICT[hp.id[3:]]
    if not hp.parents and not hp._parent_hps:
        lineage = [[]]
    else:
        lineage = []
        for p in hp.parents:#hp.parents.union(hp._parent_hps):
            p_line = get_lineage(p)
            for line in p_line:
                lineage.append(line)

    for i in range(len(lineage)):
        lineage[i] = [hp] + lineage[i]
    HPO_DICT[hp.id[3:]] = lineage
    for a in hp.alts:
        HPO_DICT[a[3:]] = lineage
    return lineage

if __name__ == '__main__':

    build_decoder()
    print(len(decoder))
    print(decoder['0000006'])
    print(decoder['0000001'])
    sys.exit()
    for code in codes:
        if not (sorted(get_ancestors_old(search_code(code))) == sorted(get_ancestors(search_code(code)))):
            print(code)
    sys.exit()
    num = input("what number?:")
    while num != '0':
        hp = search_code(num)
        if hp:
            print(hp.name)
            print(hp.alts)
            print(hp.id[3:])
            print(hp.parents)
            print(get_depth(hp.id[3:]))
        else:
            print('not found in HPO')
        num = input("what number?:")
    print("Thanks!")
    sys.exit()

    
