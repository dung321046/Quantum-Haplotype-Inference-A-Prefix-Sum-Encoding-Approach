import random

from utils import *
import os

random.seed(1)


def gen_g(h1, h2):
    ans = ''
    for i in range(len(h1)):
        if h1[i] == '1' and h2[i] == '1':
            ans += '1'
        elif h1[i] == '0' and h2[i] == '0':
            ans += '0'
        else:
            ans += '2'
    return ans


def create_haplotype(g):
    h1 = ''
    h2 = ''
    for gi in g:
        if gi == '0':
            h1 += '0'
            h2 += '0'
        elif gi == '1':
            h1 += '0'
            h2 += '0'
        else:
            if random.random() > 0.5:
                h1 += '0'
                h2 += '1'
            else:
                h1 += '1'
                h2 += '0'
    return h1, h2


def create_haplotype_from_2gen(g1, g2):
    h1 = ''
    h2 = ''
    hm = ''
    for i in range(len(g1)):
        p1 = int(g1[i])
        p2 = int(g2[i])
        if p1 + p2 == 1:
            return "", "", ""
        if p1 == p2 and p1 < 2:
            h1 += g1[i]
            h2 += g1[i]
            hm += g1[i]
        elif p1 == 2:
            if p2 < 2:
                h1 += str(1 - p2)
                h2 += g2[i]
                hm += g2[i]
            else:
                if random.random() > 0.5:
                    h1 += '0'
                    hm += '1'
                    h2 += '0'
                else:
                    h1 += '1'
                    h2 += '1'
                    hm += '0'
        elif p2 == 2 and p1 < 2:
            h1 += g1[i]
            h2 += str(1 - p1)
            hm += g1[i]
        else:
            print("ERRORRR", file=sys.stderr)
    return h1, h2, hm


def gen_genotype_set(density, hsample, rr):
    gs = set()
    for i in range(1, rr):
        for j in range(i):
            if random.random() > density:
                g = gen_g(hsample[i], hsample[j])
                gs.add(g)
    return list(gs)


def associate(h, g):
    for i, hi in enumerate(h):
        if int(g[i]) < 2 and g[i] != hi:
            return False
    return True


def find_h(h, g):
    h2 = ''
    for i, hi in enumerate(h):
        if int(g[i]) == 2:
            h2 += str(1 - int(hi))
        else:
            if g[i] != hi:
                return False, h2
            h2 += hi
    return True, h2


def single_relation(h, t, gs):
    for i in range(len(gs)):
        if i != t:
            b, ch = find_h(h, gs[i])
            if b:
                return False, ch
    return True, "Invalid"


def gen_all(hap, cur, seql):
    if len(cur) == seql:
        hap.add(cur)
        return
    gen_all(hap, cur + "0", seql)
    gen_all(hap, cur + "1", seql)


sys.stdin = open('HapSet/HapSet.inp', 'r')
folderName = "Data/"
TT = 0
while TT < 50:
    s = input()
    if not re.search("segsites", s):
        continue
    input()
    hsample = set()
    while True:
        try:
            s = input()
            if s == "":
                break
            hsample.add(s)
        except EOFError:
            break
    TT += 1
    hsample = list(hsample)
    hsample = list(dict.fromkeys(hsample))
    seql = len(hsample[0])
    init_size = len(hsample)
    sys.stdout = sys.__stdout__
    print(len(hsample))
    gs = set()
    cnt = 0
    gsize = random.randint(init_size, 2 * init_size - 1)
    hpair_set = set()
    while len(gs) < gsize:
        h1 = random.randint(0, len(hsample) - 1)
        h2 = random.randint(0, len(hsample) - 1)
        h1, h2 = min(h1, h2), max(h1, h2)
        if h1 == h2 or (h1, h2) in gs:
            continue
        gs.add(gen_g(hsample[h1], hsample[h2]))
    gs = list(gs)
    n = len(gs)
    hsample = set()
    gen_all(hsample, "", seql)
    hsample = list(hsample)
    hpairs = create_pairs_fast(hsample, n, len(hsample), seql, gs)
    selected = set()
    for K in range(gsize):
        for (u, v) in hpairs[K]:
            selected.add(u)
            selected.add(v)
    print(hpairs, file=sys.stderr)
    print(selected, file=sys.stderr)
    if not os.path.exists(folderName):
        os.mkdir(folderName)
    filename = folderName + str(TT).zfill(2) + "Gen" + str(n) + "Hap" + str(len(selected)) + ".txt"
    sys.stdout = open(filename, 'w')
    print(n, seql)
    for g in gs:
        print(g)
    print()
    print()
    print(len(selected))
    for s in selected:
        print(hsample[s])
    sys.stdout.close()