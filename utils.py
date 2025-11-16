import re
import sys
from sys import stderr
import os
import pandas as pd
import numpy as np
from collections import defaultdict

DEBUG = False


def get_basic_sol(hpairs):
    selects = set()
    for hp in hpairs:
        selects.add(hp[0][0])
        selects.add(hp[0][1])
    return list(selects)


def read_spp(filename):
    sys.stdin = open(filename, 'r')
    n, m = [int(x) for x in input().split()]
    c = [int(x) for x in input().split()]
    s = []
    for i in range(m):
        s.append([int(x) - 1 for x in input().split()])
    sys.stdin.close()
    return n, m, c, s


def read_ec(filename):
    sys.stdin = open(filename, 'r')
    n, m = [int(x) for x in input().split()]
    s = []
    for i in range(m):
        s.append([int(x) for x in input().split()])
    sys.stdin.close()
    return n, m, s


def read_cmic(filename):
    sys.stdin = open(filename, 'r')
    n, m = [int(x) for x in input().split()]
    p = [int(x) for x in input().split()]
    s = []
    for i in range(m):
        s.append([int(x) for x in input().split()])
    sys.stdin.close()
    return n, m, p, s


def read_data(filename):
    sys.stdin = open(filename, 'r')
    gsize, glen = input().split()
    gsize = int(gsize)
    glen = int(glen)
    g = []
    for i in range(gsize):
        g.append(input())
    input()
    input()
    # opt_sol = int(input())
    hsize = int(input())
    h = []
    for i in range(hsize):
        h.append(input())
    sys.stdin.close()
    return gsize, glen, g, hsize, h


def create_pairs(h, gsize, hsize, glen, g):
    hpairs = []
    for i in range(gsize):
        hpairs.append([])
    for i in range(hsize):
        for j in range(i):
            g_new = ''
            for k in range(glen):
                # print(h[i],h[j], k, glen)
                if h[i][k] == '1' and h[j][k] == '1':
                    g_new += '1'
                elif h[i][k] == '0' and h[j][k] == '0':
                    g_new += '0'
                else:
                    g_new += '2'
            if g_new in g:
                idx = g.index(g_new)
                hpairs[idx].append((i, j))
    return hpairs


def match(h, g):
    for i, hi in enumerate(h):
        if int(hi) + int(g[i]) == 1:
            return False
    return True


def create_pairs_fast(h, gsize, hsize, glen, g):
    hpairs = []
    hap_cans = []
    for i in range(gsize):
        hpairs.append([])
        hap_cans.append([])
    for i in range(gsize):
        for j in range(hsize):
            if match(h[j], g[i]):
                hap_cans[i].append(j)
    pair_set = set()
    for t in range(gsize):
        print(len(hap_cans[t]), end=" ", file=sys.stderr)
        for idx, j in enumerate(hap_cans[t]):
            for z in range(idx):
                i = hap_cans[t][z]
                if (i, j) not in pair_set:
                    pair_set.add((i, j))
                    g_new = ''
                    for k in range(glen):
                        if h[i][k] == '1' and h[j][k] == '1':
                            g_new += '1'
                        elif h[i][k] == '0' and h[j][k] == '0':
                            g_new += '0'
                        else:
                            g_new += '2'
                    if g_new in g:
                        idx = g.index(g_new)
                        hpairs[idx].append((i, j))
    print(file=sys.stderr)
    # for t in range(gsize):
    #     print(len(hpairs[t]), end="-", file=sys.stderr)
    # print(file=sys.stderr)
    return hpairs


def check(cc, pair2g, gens):
    for c in cc:
        i = 0
        prev = c
        ok = False
        while i < len(gens):
            ok = False
            for x in cc:
                if (prev, x) in pair2g and pair2g[(prev, x)] == gens[i]:
                    ok = True
                    next = prev
                    break
            if ok:
                prev = next
                i += 1
            else:
                break
        if ok:
            return True
    return False


def find_next(arr, x, pair2g, pairs):
    for pair in pairs:
        if (arr[pair[0]], x) not in pair2g:
            return False
        if pair2g[(arr[pair[0]], x)] != pair[1]:
            return False
    return True


def checkAll(cc, pair2g, nodes):
    for c in cc:
        i = 1
        arr = [c]
        ok = False
        while i < len(nodes):
            ok = False
            pairs = []
            for j in range(i):
                if (nodes[j], nodes[i]) in pair2g:
                    pairs.append((j, pair2g[(nodes[j], nodes[i])]))
            t = -1
            for x in cc:
                if find_next(arr, x, pair2g, pairs):
                    t = x
                    break
            if t == -1:
                break
            ok = True
            arr.append(t)
            i += 1
        if ok:
            # for i in range(len(nodes)):
            #     for j in range(i):
            #         if (nodes[j], nodes[i]) in pair2g:
            #             print(nodes[j], nodes[i], pair2g[(nodes[j], nodes[i])])
            #             print(arr[j], arr[i], pair2g[(arr[j], arr[i])])
            return True
    return False


def isSubset(nodes, ccs, pair2g):
    for cc in ccs:
        if checkAll(cc, pair2g, nodes):
            return True
    return False


def find_cc(hsize, gsize, hpairs):
    edges = []
    visit = []
    pair2g = dict()
    for i in range(hsize):
        edges.append([])
        visit.append(False)
    for K in range(gsize):
        for i, (u, v) in enumerate(hpairs[K]):
            edges[u].append((v, K))
            edges[v].append((u, K))
            pair2g[(u, v)] = K
            pair2g[(v, u)] = K
    ccs = []
    gsc = []
    for i in range(hsize):
        if not visit[i]:
            stack = [i]
            cc = [i]
            visit[i] = True
            gc = list()
            while len(stack) > 0:
                u = stack[0]
                stack = stack[1:]
                for (x, y) in edges[u]:
                    if not visit[x]:
                        gc.append(y)
                        stack.append(x)
                        visit[x] = True
                        cc.append(x)
            ccs.append(cc)
            gsc.append(gc)
    import numpy as np
    idx = np.argsort([-len(x) for x in ccs])
    min_ccs = []
    for id in idx:
        cc = ccs[id]
        gc = gsc[id]
        if len(cc) == len(gc) + 1:
            if not isSubset(cc, min_ccs, pair2g):
                # print(ccs[id])
                # print(gsc[id])
                min_ccs.append(cc)
            else:
                # if len(cc) > 4:
                #     print(ccs[id], "--")
                #     print(gsc[id], "--")
                # print("----Duplicate-----")
                pass
        else:
            # print(ccs[id])
            # print(gsc[id])
            min_ccs.append(cc)
    hpairs_new = []
    for i in range(gsize):
        hpairs_new.append([])
    for cc in min_ccs:
        for i in range(len(cc)):
            for j in range(i):
                if (cc[j], cc[i]) in pair2g:
                    hpairs_new[pair2g[(cc[j], cc[i])]].append((cc[j], cc[i]))
    return hpairs_new


# E = pen*(x^2 + y^2 + z^2 + xy -2xz 2-yz)
# E = 0 => min(1, x + y) = z
# z = 0 => x = y = 0 (x^2 + y^2 + xy)
# z = 1 => x^2 + y^2 + 1 + xy - 2x -2y
# x = 0 => y^2 + 1 - 2y =>y=1
# x = 1 => 1 + y^2 + 1 + y - 2 - 2y = y^2 - y
# E>0 => min(1, x + y) != z
# min(1, x + y) = z
def add_oplus(Q, x, y, z, pen):
    Q[(x, x)] += pen
    Q[(y, y)] += pen
    Q[(z, z)] += pen
    Q[(x, y)] += pen
    Q[(x, z)] -= 2 * pen
    Q[(y, z)] -= 2 * pen


def add_y_less_than_z(Q, y, z, pen):
    Q[(y, y)] += pen
    Q[(y, z)] -= pen


def add_xy0_z0(Q, x, y, z, pen):
    Q[(z, z)] += pen
    Q[(z, x)] -= pen
    Q[(z, y)] -= pen


def add_xy_less_than_2(Q, x, y, pen):
    Q[(x, y)] += pen


def add_less_than_2(Q, x, y, z, pen):
    Q[(x, x)] += pen
    Q[(y, y)] += pen
    Q[(z, z)] += pen
    Q[(x, y)] += 2 * pen
    Q[(x, z)] -= 2 * pen
    Q[(y, z)] -= 2 * pen


# E = (x - y)(1 - z)
def add_op3(Q, x, y, z, pen):
    Q[(x, x)] += pen
    Q[(y, y)] -= pen
    Q[(x, z)] -= pen
    Q[(y, z)] += pen


# (y - x_1 - x_2 - ... - x_n) ^ 2 = y + x_1 + ... + x_n - 2 y x_i + 2 x_i * x_j
def equal1(Q, y, xs, pen):
    Q[(y, y)] += pen
    n = len(xs)
    for i in range(n):
        Q[(xs[i], xs[i])] += pen
        Q[(y, xs[i])] -= 2 * pen
        for j in range(i):
            Q[(xs[i], xs[j])] += 2 * pen


def add_x_p_y_eq_z(Q, x, y, z, pen):
    Q[(z, z)] += pen
    Q[(x, x)] -= pen
    Q[(y, y)] -= pen


def equal1or2(Q, y, xs, pen):
    Q[(y, y)] += pen
    n = len(xs)
    for i in range(n):
        Q[(xs[i], xs[i])] += pen
        Q[(y, xs[i])] -= 2 * pen
        for j in range(i):
            Q[(xs[j], xs[i])] += pen


def add_2xor(Q, x, y, z, t, pen):
    Q[(x, x)] += pen
    Q[(y, y)] += pen
    Q[(z, z)] += pen
    Q[(t, t)] += pen
    Q[(x, t)] -= 2 * pen
    Q[(y, t)] -= 2 * pen
    Q[(z, t)] -= 2 * pen
    Q[(x, y)] += pen
    Q[(y, z)] += pen
    Q[(z, x)] += pen
    # Q[(x, y)] += 3 * pen
    # Q[(y, z)] += 3 * pen
    # Q[(z, x)] += 3 * pen


def p_name(u, v):
    return 'p_{' + str(v) + "," + str(u) + "}"


def s_name(u, v):
    return 's^{' + str(u) + "}_{" + str(v) + "}"


def y_name(u, v):
    return 'y^' + str(u) + "_" + str(v)


def x_name(u):
    return "x_{" + str(u) + "}"


def y1_name(u):
    return "y_{" + str(u) + "}"


def check_answer(r, vname2idx, hsize):
    for k, v in vname2idx.items():
        print(k, r.sample[vname2idx[k]])
    objv = 0
    for i in range(hsize):
        if r.sample[vname2idx[x_name(i)]] == 1:
            objv += 1
    print("Obj:", objv)


def check_answer_latex(r, hsize, gsize, h, vname2idx, hpairs, g):
    print(r"\begin{table}[ht]\begin{center}\begin{tabular}{ | c |c|c| } ")
    print(r"\hline Hanotype\newline Genotype & Variables & Value \\ \hline ")
    for i in range(hsize):
        name = x_name(i)
        print(h[i], "& $", name, "$ &", r.sample[vname2idx[name]], "\\\\")
    print("\hline ")
    for K in range(gsize):
        print(g[K])
        for (u, v) in hpairs[K]:
            name = p_name(u, v)
            print("& $", name, "$ &", r.sample[vname2idx[name]], "\\\\")
        print("\hline ")
    if hsize < 100:
        for k, v in vname2idx.items():
            if re.search("^a", k):
                print("& $", k, "$ &", r.sample[vname2idx[k]], "\\\\")
    objv = 0
    for i in range(hsize):
        if r.sample[vname2idx[x_name(i)]] == 1:
            objv += 1
    print("\hline  & Obj &", objv, " \\\\ \hline")
    print(r"\end{tabular}\end{center}")
    experiment_name = "G" + str(gsize) + "H" + str(hsize)
    print(r"\caption{\label{out:" + experiment_name + \
          "}Results of Simulated Annealing with baseline model on " + experiment_name + "}\end{table}")


def check_answer_latex_form2(r, hsize, gsize, h, vname2idx, hpairs, g):
    print(r"\begin{table}[ht]\begin{center}\begin{tabular}{ | c |c|c| } ")
    print(r"\hline Hanotype\newline Genotype & Variables & Value \\ \hline ")
    for i in range(hsize):
        name = x_name(i)
        print(h[i], "& $", name, "$ &", r.sample[vname2idx[name]], "\\\\")
    print("\hline ")
    for K in range(gsize):
        print(g[K])
        for i, (u, v) in enumerate(hpairs[K]):
            if i == 0:
                continue
            print("& $", y_name(K + 1, i), "$ &", r.sample[vname2idx[name]], "\\\\")
        print("\hline ")
    if hsize < 100:
        for k, v in vname2idx.items():
            if re.search("^y", k) and vname2idx[k] in r.sample:
                print("& $", k, "$ &", r.sample[vname2idx[k]], "\\\\")
    objv = 0
    for i in range(hsize):
        if r.sample[vname2idx[x_name(i)]] == 1:
            objv += 1
    print("\hline  & Obj &", objv, " \\\\ \hline")
    print(r"\end{tabular}\end{center}")
    experiment_name = "G" + str(gsize) + "H" + str(hsize)
    print(r"\caption{\label{out:" + experiment_name + \
          "}Results of Simulated Annealing with baseline model on " + experiment_name + "}\end{table}")


def check_constraint(constraint_set, r, vname2idx):
    numV = 0
    vioDict = defaultdict(lambda: 0)
    for constraint in constraint_set:
        if constraint[0] == "x=1":
            x = r.sample[vname2idx[constraint[1]]]
            if x != 1:
                if DEBUG:
                    print("Violate:", constraint[0], constraint[1], "=", x, file=stderr)
                numV += 1
                vioDict[constraint[0]] += 1
        elif constraint[0] == "Y=Sum2(X)":
            sum = 0
            t = 0
            for i in range(2, len(constraint)):
                sum += r.sample[vname2idx[constraint[i]]]
            if sum == 2:
                t = -1
            if sum + t != r.sample[vname2idx[constraint[1]]]:
                if DEBUG:
                    print("Violate:", constraint[0], constraint[1], r.sample[vname2idx[constraint[1]]], "!=", sum,
                          file=stderr)
                numV += 1
                vioDict[constraint[0]] += 1
        elif constraint[0] == "Y=Sum(X)":
            sum = 0
            for i in range(2, len(constraint)):
                sum += r.sample[vname2idx[constraint[i]]]
            if sum != r.sample[vname2idx[constraint[1]]]:
                if DEBUG:
                    print("Violate:", constraint[0], constraint[1], r.sample[vname2idx[constraint[1]]], "!=", sum,
                          file=stderr)
                numV += 1
                vioDict[constraint[0]] += 1
        elif constraint[0] == "xVy=1":
            x = r.sample[vname2idx[constraint[1]]]
            y = r.sample[vname2idx[constraint[2]]]
            if x == 0 and y == 0:
                if DEBUG:
                    print("Violate:", constraint[0], constraint[1], constraint[2], "=", x, y, file=stderr)
                numV += 1
                vioDict[constraint[0]] += 1
        elif constraint[0] == "x<=y":
            if vname2idx[constraint[1]] not in r.sample:
                x = 0
            else:
                x = r.sample[vname2idx[constraint[1]]]
            if vname2idx[constraint[2]] not in r.sample:
                y = 1
            else:
                y = r.sample[vname2idx[constraint[2]]]
            if x == 1 and y == 0:
                if DEBUG:
                    print("Violate:", constraint[0], constraint[1], constraint[2], "=", x, y, file=stderr)
                numV += 1
                vioDict[constraint[0]] += 1
        elif re.search("sum", constraint[0]):
            sum = 0
            for i in range(2, len(constraint)):
                sum += r.sample[vname2idx[constraint[i]]]
            if sum != constraint[1]:
                if DEBUG:
                    print("Violate:", constraint[0], constraint[1], "!=", sum, file=stderr)
                numV += 1
                vioDict[constraint[0]] += 1
        elif re.search("yk1=yk\+xu\*xv", constraint[0]):
            if vname2idx[constraint[1]] not in r.sample:
                yk1 = 1
            else:
                yk1 = r.sample[vname2idx[constraint[1]]]
            if vname2idx[constraint[2]] not in r.sample:
                yk = 0
            else:
                yk = r.sample[vname2idx[constraint[2]]]
            xu = r.sample[vname2idx[constraint[3]]]
            xv = r.sample[vname2idx[constraint[4]]]
            if yk1 - yk > min(xu, xv):
                if DEBUG:
                    print("Violate:", constraint[0], constraint[1], constraint[2], constraint[3], constraint[4], "!=",
                      yk1, yk, xu, xv, file=stderr)
                numV += 1
                vioDict[constraint[0]] += 1
        elif re.search("yk1=ykVxu=xv=0", constraint[0]):
            if vname2idx[constraint[1]] not in r.sample:
                yk1 = 1
            else:
                yk1 = r.sample[vname2idx[constraint[1]]]
            if vname2idx[constraint[2]] not in r.sample:
                yk = 0
            else:
                yk = r.sample[vname2idx[constraint[2]]]
            xu = r.sample[vname2idx[constraint[3]]]
            xv = r.sample[vname2idx[constraint[4]]]
            if yk1 != yk and min(xu, xv) != 0:
                if DEBUG:
                    print("Violate:", constraint[0], constraint[1], constraint[2], constraint[3], constraint[4], "!=",
                      yk1, yk, xu, xv, file=stderr)
                numV += 1
                vioDict[constraint[0]] += 1
        elif re.search("yk1-yk<=xi", constraint[0]):
            if vname2idx[constraint[1]] not in r.sample:
                yk1 = 1
            else:
                yk1 = r.sample[vname2idx[constraint[1]]]
            if vname2idx[constraint[2]] not in r.sample:
                yk = 0
            else:
                yk = r.sample[vname2idx[constraint[2]]]
            xu = r.sample[vname2idx[constraint[3]]]
            if yk1 != yk and xu == 0:
                if DEBUG:
                    print("Violate:", constraint[0], constraint[1], constraint[2], constraint[3], "!=",
                      yk1, yk, xu, file=stderr)
                numV += 1
                vioDict[constraint[0]] += 1
        elif constraint[0] == "x+y<=1":
            x = r.sample[vname2idx[constraint[1]]]
            y = r.sample[vname2idx[constraint[2]]]
            if x + y == 2:
                if DEBUG:
                    print("Violate:", constraint[0], constraint[1], constraint[2], "=", x, y, file=stderr)
                numV += 1
                vioDict[constraint[0]] += 1
        else:
            if DEBUG:
                print(constraint[0], file=stderr)
            x = r.sample[vname2idx[constraint[1]]]
            y = r.sample[vname2idx[constraint[2]]]
            if constraint[3] not in vname2idx:
                z = 1
            elif vname2idx[constraint[3]] not in r.sample:
                z = 1
            else:
                z = r.sample[vname2idx[constraint[3]]]
            if constraint[0] == "x*y=z":
                if x * y != z:
                    if DEBUG:
                        print("Violate:", constraint[0], constraint[1], constraint[2], constraint[3], "=", x, y, z,
                          file=stderr)
                    numV += 1
                    vioDict[constraint[0]] += 1
            elif constraint[0] == "x,y>=z":
                if not (x >= z and y >= z):
                    if DEBUG:
                        print("Violate:", constraint[0], constraint[1], constraint[2], constraint[3], "=", x, y, z,
                          file=stderr)
                    numV += 1
                    vioDict[constraint[0]] += 1
            elif constraint[0] == "x+y!=1=>z=0":
                if not (x + y != 1 and z == 1):
                    if DEBUG:
                        print("Violate:", constraint[0], constraint[1], constraint[2], constraint[3], "=", x, y, z,
                          file=stderr)
                    numV += 1
                    vioDict[constraint[0]] += 1
            elif constraint[0] == "xVy=z":
                if (x + y == 0 and z == 1) or (x + y >= 1 and z == 0):
                    if DEBUG:
                        print("Violate:", constraint[0], constraint[1], constraint[2], constraint[3], "=", x, y, z,
                          file=stderr)
                        print(vname2idx[constraint[1]], vname2idx[constraint[2]], vname2idx[constraint[3]], file=stderr)
                    numV += 1
                    vioDict[constraint[0]] += 1
            elif constraint[0] == "2>=x+y=z":
                if x + y != z:
                    if DEBUG:
                        print("Violate:", constraint[0], constraint[1], constraint[2], constraint[3], "=", x, y, z,
                          file=stderr)
                        print(vname2idx[constraint[1]], vname2idx[constraint[2]], vname2idx[constraint[3]], file=stderr)
                    numV += 1
                    vioDict[constraint[0]] += 1
            elif constraint[0] == "xTyTz=t":
                t = r.sample[vname2idx[constraint[4]]]
                if (x ^ y ^ z == 1 and t == 0) or (x + y + z == 0 and t == 1):
                    if DEBUG:
                        print("Violate:", constraint[0], constraint[1], constraint[2], constraint[3], constraint[4], "=", x,
                          y, z, t, file=stderr)
                    numV += 1
                    vioDict[constraint[0]] += 1
    return numV, vioDict


def check_exact_cover(n, m, cur_item_set, s):
    P = dict()
    for i in range(n):
        P[i] = []
    for i in range(m):
        for t in s[i]:
            P[t].append(i)
    for i in range(n):
        sum = 0
        for t in P[i]:
            sum += cur_item_set[t]
        if sum != 1:
            return 1
    return 0


def check_genotype(pairs, haps):
    for pair in pairs:
        if pair[0] in haps and pair[1] in haps:
            return
    for pair in pairs:
        if pair[0] in haps or pair[1] in haps:
            haps.add(pair[0])
            haps.add(pair[1])
            return
    haps.add(pairs[0][0])
    haps.add(pairs[0][1])
    return


def check_set_packing(states, set):
    # print(states, set)
    sum = 0
    for x in set:
        if states[x] == 1:
            sum += 1
        if sum > 1:
            sum = 1
            states[x] = 0
            # return False
    return True


def get_mean_std(values, col=1):
    vs = values[:][col]
    return np.mean(vs), np.std(vs)


def plotbox(x, y, name):
    names = []
    for i in range(len(x)):
        names.append(name + "-" + str(x[i]))
    np.set_printoptions(precision=3, suppress=True)
    y = np.asarray(y)
    # y = y.transpose()
    print(y)
    print(len(x), x)
    df = pd.DataFrame(data=y, columns=names)
    boxplot = df.boxplot()
