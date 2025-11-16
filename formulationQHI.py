from collections import defaultdict

from utils import *

pen1 = 2.0
pen3 = 2.0
pen4 = 2.0


def formulate(sizek, gsize, hsize, hpairs):
    constraint_set = []
    OFFSET = 0
    vname2idx = dict()
    Q = defaultdict(int)

    for i in range(hsize):
        idx = len(vname2idx)
        vname2idx['x_{' + str(i) + '}'] = idx
        # Add objective
        Q[(idx, idx)] += 1
    for i in range(gsize):
        for (u, v) in hpairs[i]:
            name = p_name(u, v)
            vname2idx[name] = len(vname2idx)
            x_u, x_v, x_uv = vname2idx[x_name(u)], vname2idx[x_name(v)], vname2idx[name]
            # Constraints: x_uv = 1 => x_u = x_v = 1
            # 2 p_uv - p_uv * p_u - p_uv * p_v
            Q[(x_uv, x_uv)] += 2 * pen1
            Q[(x_uv, x_u)] -= pen1
            Q[(x_uv, x_v)] -= pen1
            constraint_set.append(["x,y>=z", x_name(u), x_name(v), name])

    for K in range(gsize):
        k = len(hpairs[K])
        names_cur = [p_name(hpairs[K][i][0], hpairs[K][i][1]) for i in range(k)]
        count = 0
        while len(names_cur) >= sizek:
            cnames = ["Y=Sum2(X)"]
            idx = range(sizek)
            names = [names_cur[id] for id in idx]
            namez = 's^' + str(K + 1) + "_" + str(count)
            vname2idx[namez] = len(vname2idx)
            count += 1
            cnames.append(namez)
            for name in names:
                cnames.append(name)
                names_cur.remove(name)
            # names_cur.append(namez)
            names_cur.insert(0, namez)
            equal1or2(Q, vname2idx[namez], [vname2idx[name] for name in names], pen3)
            constraint_set.append(cnames)
        if len(names_cur) > 1:
            cnames = ["Y=Sum2(X)"]
            idx = range(len(names_cur))
            names = [names_cur[id] for id in idx]
            namez = 'a^' + str(K + 1) + "_" + str(count)
            vname2idx[namez] = len(vname2idx)
            count += 1
            cnames.append(namez)
            for name in names:
                cnames.append(name)
                names_cur.remove(name)
            # names_cur.append(namez)
            names_cur.insert(0, namez)
            equal1or2(Q, vname2idx[namez], [vname2idx[name] for name in names], pen3)
            # equal1(Q, vname2idx[namez], [vname2idx[name] for name in names], pen3)
            constraint_set.append(cnames)
        idx = []
        t = len(names_cur)
        names = ["sum" + str(t), 1]
        for name in names_cur:
            names.append(name)
            idx.append(vname2idx[name])
        constraint_set.append(names)
        # Formulation for  pen4 * (1 - idx[0] - ... - idx[t-1]) ^ 2
        # pen4 - pen4 * sum(idx[i]) + 2 * pen4 * sum(idx[i],idx[j])
        OFFSET += pen4
        for i in range(t):
            Q[(idx[i], idx[i])] -= pen4
            for j in range(i):
                Q[(idx[i], idx[j])] += 2 * pen4
    return OFFSET, Q, vname2idx, constraint_set


if __name__ == "__main__":
    files = ['HapSet/Geno4Len5.txt']
    for filename in files:
        ksize = 3
        gsize, glen, g, hsize, h = read_data(filename)
        hpairs = create_pairs(h, gsize, hsize, glen, g)
        OFFSET, Q, vname2idx, constraint_set = formulate(ksize, gsize, hsize, hpairs)
        import neal
        import time

        start = time.time()
        sampler = neal.SimulatedAnnealingSampler()
        print("#Var-#NoneZ:", len(vname2idx), '-', len(Q))
        response = sampler.sample_qubo(Q, num_reads=100, sweeps=1000)

        best = 1000000
        best_response = []
        for run in response.data():
            if run.energy < best:
                best_response.append(run)
                best = run.energy
            elif run.energy == best:
                best_response.append(run)

        for rid, r in enumerate(best_response[:10]):
            objv = 0
            for i in range(hsize):
                if r.sample[vname2idx[x_name(i)]] == 1:
                    objv += 1
            print("run:", rid, " Obj:", objv, "Energy:", r.energy + OFFSET)
            check_constraint(constraint_set, r, vname2idx)
        end = time.time()
        print("Time running", end - start)
        sys.stdout = open('output.tex', 'w')
        for r in best_response[:1]:
            check_answer_latex(r, hsize, gsize, h, vname2idx, hpairs, g)

        sys.stdout = open(filename[:-3] + "out", 'w')
        for r in best_response[:1]:
            objv = 0
            for i in range(hsize):
                if r.sample[vname2idx[x_name(i)]] == 1:
                    objv += 1
                    print("1")
                else:
                    print("0")
            print("Objective:", objv)
