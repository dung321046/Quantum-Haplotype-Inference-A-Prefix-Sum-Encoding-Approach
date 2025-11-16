from collections import defaultdict

from utils import *

pen1 = 2
pen4 = 2
pen3 = 2


def formulate(gsize, hsize, hpairs):
    constraint_set = []
    OFFSET = 0
    vname2idx = dict()
    Q = defaultdict(int)

    for i in range(hsize):
        idx = len(vname2idx)
        vname2idx['x_{' + str(i) + '}'] = idx
        Q[(idx, idx)] += 1
    for i in range(gsize):
        for (u, v) in hpairs[i]:
            name = p_name(u, v)
            vname2idx[name] = len(vname2idx)
            x_u, x_v, x_uv = vname2idx[x_name(u)], vname2idx[x_name(v)], vname2idx[name]
            Q[(x_uv, x_uv)] += 2 * pen1
            Q[(x_uv, x_u)] -= pen1
            Q[(x_uv, x_v)] -= pen1
            constraint_set.append(["x,y>=z", x_name(u), x_name(v), name])

    for K in range(gsize):
        if len(hpairs[K]) == 1:
            OFFSET += pen4
            u, v = hpairs[K][0]
            x_u, x_v = vname2idx[x_name(u)], vname2idx[x_name(v)]
            Q[(x_u, x_v)] -= pen4
            continue
        k = len(hpairs[K])
        for i in range(k - 1):
            z = len(vname2idx)
            namez = 'a^' + str(K + 1) + "_" + str(i + 1)
            vname2idx[namez] = z
            namey = p_name(hpairs[K][i + 1][0], hpairs[K][i + 1][1])
            y = vname2idx[namey]
            if i == 0:
                namex = p_name(hpairs[K][i][0], hpairs[K][i][1])
            else:
                namex = 'a^' + str(K + 1) + "_" + str(i)
            x = vname2idx[namex]
            add_oplus(Q, x, y, z, pen3)
            constraint_set.append(["xVy=z", namex, namey, namez])
        OFFSET += pen4
        qkm1 = vname2idx['a^' + str(K + 1) + "_" + str(k - 1)]
        constraint_set.append(["x=1", 'a^' + str(K + 1) + "_" + str(k - 1)])
        Q[(qkm1, qkm1)] -= pen4
    return OFFSET, Q, vname2idx, constraint_set


if __name__ == "__main__":
    files = ['HapSet/Geno4Len5.txt']
    for filename in files:
        sys.stdout = sys.__stdout__
        gsize, glen, g, hsize, h = read_data(filename)
        hpairs = create_pairs_fast(h, gsize, hsize, glen, g)
        OFFSET, Q, vname2idx, constraint_set = formulate(gsize, hsize, hpairs)
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
        sys.stdout.close()
