from collections import defaultdict

from utils import *


def get_values(selects, gsize, hsize, hpairs, vname2idx):
    state = dict()
    for i in range(hsize):
        if i in selects:
            state[vname2idx[x_name(i)]] = 1
        else:
            state[vname2idx[x_name(i)]] = 0
    for K in range(gsize):
        last = 0
        for i, (u, v) in enumerate(hpairs[K][:-1]):
            next_name = y_name(K + 1, i + 1)
            if u in selects and v in selects:
                state[vname2idx[next_name]] = 1
                last = 1
            else:
                state[vname2idx[next_name]] = last
    return state


def formulate(gsize, hsize, hpairs, A, B):
    constraint_set = []
    OFFSET = 0
    vname2idx = dict()
    Q = defaultdict(int)
    for i in range(hsize):
        vname2idx[x_name(i)] = i
        # Add objective
        Q[(i, i)] += 1
    for K in range(gsize):
        name = y_name(K + 1, 0)
        vname2idx[name] = len(vname2idx)
        for i, (u, v) in enumerate(hpairs[K]):
            next_name = y_name(K + 1, i + 1)
            vname2idx[next_name] = len(vname2idx)
            y_i, y_i1 = vname2idx[name], vname2idx[next_name]
            constraint_set.append(["x<=y", name, next_name])
            Q[(y_i, y_i)] += A
            Q[(y_i, y_i1)] -= A
            name = next_name
    for K in range(gsize):
        for i, (u, v) in enumerate(hpairs[K]):
            yk = vname2idx[y_name(K + 1, i)]
            yk1 = vname2idx[y_name(K + 1, i + 1)]
            xu = vname2idx[x_name(u)]
            xv = vname2idx[x_name(v)]
            # (y(k+1) - yk)(2 - xu - xv)
            Q[(yk1, yk1)] += 2 * B
            Q[(yk1, xu)] -= B
            Q[(yk1, xv)] -= B
            Q[(yk, yk)] -= 2 * B
            Q[(yk, xu)] += B
            Q[(yk, xv)] += B
            constraint_set.append(["yk1=yk+xu*xv", y_name(K + 1, i + 1),
                                   y_name(K + 1, i), x_name(u), x_name(v)])
    return OFFSET, Q, vname2idx, constraint_set


if __name__ == "__main__":
    import dimod
    import neal
    import time

    # files = ['Data/Geno4Len5.txt', 'Data/Geno6Len10.txt', 'Data/Geno12Len10.txt']
    files = ['Data/03Gen10Hap266.txt']
    for filename in files:
        sys.stdout = sys.__stdout__
        gsize, glen, g, hsize, h = read_data(filename)
        hpairs = create_pairs_fast(h, gsize, hsize, glen, g)
        OFFSET, Q, vname2idx, constraint_set = formulate(gsize, hsize, hpairs, 3, 1)
        bqm = dimod.BinaryQuadraticModel.from_qubo(Q, OFFSET)
        for K in range(gsize):
            bqm.fix_variable(vname2idx[y_name(K + 1, 0)], 0)
            bqm.fix_variable(vname2idx[y_name(K + 1, len(hpairs[K]))], 1)

        start = time.time()
        sampler = neal.SimulatedAnnealingSampler()
        print("#Var-#NoneZ:", len(vname2idx), '-', len(Q))
        # response = sampler.sample_qubo(Q, num_reads=100, sweeps=1000)

        initial = get_values(get_basic_sol(hpairs), gsize, hsize, hpairs, vname2idx)
        print(initial)
        print([(k, initial[v]) for k, v in vname2idx.items() if v in initial])
        response = sampler.sample(bqm, initial_states=initial, num_reads=1000, sweeps=100)

        best = 1000000
        best_response = []
        for run in response.data():
            if run.energy < best:
                best_response.append(run)
                best = run.energy
            elif run.energy == best:
                best_response.append(run)

        for rid, r in enumerate(best_response[:1]):
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
            check_answer_latex_form2(r, hsize, gsize, h, vname2idx, hpairs, g)
        sys.stdout.close()
        # sys.stdout = open(filename[:-3] + "out", 'w')
        sys.stdout = sys.__stdout__
        for r in best_response[:1]:
            objv = 0
            for i in range(hsize):
                if r.sample[vname2idx[x_name(i)]] == 1:
                    objv += 1
                    print(i, end=" ")
            print("\nObjective:", objv)
        # sys.stdout.close()
