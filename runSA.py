import time

import dimod
import neal

import HIPP_models
from utils import *

folder_pre = "output/"
folder_path = 'DataF50'
csize = 2
pen1 = 1
pen2 = 1
pen3 = 1
# H is QHI^k method (with k=csize=2)
# B is SCP method
# D is PrefixHI method
# method = "H" + str(csize) + "-SA-"
method = "D" + "-SA-P" + str(pen1) + "_" + str(pen2) + "-"
# method = "B-SA" + "-"
num_reads = 1000
sweeps = 100
config = method + str(num_reads) + "--" + str(sweeps)
print(config, file=stderr)
if not os.path.exists(folder_pre):
    os.mkdir(folder_pre)
if not os.path.exists(folder_pre + config):
    os.mkdir(folder_pre + config)
TT = 0
for f in os.listdir(folder_path):
    sys.stdout = sys.__stdout__
    TT += 1
    filenameINP = os.path.join(folder_path, f)
    if not os.path.isfile(filenameINP):
        continue
    gsize, glen, g, hsize, h = read_data(filenameINP)

    hpairs = create_pairs(h, gsize, hsize, glen, g)
    total_hpairs = 0
    disHpairs = ""
    for k in range(gsize):
        total_hpairs += len(hpairs[k])
        disHpairs += str(len(hpairs[k])) + "-"
    OFFSET, Q, vname2idx, constraint_set = HIPP_models.load_qubo_HIPP(method[:1], gsize, hsize, hpairs, pen1=pen1,
                                                                      pen2=pen2, pen3=pen3, csize=csize)
    if OFFSET is None:
        continue
    ovar = hsize + total_hpairs
    ononz = hsize + total_hpairs * 3
    sys.stdout = open(folder_pre + config + "/" + f[:-3] + "out", 'w')
    # print("ExtraV:", disHpairs, ononz, len(vname2idx) - ovar, len(Q) - ononz)
    print(f[:-4], file=stderr)
    sampler = neal.SimulatedAnnealingSampler()
    # print("#Var-#NoneZ:", len(vname2idx), '-', len(Q))
    bqm = dimod.BinaryQuadraticModel.from_qubo(Q, OFFSET)
    if method[:1] == "D" or method[:1] == "E" or method[:1] == "F" or method[:1] == "L":
        for K in range(gsize):
            bqm.fix_variable(vname2idx[y_name(K + 1, 0)], 0)
            bqm.fix_variable(vname2idx[y_name(K + 1, len(hpairs[K]))], 1)
        print("#Var-#NoneZ:", bqm.num_variables, '-', bqm.num_interactions)
    start = time.time()
    response = sampler.sample(bqm, num_reads=num_reads, sweeps=sweeps)
    end = time.time()
    postT = response.info["timing"]['postprocessing_ns'] / 1000000.0
    preT = response.info["timing"]['preprocessing_ns'] / 1000000.0
    annealT = response.info["timing"]['sampling_ns'] / 1000000.0
    print("TimeRunning", end - start, postT, preT, annealT)
    best_set = range(hsize)
    for rid, r in enumerate(response.data()):
        objv = 0
        cur_hap_set = set()
        for i in range(hsize):
            if method[:1] == "E":
                if r.sample[vname2idx[x_name(i)]] == 0:
                    objv += 1
                    cur_hap_set.add(i)
            elif r.sample[vname2idx[x_name(i)]] == 1:
                objv += 1
                cur_hap_set.add(i)
        vnum, _ = check_constraint(constraint_set, r, vname2idx)
        for k in range(gsize):
            check_genotype(hpairs[k], cur_hap_set)
        if len(best_set) > len(cur_hap_set):
            best_set = list(cur_hap_set)
        print(rid, "\t", objv, "\t", r.energy, "\t", r.num_occurrences, "\t", vnum, "\t", len(cur_hap_set))
    print("BestSolution")
    for i in best_set:
        print(x_name(i))
    print("Objective:", len(best_set))
    sys.stdout.close()
