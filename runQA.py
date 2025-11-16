import time

import dimod
from dwave.embedding.chain_strength import uniform_torque_compensation
from dwave.system.samplers import DWaveSampler

import HIPP_models
from utils import *

folder_path = 'DataT50'
out_path = "output/"
pen1 = 2
pen2 = 1
pen3 = 1
DUMP_EMB = False
pemb = "NCC50/B/"
method = "B-Q_NCC50" + "-"

# method = "B-Q-P" + str(pen1) + "_" + str(pen2) + "_" + pemb + "-"

# method = "F-Q-P" + str(pen1) + "_" + str(pen2) + "_" + str(pen3) + "-"

csize = 2
# method = "H" + str(csize) + "-Q_NCC50-"
num_reads = 1000
chain_strength = 0.20
annealing_time = 200
# chain_strength =  "default"
schedules = {1: [(0.0, 0.0), (200.0, 1.0)],
             "Pause1": [(0.0, 0.0), (50.0, 0.5), (100.0, 0.5), (200.0, 1.0)],
             "Pause2": [(0.0, 0.0), (75.0, 0.5), (125.0, 0.5), (200.0, 1.0)],
             "Pause3": [(0.0, 0.0), (100.0, 0.5), (150.0, 0.5), (200.0, 1.0)],
             "Pause4": [(0.0, 0.0), (80.0, 0.5), (280.0, 0.5), (400.0, 1.0)],
             "Quench1": [(0.0, 0.0), (180.0, 0.5), (200.0, 1.0)],
             "Quench2": [(0.0, 0.0), (190.0, 0.5), (200.0, 1.0)],
             "Quench3": [(0.0, 0.0), (195.0, 0.5), (200.0, 1.0)],
             "Quench4": [(0.0, 0.0), (160.0, 0.5), (200.0, 1.0)],
             "QP1_20": [(0.0, 0.0), (160.0, 0.5), (180.0, 0.5), (200.0, 1.0)],
             "QP2_20": [(0.0, 0.0), (170.0, 0.5), (190.0, 0.5), (200.0, 1.0)],
             "QP3_20": [(0.0, 0.0), (175.0, 0.5), (195.0, 0.5), (200.0, 1.0)],
             "QP4_20": [(0.0, 0.0), (140.0, 0.5), (160.0, 0.5), (200.0, 1.0)],
             "QP1_40": [(0.0, 0.0), (140.0, 0.5), (180.0, 0.5), (200.0, 1.0)],
             "QP2_40": [(0.0, 0.0), (150.0, 0.5), (190.0, 0.5), (200.0, 1.0)],
             "QP3_40": [(0.0, 0.0), (155.0, 0.5), (195.0, 0.5), (200.0, 1.0)],
             "QP4_40": [(0.0, 0.0), (120.0, 0.5), (160.0, 0.5), (200.0, 1.0)],
             "BangBang1": [(0.0, 0.0), (50.0, 0.4), (150.0, 0.6), (200.0, 1.0)],
             "BangBang2": [(0.0, 0.0), (50.0, 0.3), (150.0, 0.4), (200.0, 1.0)],
             "BangBang3": [(0.0, 0.0), (40.0, 0.3), (160.0, 0.4), (200.0, 1.0)],
             # "Quench3": [(0.0, 0.0), (100.0, 0.5), (150.0, 0.5), (200.0, 1.0)],
             21: [(0.0, 0.0), (50.0, 0.5), (150.0, 0.5), (200.0, 1.0)],
             22: [(0.0, 0.0), (75.0, 0.5), (125.0, 0.5), (200.0, 1.0)],
             23: [(0.0, 0.0), (75.0, 0.6), (125.0, 0.6), (200.0, 1.0)],
             3: [(0.0, 0.0), (100.0, 0.2), (150.0, 0.2), (200.0, 1.0)],
             4: [(0.0, 0.0), (130.0, 0.4), (180.0, 0.4), (200.0, 1.0)],
             12: [(0.0, 0.0), (40.0, 0.4), (180.0, 0.4), (200.0, 1.0)],
             11: [(0.0, 0.0), (40.0, 0.5), (120.0, 0.5), (200.0, 1.0)],
             13: [(0.0, 0.0), (40.0, 0.5), (130.0, 0.5), (200.0, 1.0)],
             14: [(0.0, 0.0), (30.0, 0.5), (160.0, 0.5), (200.0, 1.0)],
             15: [(0.0, 0.0), (40.0, 0.4), (160.0, 0.6), (200.0, 1.0)],
             400: [(0.0, 0.0), (400.0, 1.0)]}
sid = 1
schedule = schedules[sid]
config = method + str(num_reads) + "-" + str(chain_strength) + "s" + str(sid) + "_A" + str(annealing_time)
if not os.path.exists(out_path + config):
    os.mkdir(out_path + config)
TT = 0

for id_, f in enumerate(os.listdir(folder_path)):
    TT += 1
    sys.stdout = sys.__stdout__
    filenameINP = os.path.join(folder_path, f)
    if not os.path.isfile(filenameINP):
        continue
    gsize, glen, g, hsize, h = read_data(filenameINP)
    hpairs = create_pairs(h, gsize, hsize, glen, g)
    OFFSET, Q, vname2idx, constraint_set = HIPP_models.load_qubo_HIPP(method[:1], gsize, hsize, hpairs, pen1=pen1,
                                                                      pen2=pen2, pen3=pen3, csize=csize)
    if OFFSET is None:
        continue
    import json

    sys.stdout = open(out_path + config + "/" + f[:-3] + "out", 'w')
    token = "put-your-token-here"
    qpu_solver = DWaveSampler(token=token)
    sol_specs = qpu_solver.properties
    # print("--", sol_specs)
    print("#Var-#NoneZ:", len(vname2idx), '-', len(Q))
    start = time.time()
    from dwave.system.composites import FixedEmbeddingComposite

    with open(out_path + pemb + f, "r") as fff:
        embbed = json.load(fff)
    embbed = {int(k): v for k, v in embbed.items()}
    try:
        sampler = FixedEmbeddingComposite(qpu_solver, embedding=embbed)
    except Exception as except_error:
        from dwave.system.composites import EmbeddingComposite

        sampler = EmbeddingComposite(qpu_solver)
    bqm = dimod.BinaryQuadraticModel.from_qubo(Q, OFFSET)
    if method[:1] == "D" or method[:1] == "E" or method[:1] == "F":
        for K in range(gsize):
            bqm.fix_variable(vname2idx[y_name(K + 1, 0)], 0)
            bqm.fix_variable(vname2idx[y_name(K + 1, len(hpairs[K]))], 1)

    if chain_strength == "default":
        response = sampler.sample(bqm, num_reads=num_reads, annealing_time=annealing_time, return_embedding=True)
    else:
        try:
            c_str = uniform_torque_compensation(bqm, prefactor=chain_strength)
            # c_str = chain_strength
            print(c_str, end=" ")
            print(schedule)
            response = sampler.sample(bqm, chain_strength=c_str, num_reads=num_reads, return_embedding=True,
                                      anneal_schedule=schedule)
        except Exception as error:
            from dwave.system.composites import EmbeddingComposite

            sampler = EmbeddingComposite(qpu_solver)
            c_str = uniform_torque_compensation(bqm, prefactor=chain_strength)
            response = sampler.sample(bqm, chain_strength=c_str, num_reads=num_reads, return_embedding=True,
                                      anneal_schedule=schedule)
    chains = response.info["embedding_context"]["embedding"].values()
    if DUMP_EMB:
        with open(out_path + f, "w") as fff:
            json.dump(response.info["embedding_context"]["embedding"], fff)
    chain_length = max(len(chain) for chain in chains)
    end = time.time()
    print("TimeRunning", end - start)
    best_set = range(hsize)
    for rid, r in enumerate(response.data()):
        objv = 0
        cur_hap_set = set()
        for i in range(hsize):
            if r.sample[vname2idx[x_name(i)]] == 1:
                objv += 1
                cur_hap_set.add(i)
        vnum, _ = check_constraint(constraint_set, r, vname2idx)
        for k in range(gsize):
            check_genotype(hpairs[k], cur_hap_set)
        if len(best_set) > len(cur_hap_set):
            best_set = list(cur_hap_set)
        print(rid, "\t", objv, "\t", r.energy, "\t", vnum, "\t", len(cur_hap_set), "\t", r.num_occurrences,
              "\t", r.chain_break_fraction)
    print("BestSolution")
    for i in best_set:
        print(x_name(i))
    print("Objective:", len(best_set))
    if chain_strength == "default":
        print("ChainS-ChainL", response.info["embedding_context"]["chain_strength"], chain_length)
    else:
        print("ChainS-ChainL", chain_strength, chain_length)
    print("Info", response.info["timing"])
    print("Emb", response.info["embedding_context"])
    emb = response.info["embedding_context"]["embedding"]
    for name, id in vname2idx.items():
        if id in emb:
            print(name, len(emb[id]))
    sys.stdout.close()
