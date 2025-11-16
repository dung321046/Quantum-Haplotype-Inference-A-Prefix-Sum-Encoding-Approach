from ortools.linear_solver import pywraplp

from utils import *


def HIsolver(gsize, hsize, hpairs):
    """Linear programming sample."""
    # solver = pywraplp.Solver.CreateSolver('CBC')
    # Create the mip solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SAT')
    if not solver:
        print("Solver not found!")
        return
    X = []
    for i in range(hsize):
        X.append(solver.IntVar(0, 1, x_name(i)))
    P = dict()
    for i in range(gsize):
        for (u, v) in hpairs[i]:
            name = p_name(u, v)
            P[(u, v)] = solver.IntVar(0, 1, name)
            solver.Add(X[u] >= P[(u, v)])
            solver.Add(X[v] >= P[(u, v)])
            # print("(", u, v, ")", end=" ")
        # print("G", i)
        solver.Add(sum(P[(u, v)] for (u, v) in hpairs[i]) >= 1)
    print(solver.NumVariables(), hsize, sum([len(a) for a in hpairs]), file=sys.stderr)
    print('G-H-#V-#C', gsize, hsize, solver.NumVariables(), solver.NumConstraints())
    objective = solver.Objective()
    for i in range(hsize):
        objective.SetCoefficient(X[i], 1)
    objective.SetMinimization()
    # [END objective]

    # Solve the system.
    # [START solve]
    status = solver.Solve()
    # [END solve]

    # [START print_solution]
    if status == pywraplp.Solver.OPTIMAL:
        print('OptimalSolution', solver.Objective().Value())
        for i in range(hsize):
            if X[i].solution_value() == 1:
                print(x_name(i), end=",")
        print()
        # for k, v in P.items():
        #     print(p_name(k[0], k[1]), '=', v.solution_value())
    else:
        print('No Optimal')
    # [END print_solution]

    print('TimeRunning %f' % (solver.wall_time() / 1000.0))
    print('Problem solved in %d iterations' % solver.iterations())


if __name__ == "__main__":
    # files = ['Data2/Geno2Len4-preprocess.txt', 'DataM/00Gen6Hap128.txt', 'Data2/Geno5Len4-preprocess.txt']
    # for filename in files[1:2]:
    #     gsize, glen, g, hsize, h = read_data(filename)
    #     hpairs = create_pairs(h, gsize, hsize, glen, g)
    #     HIsolver(gsize, hsize, hpairs)
    import os

    folder_path = 'HapSet/DataPreprocess/'
    outputp = "output/"
    method = "SCIP"
    config = method
    if not os.path.exists(outputp):
        os.mkdir(outputp)
    if not os.path.exists(outputp + config):
        os.mkdir(outputp + config)
    for f in os.listdir(folder_path):
        filenameINP = os.path.join(folder_path, f)
        if not os.path.isfile(filenameINP):
            continue
        gsize, glen, g, hsize, h = read_data(filenameINP)
        hpairs = create_pairs(h, gsize, hsize, glen, g)
        sys.stdout = open(outputp + config + "/" + f[:-3] + "out", 'w')
        HIsolver(gsize, hsize, hpairs)
        sys.stdout.close()
