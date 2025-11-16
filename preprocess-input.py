import os

import numpy as np

from utils import *

folder_path = 'DataT50/'
out_folder = "DataF50/"
if not os.path.exists(out_folder):
    os.mkdir(out_folder)
for f in os.listdir(folder_path):
    filenameINP = os.path.join(folder_path, f)
    if os.path.isfile(filenameINP):
        gsize, glen, g, hsize, h = read_data(filenameINP)
        hpairs = []
        for i in range(gsize):
            hpairs.append([])
        deg = np.zeros(hsize)
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
                    deg[i] += 1
                    deg[j] += 1
        imp_idx = set()
        for k in range(gsize):
            npair = 0
            for (u, v) in hpairs[k]:
                if deg[u] == 1 and deg[v] == 1:
                    continue
                npair += 1
                imp_idx.add(u)
                imp_idx.add(v)
            if npair == 0:
                imp_idx.add(hpairs[k][0][0])
                imp_idx.add(hpairs[k][0][1])
        filenameOUT = out_folder + f[:2] + "Gen" + str(gsize) + "Hap" + str(len(imp_idx)) + ".txt"
        sys.stdout = open(filenameOUT, 'w')
        print(gsize, glen)
        for gi in g:
            print(gi)
        print()
        print()
        print(len(imp_idx))
        for i in imp_idx:
            print(h[i])
        sys.stdout.close()
