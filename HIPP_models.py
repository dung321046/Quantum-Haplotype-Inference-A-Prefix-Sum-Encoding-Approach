import formulationPrefixHI
import formulationQHI
import formulationSCP
from sys import stderr


def load_qubo_HIPP(method, gsize, hsize, hpairs, pen1, pen2, pen3=None, csize=2):
    if method == "B":
        OFFSET, Q, vname2idx, constraint_set = formulationSCP.formulate(gsize, hsize, hpairs)
    elif method == "D":
        OFFSET, Q, vname2idx, constraint_set = formulationPrefixHI.formulate(gsize, hsize, hpairs, pen1, pen2)

    elif method == "H":
        OFFSET, Q, vname2idx, constraint_set = formulationQHI.formulate(csize, gsize, hsize, hpairs)
    else:
        print("ERRORR", file=stderr)
        return None, None, None, None
    return OFFSET, Q, vname2idx, constraint_set
