# GetBaseGraphs.py
from sage.all import *
import pdb
import numpy as np
from Classes import EGraph
import sys

configs = {}
if (sys.argv[0] == None):
    raise TypeError
else:
    configs["edges"] = int(sys.argv[1])
if (sys.argv[0] == None):
    raise TypeError
else:
    configs["verts"] = int(sys.argv[2])
if (sys.argv[0] == None):
    raise TypeError
else:
    configs["loops"] = int(sys.argv[3])
# --------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------Functions for Generating Base Graphs-------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------
# set the number of maximum number of multiple edges between two verts
# this is the number of times the same edge can be repeated in the base graph
max_repeat = int(float(2 * configs["edges"]) / configs["verts"])
if (2 * configs["edges"] % configs["verts"] > 1):
    max_repeat += 1


def reset_vector(vector, i):
    repeated = 1
    if (i == len(vector) - 1):
        return
    else:
        temp = vector[i]
        # set repeated
        if (i > 0):
            for k in range(1, min(i, max_repeat)):
                if (vector[i - k] == temp):
                    repeated += 1
                else:
                    break
        for index in range(i + 1, len(vector)):
            if (repeated == max_repeat):
                temp += 1
                repeated = 1
            else:
                repeated += 1
            vector[index] = temp


def Iterate_vector(vector, m):
    for i in range(0, len(vector)):
        if (vector[(len(vector) - 1) - i] != m - int(float(i) / max_repeat)):
            vector[(len(vector) - 1) - i] += 1
            reset_vector(vector, len(vector) - 1 - i)
            return True
    return False


def Construct_EGraph_no_loops(APE, IMA):
    G = DiGraph(configs["verts"], multiedges=true)
    IM = np.zeros((configs["verts"], configs["edges"]), dtype=np.int)
    for j in range(configs["edges"]):
        G.add_edge(APE[IMA[j]][0], APE[IMA[j]][1], j)
        IM[APE[IMA[j]][0]][j] = -1
        IM[APE[IMA[j]][1]][j] = 1
    E = EGraph(G, IM, configs)
    return E


def Construct_EGraph(LIMA, APE, IMA):
    G = DiGraph(configs["verts"], loops=true, multiedges=true)
    edge_count = 0
    IM = np.zeros((configs["verts"], configs["edges"]), dtype=np.int)
    for k in range(len(LIMA)):
        for j in range(LIMA[k]):
            G.add_edge(k, k, edge_count)
            edge_count += 1
    for j in range(len(IMA)):
        G.add_edge(APE[IMA[j]][0], APE[IMA[j]][1], edge_count)
        IM[APE[IMA[j]][0]][edge_count] = -1
        IM[APE[IMA[j]][1]][edge_count] = 1
        edge_count += 1
    E = EGraph(G, IM, configs)
    return E


# ----------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------Generate All Base Graphs---------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------
LA = Partitions(configs["loops"], max_length=configs["loops"]).list()
All_Edges = []
for i in range(configs["verts"]):
    for j in range(i + 1, configs["verts"]):
        A = [0] * 2
        A[0] = i
        A[1] = j
        All_Edges.append(A)

max_value = len(All_Edges) - 1
All_EGraphs = []
All_Graphs = []

Output_file = str(configs["verts"]) + "V" + str(configs["edges"]) + "E" + str(configs["loops"]) + "LBaseGraphs.txt"
with open(Output_file, "w") as f:
    Incidence_Matrix_Array = [0] * (configs["edges"] - configs["loops"])
    reset_vector(Incidence_Matrix_Array, 0)
    if (configs["loops"] > 0):
        for la in LA:
            while (True):
                E = Construct_EGraph(la, All_Edges, Incidence_Matrix_Array)
                if (E.G.is_connected()):
                    UG = E.G.to_undirected()
                    is_new = True
                    for g in All_Graphs:
                        if (UG.is_isomorphic(g)):
                            is_new = False
                            break
                    if (is_new):
                        All_Graphs.append(UG)
                        All_EGraphs.append(E)
                        f.write(str(la) + "\n")
                        f.write(str(Incidence_Matrix_Array) + "\n\n")
                if (not Iterate_vector(Incidence_Matrix_Array, max_value)):
                    break
            Incidence_Matrix_Array = [0] * (configs["edges"] - configs["loops"])
            reset_vector(Incidence_Matrix_Array, 0)
    # No Loops
    else:
        la = ()
        while (True):
            E = Construct_EGraph_no_loops(All_Edges, Incidence_Matrix_Array)
            if (E.G.is_connected()):
                UG = E.G.to_undirected()
                is_new = True
                for g in All_Graphs:
                    if (UG.is_isomorphic(g)):
                        is_new = False
                        break
                if (is_new):
                    All_Graphs.append(UG)
                    All_EGraphs.append(E)
                    f.write(str(la) + "\n")
                    f.write(str(Incidence_Matrix_Array) + "\n\n")

            if (not Iterate_vector(Incidence_Matrix_Array, max_value)):
                break

print(len(All_EGraphs))
