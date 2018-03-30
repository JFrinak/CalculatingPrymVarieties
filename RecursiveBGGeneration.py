# RecursiveBGGeneration.py
from sage.all import *
import pdb
import numpy as np
from Classes import EGraph, GraphShell
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


# ----------------------------------------------------------------------------------------------------------
def LoadGraphs(verts, edges, loops):
    # type: (int, int, int) -> [GraphShell]
    input_file = str(verts) + "V" + str(edges) + "E" + str(loops) + "LBaseGraphs.txt"
    All_GraphShells = []
    with open(input_file, "r") as inF:
        row = 0
        LIMA = []
        IMA = []
        for line in inF:
            if (line == "\n"):
                gs = GraphShell(LIMA, IMA, configs)
                All_GraphShells.append(gs)
                row = 0
                LIMA = []
                IMA = []
            elif (row == 0):
                if (loops > 0):
                    line = line.rstrip()
                    num = ""
                    for ch in line:
                        if (ch not in ['(', ',', ')', '[', ']']):
                            num += ch
                        elif (ch in [',', ')', ']']):
                            if num != "":
                                integer = int(num)
                                LIMA.append(integer)
                                num = ""
                row += 1
            elif (row == 1):
                IMA = []
                line = line.rstrip()
                num = ""
                for ch in line:
                    if (ch not in ['[', ',', ']']):
                        num += ch
                    elif (ch in [',', ']']):
                        if num != "":
                            integer = int(num)
                            IMA.append(integer)
                            num = ""

    return All_GraphShells


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
def swapVertices(APE,IMA,vert1,vert2):
    #Adjust IMA
    newIMA = []
    for e in range(configs["edges"]-configs["loops"]):
        edge = [x for x in APE[IMA[e]]]
        newEdge=[0]*2
        if bool(vert1 in edge) != bool(vert2 in edge):
            if(edge[0] == vert1 and edge[1]< vert2):
                newEdge = [edge[1],vert2]
            elif(edge[0] == vert1 and vert2<edge[1]):
                newEdge = [vert2,edge[1]]
            elif(edge[0] == vert2):
                newEdge = [vert1,edge[1]]
            elif(edge[1] == vert1):
                newEdge = [edge[0], vert2]
            elif(edge[1] == vert2 and edge[0]<vert1):
                newEdge = [edge[0],vert1]
            elif(edge[1]==vert2 and vert1<edge[0]):
                newEdge = [vert1,edge[0]]
            newEdgeIndex = APE.index(newEdge)

        else:
            newEdgeIndex = IMA[e]

        if(len(newIMA)==0):
            newIMA.append(newEdgeIndex)
        else:
            for i in range(len(newIMA)+1):
                if(i == len(newIMA)):
                    newIMA.append(newEdgeIndex)
                else:
                    if(newEdgeIndex<= newIMA[i]):
                        newIMA.insert(i,newEdgeIndex)
                        break
    return newIMA


# ----------------------------------------------------------------------------------------------------------
All_Edges = []
for i in range(configs["verts"]):
    for j in range(i + 1, configs["verts"]):
        A = [0] * 2
        A[0] = i
        A[1] = j
        All_Edges.append(A)

max_value = len(All_Edges) - 1
IMA = [0] * (configs["edges"]-configs["loops"])
la = ()
All_Graphs = []
All_EGraphs = []

Output_file = str(configs["verts"]) + "V" + str(configs["edges"]) + "E" + str(configs["loops"]) + "LBaseGraphs.txt"
with open(Output_file, "w") as f:
    count = 0
    if configs["edges"]<configs["verts"]:
        lowerDimensionAllEdges = []
        for i in range(configs["verts"]-1):
            for j in range(i+1,configs["verts"]-1):
                A = [0]*2
                A[0] = i
                A[1] = j
                lowerDimensionAllEdges.append(A)

        lowerDimensionGraphShells = LoadGraphs(configs["verts"]-1,configs["edges"]-1,0);
        for shell in lowerDimensionGraphShells:
            oldIMA = [x for x in shell.IMA]
            #Need to adjust oldIMA to the context where we have the right number of vertices.
            newIMA = [x + lowerDimensionAllEdges[x][0] for x in oldIMA]
            possibleNewEdges = []
            for x in range(1, configs["verts"]):
                possibleNewEdges.append(x * configs["verts"] - x * (x + 1) / 2 - 1)
            for k in range(len(possibleNewEdges)):
                IMA = [x for x in newIMA]
                IMA.append(possibleNewEdges[k])
                IMA.sort()
                E = Construct_EGraph_no_loops(All_Edges, IMA)
                UG = E.G.to_undirected()
                if E.G.is_connected():
                    is_new = True
                    for g in All_Graphs:
                        if (UG.is_isomorphic(g)):
                            is_new = False
                            break
                    if (is_new):
                        All_Graphs.append(UG)
                        All_EGraphs.append(E)
                        f.write(str(la) + "\n")
                        f.write(str(IMA) + "\n\n")
            count += 1
            num = float(count) / len(lowerDimensionGraphShells)
            percent = int(num * 100)
            string = str(percent) + "%"
            print(string)
        f.write(str(len(All_EGraphs)) + " Base Graphs")

    elif configs["loops"] == 0:
        lowerDimensionGraphShells = LoadGraphs(configs["verts"], configs["edges"] - 1, configs["loops"])
        for shell in lowerDimensionGraphShells:
            for edge in range(len(All_Edges)):
                IMA = [x for x in shell.IMA]
                IMA.append(edge)
                IMA.sort()
                E = Construct_EGraph_no_loops(All_Edges, IMA)
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
                    f.write(str(IMA) + "\n\n")
            count += 1
            num = float(count) / len(lowerDimensionGraphShells)
            percent = int(num * 100)
            string = str(percent) + "%"
            print(string)
        f.write(str(len(All_EGraphs)) + " Base Graphs")

    else:
        lowerDimensionGraphShells = LoadGraphs(configs["verts"], configs["edges"] - 1, configs["loops"] - 1)
        for shell in lowerDimensionGraphShells:
            IMA = [x for x in shell.IMA]
            loopList = [x for x in shell.LIMA]

            for v in range(configs["verts"]):
                needIMAChange = v>0 and (len(loopList)==0 or (v<len(loopList) and loopList[v]==loopList[v-1])
                                or v>len(loopList))
                if needIMAChange:
                    firstIndex = 0
                    if len(loopList)==0:
                        loopList.append(1)
                        firstindex = 0
                    elif v<len(loopList) and loopList[v] == loopList[v-1]:
                        firstIndex = loopList.index(loopList[v]) #firstIndex might not be v-1, i.e. lL = (1,1,1) and v=2
                        loopList[firstIndex]+=1
                    elif v>len(loopList):
                        firstIndex = len(loopList)
                        loopList.append(1)
                    else:
                        raise ValueError('We have not handled all IMA change cases')
                    newIMA = swapVertices(All_Edges,IMA,firstIndex,v)
                    E = Construct_EGraph(loopList, All_Edges, newIMA)
                    UG = E.G.to_undirected()
                    is_new = True
                    for g in All_Graphs:
                        if (UG.is_isomorphic(g)):
                            is_new = False
                            break
                    if (is_new):
                        All_Graphs.append(UG)
                        All_EGraphs.append(E)
                        f.write(str(tuple(loopList)) + "\n")
                        f.write(str(newIMA) + "\n\n")
                    loopList = [x for x in shell.LIMA]
                    IMA = [x for x in shell.IMA]

                else:
                    if v==0 and len(loopList)==0:
                        loopList.append(1)
                    elif v == 0:
                        loopList[v] += 1
                    elif v<len(loopList) and loopList[v]<loopList[v-1]:
                        loopList[v] +=1
                    elif v == len(loopList):
                        loopList.append(1)
                    else:
                        raise ValueError('We have not handled all non-IMA change cases')
                    E = Construct_EGraph(loopList, All_Edges, IMA)
                    UG = E.G.to_undirected()
                    is_new = True
                    for g in All_Graphs:
                        if (UG.is_isomorphic(g)):
                            is_new = False
                            break
                    if (is_new):
                        All_Graphs.append(UG)
                        All_EGraphs.append(E)
                        f.write(str(tuple(loopList)) + "\n")
                        f.write(str(IMA) + "\n\n")
                    loopList = [x for x in shell.LIMA]
                    IMA = [x for x in shell.IMA]

            count += 1
            num = float(count) / len(lowerDimensionGraphShells)
            percent = int(num * 100)
            string = str(percent) + "%"
            print(string)
        f.write(str(len(All_EGraphs)) + " Base Graphs")
