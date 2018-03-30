# GetCovers.py
from sage.all import *
import pdb
import numpy as np
from Classes import *
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
# ----------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------Function for Loading Base Graphs from Text file----------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------
All_Edges = []

for i in range(configs["verts"]):
	for j in range(i + 1, configs["verts"]):
		A = [0] * 2
		A[0] = i
		A[1] = j
		All_Edges.append(A)


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


def LoadGraphs():
	input_file = str(configs["verts"]) + "V" + str(configs["edges"]) + "E" + str(configs["loops"]) + "LBaseGraphs.txt"
	All_EGraphs = []
	with open(input_file, "r") as inF:
		row = 0
		for line in inF:
			if (line == "\n"):
				if (configs["loops"] > 0):
					E = Construct_EGraph(LIMA, All_Edges, IMA)
				else:
					E = Construct_EGraph_no_loops(All_Edges, IMA)
				if 1 in E.G.degree():
					All_EGraphs.append(E)
				row = 0
			elif (row == 0):
				if configs["loops"] > 0:
					LIMA = []
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
						integer = int(num)
						IMA.append(integer)
						num = ""
	return All_EGraphs


# ----------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Functions for Generating Cover Graphs---------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------
def reset_Boolean_List(B, col):
	if (col == len(B) - 1):
		return
	for i in range(col + 1, len(B)):
		B[i] = False


def iterate_Boolean_List(B):
	for i in range(len(B)):
		if (B[len(B) - 1 - i] == False):
			B[len(B) - 1 - i] = True
			reset_Boolean_List(B, len(B) - 1 - i)
			return True
	return False


def get_Possible_Cover_IM_no_CCEdges(E, CV, CE):
	Possible_Cover_IM = np.zeros((2 * configs["verts"], 2 * configs["edges"]), dtype=np.int)
	# copy base graph
	for v in range(configs["verts"]):
		for e in range(configs["edges"]):
			Possible_Cover_IM[2 * v][2 * e] = E.IM[v][e]

	# Covering Information
	for e in range(configs["edges"]):
		# Edge is covered
		if (CE[e]):
			if (CV[E.start(e)]):
				Possible_Cover_IM[2 * E.start(e) + 1][2 * e + 1] = -1
			else:
				Possible_Cover_IM[2 * E.start(e)][2 * e + 1] = -1
			if (CV[E.end(e)]):
				Possible_Cover_IM[2 * E.end(e) + 1][2 * e + 1] = 1
			else:
				Possible_Cover_IM[2 * E.end(e)][2 * e + 1] = 1

	# No loops in this base graph
	return Possible_Cover_IM


def get_Possible_Cover_IM(E, CV, CE, cc_edges, cc_edges_iscrossed, cc_edges_locations):
	Possible_Cover_IM = np.zeros((2 * configs["verts"], 2 * configs["edges"]), dtype=np.int)
	# copy base graph
	for v in range(configs["verts"]):
		for e in range(configs["edges"]):
			Possible_Cover_IM[2 * v][2 * e] = E.IM[v][e]
	# e is a loop and must be crossed (no loops in cover graphs)
	loop_num = 0
	graph_loops = E.G.loop_edges()
	for loop in graph_loops:
		vertex = loop[0]
		Possible_Cover_IM[2 * vertex][2 * loop_num] = -1
		Possible_Cover_IM[2 * vertex + 1][2 * loop_num + 1] = -1
		Possible_Cover_IM[2 * vertex][2 * loop_num + 1] = 1
		Possible_Cover_IM[2 * vertex + 1][2 * loop_num] = 1
		loop_num += 1
	# Covering Information
	for e in range(configs["loops"], configs["edges"]):
		edge_start = E.start(e)
		edge_end = E.end(e)
		# Edge is part of ccl and should be crossed
		if (cc_edges[e] and cc_edges_iscrossed[cc_edges_locations.index(e)]):
			Possible_Cover_IM[2 * edge_start][2 * e] = -1
			Possible_Cover_IM[2 * edge_start + 1][2 * e + 1] = -1
			Possible_Cover_IM[2 * edge_end][2 * e + 1] = 1
			Possible_Cover_IM[2 * edge_end + 1][2 * e] = 1
			Possible_Cover_IM[2 * edge_end][2 * e] = 0

		elif (CE[e]):
			if (CV[edge_start]):
				Possible_Cover_IM[2 * edge_start + 1][2 * e + 1] = -1
			else:
				Possible_Cover_IM[2 * edge_start][2 * e + 1] = -1
			if (CV[edge_end]):
				Possible_Cover_IM[2 * edge_end + 1][2 * e + 1] = 1
			else:
				Possible_Cover_IM[2 * edge_end][2 * e + 1] = 1
	return Possible_Cover_IM


def edge_start(M, e):
	for v in range(2 * configs["verts"]):
		if (M[v][e] == -1):
			return v
	return -1


def edge_end(M, e):
	for v in range(2 * configs["verts"]):
		if (M[v][e] == 1):
			return v
	return -1


def is_Cover_no_CCEdges(PCIM, CVerts, CEdges):
	for e in range(configs["edges"]):
		start1 = -1
		start2 = -1
		end1 = -1
		end2 = -1
		if (e >= configs["loops"] and CEdges[e]):
			start1 = edge_start(PCIM, 2 * e + 1)
			end1 = edge_end(PCIM, 2 * e + 1)
			if (CVerts[edge_start(PCIM, 2 * e) / 2]):
				start2 = edge_start(PCIM, 2 * e) + 1
			else:
				start2 = edge_start(PCIM, 2 * e)
			if (CVerts[edge_end(PCIM, 2 * e) / 2]):
				end2 = edge_end(PCIM, 2 * e) + 1
			else:
				end2 = edge_end(PCIM, 2 * e)
			if (start1 != start2 or end1 != end2):
				return false
		elif (e >= configs["loops"]):
			start1 = edge_start(PCIM, 2 * e)
			end1 = edge_end(PCIM, 2 * e)
			if (CVerts[edge_start(PCIM, 2 * e) / 2]):
				start2 = edge_start(PCIM, 2 * e) + 1
			else:
				start2 = edge_start(PCIM, 2 * e)
			if (CVerts[edge_end(PCIM, 2 * e) / 2]):
				end2 = edge_end(PCIM, 2 * e) + 1
			else:
				end2 = edge_end(PCIM, 2 * e)
			if (start1 != start2 or end1 != end2):
				return false
	# We have constructed the loops such that the will always satisfy the cover conditions.
	return true


def is_Cover(PCIM, CVerts, CEdges, cc_edge_locations, cc_edges_isCrossed):
	# create isCrossed
	isCrossed = [false] * configs["edges"]
	for e in range(configs["edges"]):
		if (e in cc_edge_locations):
			if (cc_edges_isCrossed[cc_edge_locations.index(e)]):
				isCrossed[e] = True
	for e in range(configs["edges"]):
		start1 = -1
		start2 = -1
		end1 = -1
		end2 = -1
		if (not isCrossed[e] and e >= configs["loops"]):
			if (CEdges[e]):
				start1 = edge_start(PCIM, 2 * e + 1)
				end1 = edge_end(PCIM, 2 * e + 1)
				if (CVerts[edge_start(PCIM, 2 * e) / 2]):
					start2 = edge_start(PCIM, 2 * e) + 1
				else:
					start2 = edge_start(PCIM, 2 * e)
				if (CVerts[edge_end(PCIM, 2 * e) / 2]):
					end2 = edge_end(PCIM, 2 * e) + 1
				else:
					end2 = edge_end(PCIM, 2 * e)
				if (start1 != start2 or end1 != end2):
					return false
			else:
				if (CVerts[edge_start(PCIM, 2 * e) / 2]):
					return false
				if (CVerts[edge_end(PCIM, 2 * e) / 2]):
					return false
	# We have constructed the loops such that the will always satisfy the cover conditions.
	# Crossed edges will also satisfy the cover criteria by construction
	return true

# def vertCoverIsAdmissible(E,CVerts):
# 	for v in range(configs["verts"]):
# 		if E.G.degree()[v]==2 and CVerts[v]==True:
# 			return False
# 	return True

# ----------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------Generate All Cover Graphs--------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------

# Load all base Graphs
All_EGraphs = LoadGraphs()
print("Done Loading Base Graphs...")
print("Total number of base graphs: " + str(len(All_EGraphs)))

APEC = []
CEdges = [false] * configs["edges"]
for l in range(configs["loops"]):
	CEdges[l] = True
while (True):
	CE = list(CEdges)
	APEC.append(CE)
	if (not iterate_Boolean_List(CEdges)):
		break
APVC = []
CVerts = [False] * configs["verts"]
#If there are loops vertex 0 will always need to be covered
if(configs["loops"]>0):
	CVerts[0]=True
while (True):
	CV = list(CVerts)
	APVC.append(CV)
	if (not iterate_Boolean_List(CVerts)):
		break

totalNumberCGs = 0
output_file = str(configs["verts"]) + "V" + str(configs["edges"]) + "E" + str(configs["loops"]) + "L_FS_output.txt"
second_output = str(configs["verts"]) + "V" + str(configs["edges"]) + "E" + str(configs["loops"]) + "L_FS_Graphs.txt"

with open(output_file, "w") as f:
	with open(second_output, "w") as g:
		count = 0
		for E in All_EGraphs:
			All_CGS_For_G = []
			for vc in APVC:
				for ec in APEC:
					E.setNumberCCEdges(vc, ec)
					# There cannot be any loops in cover graph therefore all loop vertices must be covered.
					# Recall all loops are already covered by construction of APEC
					if (E.are_LVs_covered):
						if (E.numberCCEdges > 0):
							# Determine which edges are completely covered
							cc_edge_locations = []
							# Get locations of completely covered edges
							for e in range(configs["edges"]):
								if (E.ccEdges[e]):
									cc_edge_locations.append(e)
							cc_edges_iscrossed = [False] * E.numberCCEdges
							# make sure loops are always crossed
							for i in range(len(cc_edge_locations)):
								if (cc_edge_locations[i] < configs["loops"]):
									#loops come first in cc_edge_locations
									cc_edges_iscrossed[i] = True
							while (True):
								numberCrossedCCEdges = cc_edges_iscrossed.count(True)
								crossingThreshold = E.numberCCEdges+configs["loops"]-2*numberCrossedCCEdges
								if(crossingThreshold>=0):
									Possible_IM = get_Possible_Cover_IM(E, vc, ec, E.ccEdges, cc_edges_iscrossed,
																		cc_edge_locations)
									if (is_Cover(Possible_IM, vc, ec, cc_edge_locations, cc_edges_iscrossed)):
										CG = CoverGraph(Possible_IM, ec, vc, configs)
										if (CG.G.is_connected() and CG.FS and CG.Basis_Linear_Forms.nrows() > 0):
											is_new = True
											for cg in All_CGS_For_G:
												if (CG.G.is_isomorphic(cg.G)):
													is_new = False
													break
											if (is_new):
												All_CGS_For_G.append(CG)
												f.write(str(CG.Basis_Linear_Forms))
												f.write("\n\n")
												g.write(CG.output())
												totalNumberCGs+=1

								if (not iterate_Boolean_List(cc_edges_iscrossed)):
									break

						# No CCEdges
						else:
							Possible_IM = get_Possible_Cover_IM_no_CCEdges(E, vc, ec)
							if (is_Cover_no_CCEdges(Possible_IM, vc, ec)):
								CG = CoverGraph(Possible_IM, ec, vc, configs)
								if (CG.G.is_connected() and CG.FS and CG.Basis_Linear_Forms.nrows() > 0):
									is_new = True
									for cg in All_CGS_For_G:
										if (CG.G.is_isomorphic(cg.G)):
											is_new = False
											break
									if (is_new):
										All_CGS_For_G.append(CG)
										f.write(str(CG.Basis_Linear_Forms))
										f.write("\n\n")
										g.write(CG.output())
										totalNumberCGs+=1
			count += 1
			num = float(count) / len(All_EGraphs)
			percent = int(num * 100)
			string = str(percent) + "%"
			print(string)

		if(totalNumberCGs>0):
			g.write("\n\n")
			g.write("Total Number of Cover Graphs: " + str(totalNumberCGs))
			f.write("\n\n")
			f.write("Total Number of Cover Graphs: " + str(totalNumberCGs))
print(str(totalNumberCGs) + " Cover Graphs\n")
