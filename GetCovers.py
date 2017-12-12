#GetCovers.py
from sage.all import *
import pdb
import numpy as np
from Classes import *
import sys

configs={}
if(sys.argv[0] == None):
	raise TypeError
else:
	configs["edges"] = int(sys.argv[1])
if(sys.argv[0] == None):
	raise TypeError
else:
	configs["verts"] = int(sys.argv[2])
if(sys.argv[0] == None):
	raise TypeError
else:
	configs["loops"] = int(sys.argv[3])
#----------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------Function for Loading Base Graphs from Text file----------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------
All_Edges = []

for i in range(configs["verts"]):
    for j in range(i+1,configs["verts"]):
        A = [0]*2
        A[0] = i
        A[1] = j
        All_Edges.append(A)

def Construct_EGraph_no_loops(APE,IMA):
	G = DiGraph(configs["verts"],multiedges=true)
	IM = np.zeros((configs["verts"],configs["edges"]),dtype=np.int)
    	for j in range(configs["edges"]):
        	G.add_edge(APE[IMA[j]][0],APE[IMA[j]][1],j) 
		IM[APE[IMA[j]][0]][j] = -1
		IM[APE[IMA[j]][1]][j] = 1
    	E = EGraph(G,IM,configs)
    	return E
    
def Construct_EGraph(LIMA,APE,IMA):
	G = DiGraph(configs["verts"],loops=true, multiedges=true)
	edge_count = 0
	IM = np.zeros((configs["verts"],configs["edges"]),dtype=np.int)
    	for k in range(len(LIMA)):
        	for j in range(LIMA[k]):
           		G.add_edge(k,k,edge_count)
           		edge_count+=1
    	for j in range(len(IMA)):
        	G.add_edge(APE[IMA[j]][0],APE[IMA[j]][1],edge_count)
		IM[APE[IMA[j]][0]][edge_count] = -1
		IM[APE[IMA[j]][1]][edge_count] = 1
		edge_count+=1 
    	E = EGraph(G,IM,configs)
    	return E
def LoadGraphs():
	input_file = str(configs["verts"])+"V"+str(configs["edges"])+"E"+str(configs["loops"])+"LBaseGraphs.txt"
	All_EGraphs = []
	with open(input_file,"r") as inF:
		row = 0
		for line in inF:
			if(line == "\n"):
				if(configs["loops"]>0):
            				E = Construct_EGraph(LIMA,All_Edges,IMA)
				else:
        				E = Construct_EGraph_no_loops(All_Edges,IMA)
				All_EGraphs.append(E)
				row = 0
			elif(row == 0):
				if(configs["loops"]>0):
					LIMA = []
					line=line.rstrip()
					num =""
					for ch in line:
						if(ch not in ['(',',',')','[',']']):
							num += ch
						elif(ch in [',',')',']']):
							integer = int(num)
							LIMA.append(integer)
							num = ""
				row+=1
			elif(row == 1):
				IMA = []
				line = line.rstrip()
				num =""
				for ch in line:
					if(ch not in ['[',',',']']):
						num += ch
					elif(ch in [',',']']):
						integer = int(num)
						IMA.append(integer)
						num = ""
	return All_EGraphs
#----------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------Functions for Generating Cover Graphs---------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------
def is_true_after(CE,col):
    for i in range(col+1,len(CE)):
        if(not CE[i]):
            return False
    return True
    
def Iterate_Boolean_List(B,col,size):
    if(B[col]):
        if(col==size-1):
            return
        else:
            Iterate_Boolean_List(B,col+1,size)
    else:
        if(col==size-1):
            B[col]=True
            return
        else:
            if(is_true_after(B,col)):
                B[col]=True
                for i in range(col+1,size):
                    B[i]=False
                return
            else:
                Iterate_Boolean_List(B,col+1,size)

                
def get_Possible_Cover_IM_no_loops(E,CV,CE):
    Possible_Cover_IM = np.zeros((2*configs["verts"],2*configs["edges"]), dtype=np.int)
    #copy base graph
    for v in range(configs["verts"]):
        for e in range(configs["edges"]):
            Possible_Cover_IM[2*v][2*e]=E.IM[v][e]
            
    #Covering Information
    for e in range(configs["edges"]):
        # Edge is covered
        if(CE[e]):
	    if(CV[E.start(e)]):
                Possible_Cover_IM[2*E.start(e)+1][2*e+1]=-1
	    else:
		Possible_Cover_IM[2*E.start(e)][2*e+1]=-1
	    if(CV[E.end(e)]):
	        Possible_Cover_IM[2*E.end(e)+1][2*e+1]=1
	    else:
		Possible_Cover_IM[2*E.end(e)][2*e+1]=1

    #No loops in this base graph
    return Possible_Cover_IM

def get_Possible_Cover_IM(E, CV, CE, ccl_edges, ccl_edges_iscrossed, ccl_edges_locations):
    	Possible_Cover_IM = np.zeros((2*configs["verts"],2*configs["edges"]), dtype=np.int)
    	#copy base graph
    	for v in range(configs["verts"]):
        	for e in range(configs["edges"]):
       			Possible_Cover_IM[2*v][2*e]=E.IM[v][e]
       	# e is a loop and must be crossed (no loops in cover graphs)
        loop_num = 0    
	graph_loops = E.G.loop_edges()
        for loop in graph_loops:
        	vertex = loop[0]
                Possible_Cover_IM[2*vertex][2*loop_num]=-1
		Possible_Cover_IM[2*vertex+1][2*loop_num+1]=-1
		Possible_Cover_IM[2*vertex][2*loop_num+1]=1
	        Possible_Cover_IM[2*vertex+1][2*loop_num]=1
		loop_num+=1
    	#Covering Information
    	for e in range(configs["loops"],configs["edges"]):
        	edge_start = E.start(e)
        	edge_end = E.end(e)
        	# Edge is part of ccl and should be crossed
        	if(ccl_edges[e] and ccl_edges_iscrossed[ccl_edges_locations.index(e)]):
	    		Possible_Cover_IM[2*edge_start][2*e]=-1
	    		Possible_Cover_IM[2*edge_start+1][2*e+1]=-1
	    		Possible_Cover_IM[2*edge_end][2*e+1]=1
	    		Possible_Cover_IM[2*edge_end+1][2*e]=1
	    		Possible_Cover_IM[2*edge_end][2*e]=0

        	elif(CE[e]):
    	    		if(CV[edge_start]):
                		Possible_Cover_IM[2*edge_start+1][2*e+1]=-1
            		else:
                		Possible_Cover_IM[2*edge_start][2*e+1]=-1
            		if(CV[edge_end]):
                		Possible_Cover_IM[2*edge_end+1][2*e+1]=1
            		else:
        			Possible_Cover_IM[2*edge_end][2*e+1]=1
    	return Possible_Cover_IM 


def edge_start(M,e):
	for v in range(2*configs["verts"]):
        	if(M[v][e]==-1):
                	return v
        return -1
def edge_end(M,e):
        for v in range(2*configs["verts"]):
        	if(M[v][e]==1):
                	return v
        return -1

def is_Cover_no_loops(PCIM, CVerts, CEdges):
	start1=-1
	start2=-1
	end1=-1
	end2=-1
	for e in range(configs["edges"]):
		if(e>=configs["loops"] and CEdges[e]):
			start1=edge_start(PCIM,2*e+1)
			end1=edge_end(PCIM,2*e+1)
			if(CVerts[edge_start(PCIM,2*e)/2]):
				start2=edge_start(PCIM,2*e)+1
			else:
				start2=edge_start(PCIM,2*e)
			if(CVerts[edge_end(PCIM,2*e)/2]):
				end2=edge_end(PCIM,2*e)+1
			else:
				end2=edge_end(PCIM,2*e)
			if(start1!=start2 or end1!=end2):
				return false
		elif(e>=configs["loops"]):
			start1=edge_start(PCIM,2*e)
			end1=edge_end(PCIM,2*e)
			if(CVerts[edge_start(PCIM,2*e)/2]):
				start2=edge_start(PCIM,2*e)+1
			else:
				start2=edge_start(PCIM,2*e)
			if(CVerts[edge_end(PCIM,2*e)/2]):
				end2=edge_end(PCIM,2*e)+1
			else:
				end2=edge_end(PCIM,2*e)
			if(start1!=start2 or end1!=end2):
				return false
		start1 = -1
		start2=-1
		end1=-1
		end2=-1
		#We have constructed the loops such that the will always satisfy the cover conditions. 
	return true

def is_Cover(PCIM, CVerts, CEdges, CCL_edge_Locations, ccl_isCrossed):
	start1=-1
	start2=-1
	end1=-1
	end2=-1
	#create isCrossed
	isCrossed = [false]*configs["edges"]
	for e in range(configs["edges"]):
		if(e in CCL_edge_Locations):
			if(ccl_isCrossed[CCL_edge_Locations.index(e)]):
				isCrossed[e]=True
	for e in range(configs["edges"]):
		if(not isCrossed[e] and e>=configs["loops"]):
			if(CEdges[e]):
				start1=edge_start(PCIM,2*e+1)
				end1=edge_end(PCIM,2*e+1)
				if(CVerts[edge_start(PCIM,2*e)/2]):
					start2=edge_start(PCIM,2*e)+1
				else:
					start2=edge_start(PCIM,2*e)
				if(CVerts[edge_end(PCIM,2*e)/2]):
					end2=edge_end(PCIM,2*e)+1
				else:
					end2=edge_end(PCIM,2*e)
				if(start1!=start2 or end1!=end2):
					return false
			else:
				start1=edge_start(PCIM,2*e)
				end1=edge_end(PCIM,2*e)
				if(CVerts[edge_start(PCIM,2*e)/2]):
					start2=edge_start(PCIM,2*e)+1
				else:
					start2=edge_start(PCIM,2*e)
				if(CVerts[edge_end(PCIM,2*e)/2]):
					end2=edge_end(PCIM,2*e)+1
				else:
					end2=edge_end(PCIM,2*e)
				if(start1!=start2 or end1!=end2):
					return false

		start1=-1
		start2=-1
		end1=-1
		end2=-1
		#We have constructed the loops such that the will always satisfy the cover conditions. 
		#Crossed this.edges will also satisfy the cover criteria by construction
	return true


#----------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------Generate All Cover Graphs--------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

#Load all base Graphs
All_EGraphs = LoadGraphs()
print("Done Loading Base Graphs...\n")

APEC = []
FE = [false]*configs["edges"]
for l in range(configs["loops"]):
    FE[l] = True
APEC.append(FE)
CEdges = list(FE)
while not (is_true_after(CEdges,configs["loops"]-1)):
    Iterate_Boolean_List(CEdges,0,configs["edges"])
    CE = list(CEdges)
    APEC.append(CE)
    
APVC = []
FV = [False]*configs["verts"]
APVC.append(FV)
CVerts = list(FV)
while not (is_true_after(CVerts,-1)):
    Iterate_Boolean_List(CVerts,0,configs["verts"])
    CV = list(CVerts)
    APVC.append(CV)
    
All_CGs = []
output_file = str(configs["verts"])+"V"+str(configs["edges"])+"E"+str(configs["loops"])+"L_FS_output.txt"
second_output = str(configs["verts"])+"V"+str(configs["edges"])+"E"+str(configs["loops"])+"L_FS_Graphs.txt"

with open(output_file, "w") as f:
	with open(second_output, "w") as g:
		count = 0
		for E in All_EGraphs:
			All_CGS_For_G=[]
			for vc in APVC:
			    for ec in APEC:
				E.generate_CC_loops(vc,ec)
				#There cannot be any loops in cover graph therefore all loop vertices must be covered. 
				#Recall all loops are already covered by construction of APEC
				if(E.are_LVs_covered):
					if(len(E.CC_loops)>0):
						#Determine which edges are a part of CC_loops
						ccl_edges = [False]*configs["edges"]
						ccl_edge_locations =[]
						for bv in E.CC_loops:
							for e in range(configs["edges"]):
								if(bv[e]!=0 and not ccl_edges[e]):
									ccl_edges[e] = True
									ccl_edge_locations.append(e)
					    
						ccl_edges_iscrossed = [false]*len(ccl_edge_locations)
						#make sure loops are always crossed
						for i in range(len(ccl_edge_locations)):
							if(ccl_edge_locations[i]<configs["loops"]):
								#We are assuming that loops come first in ccl_edge_locations
								ccl_edges_iscrossed[i] = True
						Possible_IM = get_Possible_Cover_IM(E, vc, ec, ccl_edges, ccl_edges_iscrossed, ccl_edge_locations)
						if(is_Cover(Possible_IM,vc,ec,ccl_edge_locations,ccl_edges_iscrossed)):
							CG = CoverGraph(Possible_IM,ec,vc,configs)
							if(CG.G.is_connected() and CG.FS and CG.Basis_Linear_Forms.nrows()>0):
								is_new = True
								for cg in All_CGS_For_G:
									if(CG.G.is_isomorphic(cg.G)):
										is_new = False
										break
								if(is_new):
									All_CGS_For_G.append(CG)
									f.write(str(CG.Basis_Linear_Forms))
									f.write("\n\n")
									g.write(CG.output())
						while(not is_true_after(ccl_edges_iscrossed,-1)):
							Iterate_Boolean_List(ccl_edges_iscrossed,-1,len(ccl_edges_iscrossed))
							#make sure there are less than or equal to number of crossings per numer of ccls
							if(ccl_edges_iscrossed.count(True)<= len(E.CC_loops)):
								Possible_IM = get_Possible_Cover_IM(E, vc,ec, ccl_edges, ccl_edges_iscrossed, ccl_edge_locations)
								if(is_Cover(Possible_IM,vc,ec,ccl_edge_locations,ccl_edges_iscrossed)):
									CG = CoverGraph(Possible_IM,ec,vc,configs)
									if(CG.G.is_connected() and CG.FS and CG.Basis_Linear_Forms.nrows()>0):
										is_new = True
										for cg in All_CGS_For_G:
											if(CG.G.is_isomorphic(cg.G)):
												is_new = False
												break
										if(is_new):
											All_CGS_For_G.append(CG)
											f.write(str(CG.Basis_Linear_Forms))
											f.write("\n\n")
											g.write(CG.output())
					#No Loops
					else:
						Possible_IM = get_Possible_Cover_IM_no_loops(E,vc,ec)
						if(is_Cover_no_loops(Possible_IM,vc,ec)):
							CG = CoverGraph(Possible_IM,ec,vc,configs)
							if(CG.G.is_connected() and CG.FS and CG.Basis_Linear_Forms.nrows()>0):
								is_new = True
								for cg in All_CGS_For_G:
									if(CG.G.is_isomorphic(cg.G)):
										is_new = False
										break
								if(is_new):
									All_CGS_For_G.append(CG)
									f.write(str(CG.Basis_Linear_Forms))
									f.write("\n\n")
									g.write(CG.output())
			All_CGs.extend(All_CGS_For_G)
			count+=1
			num = float(count)/len(All_EGraphs)
			percent = int(num*100)
			string = str(percent) + "%"
			print(string)
print(str(len(All_CGs)) + " Cover Graphs\n")
