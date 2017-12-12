from sage.all import *
import pdb
import numpy as np
#--------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------Classes---------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------   
 
class EGraph:
    def __init__(self,G,IM,configs):
        self.G = G
	self.IM = IM
	self.configs = configs
        self.Homology_Basis = Matrix(self.IM).right_kernel().basis()
        self.CC_loops = []
        self.are_LVs_covered = True

    def generate_CC_loops(self,CV,CE):
        loop_verts = self.G.loop_vertices()
        self.are_LVs_covered = True
        for v in loop_verts:
            if(not CV[v]):
                self.are_LVs_covered = False
                break
        flag = True
        #For each basis vector in the homology basis
        for bv in self.Homology_Basis:
            flag = True
            #For each edge in the basis vector
            for e in range(self.configs["edges"]):
                #if the edge is actually a part of the cycle
                if(bv[e]!= 0):
                    #verify that the edge is covered 
                    condition_1 = not CE[e]
                    #verify that if the edge is a loop then the vertex is covered
                    condition_2 = (e<self.configs["loops"]) and not self.are_LVs_covered
                    #verify that if edge not a loop then the starting vertex of the edge is covered
                    condition_3 = (e>=self.configs["loops"]) and not CV[self.start(e)]
                    #verify that if edge not a loop then the ending vertex of the edge is covered
                    condition_4 = (e>=self.configs["loops"]) and not CV[self.end(e)]
                    if(condition_1 or condition_2 or condition_3 or condition_4):
                        flag = False
                        break;
            if(flag):
                self.CC_loops.append(bv)
                
    def start(self,e):
        for v in range(self.configs["verts"]):
            if(self.IM[v][e]==-1):
                return v
        return -1
    def end(self,e):
        for v in range(self.configs["verts"]):
            if(self.IM[v][e]==1):
                return v
        return -1

#--------------------------------------------------------------------------------------------------------------------------------------------------
class CoverGraph:
	def __init__(self,IM,CE,CV,configs):
		self.configs = configs
		self.IM = IM
		self.CE = CE
		self.CV = CV
		self.sage_IM = np.copy(self.IM)
		#remove zero rows
		self.sage_IM=self.sage_IM[~np.all(self.sage_IM==0,axis=1)]
		#remove zero cols
		self.sage_IM=self.sage_IM[:,~np.all(self.sage_IM==0,axis=0)]
                self.G = Graph(Matrix(self.sage_IM))   
		self.get_Homology_basis()
		self.get_Image_CHB()
		self.find_basis_Linear_Forms()
		self.check_FS()
	def get_Homology_basis(self):
		self.Homology_Basis = list(Matrix(self.IM).right_kernel().basis())
		for bv in self.Homology_Basis:
			count = 0 
			for e in range(2*self.configs["edges"]):
				if(bv[e]==0):
					count+=1
			#remove uncovered edges
			if(count==2*self.configs["edges"]-1):
				self.Homology_Basis.remove(bv)
	def get_Image_CHB(self):
		#get identity - involution matrix
		Identity_Minus_Edge_Involution = np.zeros((2*self.configs["edges"],2*self.configs["edges"]),dtype=np.int)
		for i in range(self.configs["edges"]):
			if(self.CE[i]):
				Identity_Minus_Edge_Involution[2*i+1][2*i]=-1
				Identity_Minus_Edge_Involution[2*i][2*i+1]=-1
				Identity_Minus_Edge_Involution[2*i][2*i]=1
				Identity_Minus_Edge_Involution[2*i+1][2*i+1]=1
		HB_matrix = np.matrix(self.Homology_Basis)
		Unreduced_Image = np.dot(HB_matrix,Identity_Minus_Edge_Involution)
		self.Image_CHB = Matrix(ZZ,Unreduced_Image).echelon_form()
		zero_rows = []
		for row in range(self.Image_CHB.nrows()):
			is_zero = True
			for col in range(2*self.configs["edges"]):
				if(self.Image_CHB[row,col]!=0):
					is_zero = False
					break
			if(is_zero):
				zero_rows.append(row)
		if(len(zero_rows)>0):
			self.Image_CHB = self.Image_CHB.delete_rows(zero_rows)
	def find_basis_Linear_Forms(self):
		self.Basis_Linear_Forms = []
		for j in range(self.Image_CHB.nrows()):
			b_vector = [0]*self.configs["edges"]
			for i in range(self.configs["edges"]):
				if(self.CE[i]):
					b_vector[i]=self.Image_CHB[j][2*i]-self.Image_CHB[j][2*i+1]
				else:
					b_vector[i]=self.Image_CHB[j][2*i]
			self.Basis_Linear_Forms.append(b_vector)
		#Make vectors primitive, Divide by two if possible
		for e in range(self.configs["edges"]):
			while(self.is_divisible_by_two(e)):
				for z in self.Basis_Linear_Forms:
					z[e]/=2
		self.Basis_Linear_Forms = Matrix(self.Basis_Linear_Forms)
	def is_divisible_by_two(self, e):
		flag= False
		for z in self.Basis_Linear_Forms:
			x = z[e]
			if(x%2!=0):
				return False
			#Need to make sure row is not all zeros to avoid infinite loop
			if(x!=0):
				if(not flag):
					flag=true
		return flag

	def check_FS(self):
		FS = False
		rank = self.Basis_Linear_Forms.rank()
		minors = self.Basis_Linear_Forms.minors(rank)
		units = [0,-1,1] #units and zero of ZZ
		for m in minors:
			if(m not in units):
				FS = true
		if(FS and rank>3):
			self.FS = True
		else:
			self.FS = False
	def output(self):
		string = "Incidence Matrix: \n"
		string += str(self.IM) + "\n"
		string += "Basis of Linear Forms: \n"
		string += str(self.Basis_Linear_Forms) + "\n"
		return string
