from Classes import *

configs = {"edges":5, "verts":3, "loops": 0}

Possible_Cover_IM = np.zeros((6, 10), dtype=np.int)
Possible_Cover_IM[0][0] = -1
Possible_Cover_IM[0][1] = 1
Possible_Cover_IM[1][0] = 1
Possible_Cover_IM[1][1] = -1
Possible_Cover_IM[0][2] = -1
Possible_Cover_IM[3][2] = 1
Possible_Cover_IM[1][3] = -1
Possible_Cover_IM[2][3] = 1
Possible_Cover_IM[0][4] = -1
Possible_Cover_IM[4][4] = 1
Possible_Cover_IM[1][5] = -1
Possible_Cover_IM[5][5] = 1
Possible_Cover_IM[2][6] = -1
Possible_Cover_IM[5][6] = 1
Possible_Cover_IM[3][7] = -1
Possible_Cover_IM[4][7] = 1
Possible_Cover_IM[2][8] = -1
Possible_Cover_IM[4][8] = 1
Possible_Cover_IM[3][9] = -1
Possible_Cover_IM[5][9] = 1

print Possible_Cover_IM

CV = [True,True, True]
CE = [True, True, True, True, True]

CG = CoverGraph(Possible_Cover_IM, CE, CV, configs)

print "Homology Basis: "
print CG.Homology_Basis
print "Image_chb"
print CG.Image_CHB
print "Basis_Linear_Forms"
print CG.Basis_Linear_Forms
print "Check FS:"
print CG.FS
