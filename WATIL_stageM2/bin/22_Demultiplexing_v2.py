#!/usr/bin/env python
# -*-coding:Latin-1 -*

import re
import os, sys
import glob
from collections import defaultdict
import collections

primer = [line for i, line in enumerate(open(sys.argv[2])) if i % 2 == 1 ]
for i in range(len(primer)):
	primer[i] = primer[i].upper()

filename=os.path.splitext(sys.argv[1])[0]

#print(filename) 
	
lol=[]
for i in primer:
	lol.append(i.strip())
#print(lol)
	
Glu_HPM_1Ax=[]
Glu_HPM_1Bx=[]
Glu_HPM_1Dx=[]
Glu_HPM_1Ay=[]
Glu_HPM_1By=[]
Glu_HPM_1Dy=[]
Glu_FPM_i_1_A=[]
Glu_FPM_i_2_A=[]
Glu_FPM_m_3_A=[]
Glu_FPM_m_4_B=[]
Glu_FPM_m_5_B=[]
Glu_FPM_m_1_D=[]
Glu_FPM_m_3_D=[]
Glu_FPM_m_4_D=[]
Glu_FPM_m_5_D=[]
Glu_FPM_m_6_D=[]
Glu_FPM_m_7_D=[]
Glu_FPM_m_8_D=[]

with open(sys.argv[1],"r") as my_file:
	data=my_file.read()
	data_split=data.split("\n")
	for i in range(len(data_split)):
		if i % 2:
			line = [data_split[i-1], data_split[i]]
			start=line[1][0:55]
			tail=line[1][-55:]
			if ((lol[0] in start and lol[2] in tail) or (lol[1] in start and lol[3] in tail)):
				Glu_HPM_1Ax.append(line)					
			elif ((lol[4] in start and lol[6] in tail) or (lol[5] in start and lol[7] in tail)):
				Glu_HPM_1Bx.append(line)
			elif ((lol[8] in start and lol[10] in tail) or (lol[9] in start and lol[11] in tail)):		
  				Glu_HPM_1Dx.append(line)
			elif ((lol[12] in start and lol[14] in tail) or (lol[13] in start and lol[15] in tail)):
				Glu_HPM_1Ay.append(line)
			elif ((lol[16] in start and lol[18] in tail) or (lol[17] in start and lol[19] in tail)):
				Glu_HPM_1By.append(line)
			elif ((lol[20] in start and lol[22] in tail) or (lol[21] in start and lol[23] in tail)):
				Glu_HPM_1Dy.append(line)
			elif ((lol[24] in start and lol[26] in tail) or (lol[25] in start and lol[27] in tail)):
				Glu_FPM_i_1_A.append(line)
			elif ((lol[28] in start and lol[30] in tail) or (lol[29] in start and lol[31] in tail)):
				Glu_FPM_i_2_A.append(line)
			elif ((lol[32] in start and lol[34] in tail) or (lol[33] in start and lol[35] in tail)):
				Glu_FPM_m_3_A.append(line)
			elif ((lol[36] in start and lol[38] in tail) or (lol[37] in start and lol[39] in tail)):
				Glu_FPM_m_4_B.append(line)
			elif ((lol[40] in start and lol[42] in tail) or (lol[41] in start and lol[43] in tail)):
				Glu_FPM_m_5_B.append(line)
			elif ((lol[44] in start and lol[46] in tail) or (lol[45] in start and lol[47] in tail)):
				Glu_FPM_m_1_D.append(line)
			elif ((lol[48] in start and lol[50] in tail) or (lol[49] in start and lol[51] in tail)):
				Glu_FPM_m_3_D.append(line)
			elif ((lol[52] in start and lol[54] in tail) or (lol[53] in start and lol[55] in tail)):
				Glu_FPM_m_4_D.append(line)
			elif ((lol[56] in start and lol[58] in tail) or (lol[57] in start and lol[59] in tail)):
				Glu_FPM_m_5_D.append(line)
			elif ((lol[60] in start and lol[62] in tail) or (lol[61] in start and lol[63] in tail)):
				Glu_FPM_m_6_D.append(line)
			elif ((lol[64] in start and lol[66] in tail) or (lol[65] in start and lol[67] in tail)):
				Glu_FPM_m_7_D.append(line)
			elif ((lol[68] in start and lol[70] in tail) or (lol[69] in start and lol[71] in tail)):
				Glu_FPM_m_8_D.append(line)

dico={}

dico[filename + "_Glu_HPM_1Ax"]=Glu_HPM_1Ax
dico[filename + "_Glu_HPM_1Bx"]=Glu_HPM_1Bx
dico[filename + "_Glu_HPM_1Dx"]=Glu_HPM_1Dx
dico[filename + "_Glu_HPM_1Ay"]=Glu_HPM_1Ay
dico[filename + "_Glu_HPM_1By"]=Glu_HPM_1By
dico[filename + "_Glu_HPM_1Dy"]=Glu_HPM_1Dy
dico[filename + "_Glu_FPM_i_1_A"]=Glu_FPM_i_1_A
dico[filename + "_Glu_FPM_i_2_A"]=Glu_FPM_i_2_A
dico[filename + "_Glu_FPM_m_3_A"]=Glu_FPM_m_3_A
dico[filename + "_Glu_FPM_m_4_B"]=Glu_FPM_m_4_B
dico[filename + "_Glu_FPM_m_5_B"]=Glu_FPM_m_5_B
dico[filename + "_Glu_FPM_m_1_D"]=Glu_FPM_m_1_D
dico[filename + "_Glu_FPM_m_3_D"]=Glu_FPM_m_3_D
dico[filename + "_Glu_FPM_m_4_D"]=Glu_FPM_m_4_D
dico[filename + "_Glu_FPM_m_5_D"]=Glu_FPM_m_5_D
dico[filename + "_Glu_FPM_m_6_D"]=Glu_FPM_m_6_D
dico[filename + "_Glu_FPM_m_7_D"]=Glu_FPM_m_7_D
dico[filename + "_Glu_FPM_m_8_D"]=Glu_FPM_m_8_D


dico = {key:val for key, val in dico.items()}
#print(dico)

for k,v in dico.items():
        print(k, len([item for item in v]))


for (k,v) in dico.items():
	with open("{}_test.fasta".format(k), "w") as f:
		csv='\n'.join(['\n'.join(t) for t in v])
		f.write("%s\n" % (csv))  

		
