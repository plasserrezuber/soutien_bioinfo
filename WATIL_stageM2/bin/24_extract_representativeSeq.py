#!/bin/python3

from collections import defaultdict
import re,argparse

clust=[]
seq=[]
dseq=defaultdict(list)
dclust=defaultdict(int)

def parseOptions():
	#use argparse to setup/get options from shell cmd line
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", help="fichier fasta", required=True)
	return(parser.parse_args())
parameters=parseOptions()


#enregistrement dico cluster_nb_de_seq et cluster_noms_des_seq_du_cluster
with open(parameters.input,"r") as f:
	for li in f:
		l=li.rstrip("\n")
		if l.startswith(">"):
			#print("###################")
			clust=l.lstrip(">")
			#print(clust)

		elif not l.startswith(">"):
			seq=l.split()[2].lstrip(">").rstrip("...")
			#print(seq)
			dseq[clust].append(seq)
			
			dclust[clust]=l.split()[0]
			
			#print(dclust)
			#print(dseq)
            
#identification cluster majoritaire
nbseq=0
for k,v in dclust.items():
	if int(v) > int(nbseq):
		major_clust=k
		nbseq=v
print(major_clust)

#extraction liste des noms des sequences de ce cluster majoritaire
with open(re.sub(".clstr", ".txt", parameters.input),"w") as fo:
	for s in dseq[major_clust]:
		fo.write(s+'\n')
