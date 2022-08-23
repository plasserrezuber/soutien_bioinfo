#!/bin/python
#sbatch --export=ALL -n1 --wrap='./fasta_per_genome_to_fasta_per_gene.py'
from Bio import SeqIO
import re

# pour ajouter Chinese-spring aux 10 genomes
# sed 's/>lcl/>/' gli_glu_iwgsc_refseqv1.0_FuncAnnot_v1.1_completed_v1.0.fna > gli_glu_refseqv1.0_FuncAnnot_gmap_sur_CS.fna

genomes=['CS', 'arinaLrFor', 'jagger', 'julius', 'Lancer', 'landmark', 'mace', 'Norin61', 'spet', 'stanley', 'SY_mattis']

with open('/home/palasser/data/gli_glu_iwgsc_refseqv1.0_FuncAnnot_v1.1_completed_v1.0.TAB', 'r') as f1:
    for li in f1.readlines():
        li=li.rstrip()
        match=re.search('(^Traes\w+)\s[T-].+', li)
        if match:
            genec=match.group(1)+'.1'
            with open('/home/palasser/results/bedtools/tenGenomes/bygene/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_10Genomes_'+genec+'.fna', 'w') as FNAoutput:
                for g in genomes:
                    listgenes=[]
                    for record in SeqIO.parse('/home/palasser/results/bedtools/tenGenomes/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_'+g+'.fna', 'fasta'):
                        if record.id==genec and genec not in listgenes:
                            record.id=g+'|'+record.id
                            SeqIO.write(record, FNAoutput, "fasta")
                            listgenes.append(genec)
