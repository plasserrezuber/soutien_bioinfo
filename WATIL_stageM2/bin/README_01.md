##Mail Vincent Pailler 23 fevrier 2022 
<br>
Les commandes suivantes sont utilisables avec la suite d'outils smrttools/10.1.0.119588 disponible sur HPC2 .
Doc disponible [https://www.pacb.com/wp-content/uploads/SMRT_Tools_Reference_Guide_v10.1.pdf] {: #https://www.pacb.com/wp-content/uploads/SMRT_Tools_Reference_Guide_v10.1.pd}
On est parti d'un pool d'amplicons/gènes multiplexés sur 96 cultivars de blé. A partir des données brutes (=subreads).

###CCS
J'ai généré les données CCS : avec le fichiers de subreads en input et le fichier ccs.bam en output:
```
module load smrttools/10.1.0.119588
ccs -j 128 --min-passes 3 *.subreads.bam ccs.bam
```
J'ai obtenu:
number of CCS:     2129864
total CCS length:  13198910525
mean CCS size:      6197.07

Tu peux télécharger sur ton home avec ce lien : 
[https://gentyane-share.s3.mesocentre.uca.fr/E779/ccs-E779-amplicons.bam?AWSAccessKeyId=33e6912b6a9f4ee7a0f48866691375ba&Expires=1646228054&Signature=PZSKDBTh%2FRgv3Pe0hZfUP%2F5i%2BJI%3D]
{: #https://gentyane-share.s3.mesocentre.uca.fr/E779/ccs-E779-amplicons.bam?AWSAccessKeyId=33e6912b6a9f4ee7a0f48866691375ba&Expires=1646228054&Signature=PZSKDBTh%2FRgv3Pe0hZfUP%2F5i%2BJI%3D}
<br>
###LIMA
J'ai ensuite démultiplexé les données selon les 96 cultivars avec en input ton fichier de reads CCS , ton fichier de barcodes .fasta (n=96) et en output tes 96 potentiels échantillons (demux.bc20XX--bc20XX.bam). Pour les paramètres, je te laisse regarder la doc.

```
module load smrttools/10.1.0.119588
lima --split-bam-named -j 32 --ccs --peek-guess  ccs.bam /home/vipailler/Barcodes/SMRTbell_Barcoded_Adapter_Plate_3.0_bc2001-bc2096.fasta demux.bam
```

Le premier niveau de démultiplexage selon les 96 cultivars (94 réellement retrouvés) est fait. Maintenant on s'intéresse aux amplicons/gènes ; soit 18 amplicons * 96 cultivars . Mireille m'a donné la liste des primers qu'elle a utilisé pour les PCRs.
<br>
###Python-EMBOSS ou PBAA

Methode Python-EMBOSS: Ici j'ai fait un script python (que tu pourras largement modifier pour mieux l'automatiser si tu veux) , qui va chercher en 5' et / ou 3' la séquence des primers (+ taille attendue de l'amplicon) sur les différents reads de tes 94 échantillons. Tu as juste à renseigner en arguments le fichier de reads fasta de chaque échantillon + le fichier de primers. C'est le résultat du tableur Excel de ce matin. 

Methode PBAA
