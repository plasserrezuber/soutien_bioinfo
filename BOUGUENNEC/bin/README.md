# TE annotation

Transposable elements annotation of *Secale cereale* genome, in collaboration with A. BOUGUENNEC (UMR GDEC, INRAE).  
from /home/login/projects/soutien_bioinfo/BOUGUENNEC/ on HPC2 cluster (Mésocentre, Université Clermont-Auvergne)

## Material & Methods  

Reference genome: [Rabanus-Wallace et al, 2021](https://doi.org/10.1038/s41588-021-00807-0)  

### CLARI-TE_ sh

A shell pipeline for TE annotation working on HPC cluster (slurm) with job dependencies + job-array and using RepeatMasker and ClariTE tools.

### CLARI-TE_ smk

A snakemake pipeline for TE annotation working on HPC cluster (slurm).  
Are required: gcc/8.1.0  python/3.7.1 snakemake/5.25.0  

Fill in the "config.yml" file with your own genome informations, and the library that RepeatMasker will use if you desire to use another one.  
Eventually custom the cluster parameters in "hpc2_ressources.json" file according to the HPC cluster used.  

To launch the smk pipeline: 
```console
snakemake --use-conda -j 20 --cluster-config hpc2_ressources.json --cluster --latency-wait 30
```

A visualization of the rules executed by the smk pipeline:
![rulegraph](rulegraph.png)


## Results  

![Table](table.png)

## Support  

pauline.lasserre-zuber@inrae.fr, frederic.choulet@inrae.fr  

## Roadmap  

Describe TE length proportion changes for 10 top TE-families compared to reference sub-genomes.

## Authors and acknowledgment  

Pauline LASSERRE-ZUBER (INRAe), Frederic CHOULET (INRAe)  