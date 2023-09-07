onsuccess:
    print("Workflow finished, no error")

CHROM = ["chr1R", "chr2R", "chr3R", "chr4R", "chr5R", "chr6R", "chr7R", "chrUn"]


onstart:
    print("##### ClariTE annotation Workflow #####\n") 
    print("## Creating output folders ##\n")
    shell('mkdir -p logs')
    shell('mkdir -p results')
    shell('mkdir -p results/expand("{chrom}", chrom=CHROM)"')


rule all:
    input:
        expand("results/{chrom}/{chrom}.gff3", chrom=CHROM)
        config["TeIDsPrefix"]+"_clariTE.gff3"


rule merge_chrom_gff3:
    output:
        config["TeIDsPrefix"]+"_clariTE.gff3"
    input:
        expand("results/{chrom}/{chrom}.gff3", chrom=CHROM)
    log:
        err="logs/gt_gff3_tidy.err"
    conda:
        "env.yml"
    shell: "gt gff3 -sort -tidy -retainids {input} 1> {output} 2> {log.err}" 


rule check_gff3:
    output:
        "logs/{chrom}_gff3validator.std"
    input:
        "results/{chrom}/{chrom}.gff3"
    log:
        std="logs/{chrom}_gff3validator.std"
    conda:
        "env.yml"
    shell:
        """
        ##sed -i "s/=TraesRe_chr/='config['TeIDsPrefix']'_chr/g" {input}
        gt gff3validator {input} 1> {log.std}
        """


rule embl_to_gff3:
    output:
        "results/{chrom}/{chrom}.gff3"
    input: ID, = glob_wildcards("results/{chrom}/{chrom}.{id}.fa.out_anno.embl")
    conda:
        "env.yml"
    shell:
        """
        bin/embl_to_gff3.sh {chrom} config['genomeFasta'] config['TeIDsPrefix']
        """


rule clariTE:
    output: EMBL, = glob_wildcards("results/{chrom}/{chrom}.{embl}.fa.out_anno.embl")
    input: ID, = glob_wildcards("results/{chrom}/{chrom}.{id}.fa")
    conda: "perl-bioperl=1.7"
    shell:
        """
        sed -i 's/#Unspecified/#Unknown/' {input}.out.xm
        bin/clariTE_1.0/bin/clariTE.pl -dir results/{chrom}/ -LTR bin/CLARIwheat.LTR_position -classi bin/CLARIwheat.classification -fasta {input} {input}.out.xm
        """


rule repeatmasker:
    output: XM, = glob_wildcards("results/{chrom}/{chrom}.{xm}.fa.out.xm")
    input: ID, = glob_wildcards("results/{chrom}/{chrom}.{id}.fa")
    threads: 16
    conda:
        "env.yml"
    shell:
        """
        RepeatMasker -e crossmatch -lib config['clariTElib'] -xsmall -nolow -xm -pa {threads} -q {input}
        """


rule chunks_fasta:
    output:
        ID, = glob_wildcards("results/{chrom}/{chrom}.{id}.fa")
    input: 
        "results/{chrom}/{chrom}.windows.fasta"
    conda:
        "env.yml"
    shell: "fastaexplode -f {input} -d results/{chrom}"


rule chunks_multi_fasta:
    output:
        "results/{chrom}/{chrom}.windows.fasta"
    input:
        "results/{chrom}/{chrom}.windows.bed"
    conda:
        "env.yml"
    shell: "bedtools getfasta -bed {input} -fi results/{chrom}/{chrom}.fasta > {output}"


rule chunks_bed:
    output:
        "results/{chrom}/{chrom}.windows.bed"
    input:
        "results/{chrom}/{chrom}.fasta.fai"
    conda:
        "env.yml"
    shell: "bedtools makewindows -g {input} -w 10000000 > {output}"


rule chrom_fai:
    output:
        "results/{chrom}/{chrom}.fasta.fai"
    input:
        "results/{chrom}/{chrom}.fasta"
    conda:
        "env.yml"
    shell: "samtools faidx {input}"


rule chrom_fasta:
    output:
        "results/{chrom}/{chrom}.fasta"
    input:
        config["genomeFasta"]
    conda:
        "env.yml"
    shell: "samtools faidx {input} {chrom} > {output}"


rule genome_fai:
    output:
        config["genomeFasta"]+".fai"
    input:
        config["genomeFasta"]
    conda:
        "env.yml"
    shell: "samtools faidx {input}"
