#!/bin/bash

OUTPUT='/home/palasser/soutien_bioinfo/BOUCHET'
mkdir $OUTPUT

## DEMANDE: positions physiques approximatives correspondante sur CS

# arianalrfor 7B  9583954 /home/palasser/results/bwa/isbp_RENAN_PSEUDOV2_vs_10GENOMES/ISBPS_RENAN_PSEUDOV2_vs_arianaLrFor_input_dotplot.txt
# landmark    1A 101110372    /home/palasser/results/bwa/isbp_RENAN_PSEUDOV2_vs_10GENOMES/ISBPS_RENAN_PSEUDOV2_vs_CDC_Landmark_input_dotplot.txt
# RENAN   3B    8000000   /home/palasser/results/bwa/isbp_REFSEQV2_vs_RENAN_PSEUDOV2/ISBPS_REFSEQV2_RENAN_PSEUDOV2_input_dotplot.txt
# RENAN   2D    22107109  /home/palasser/results/bwa/isbp_REFSEQV2_vs_RENAN_PSEUDOV2/ISBPS_REFSEQV2_RENAN_PSEUDOV2_input_dotplot.txt


echo -e "accession\tchrom\tposition\tisbp\tchrom_REFSEQV2\tposition_REFSEQV2\tchrom_RENAN\tposition_RENAN" > corresp_position_phys_genomes_vs_CS_RESULTS_220530.tsv

isbp_arina=$(gawk '{if ($2~"7B" && $5>=9583954-100000 && $5<=9583954+100000) print $1 }' /home/palasser/results/bwa/isbp_RENAN_PSEUDOV2_vs_10GENOMES/ISBPS_RENAN_PSEUDOV2_vs_arianaLrFor_input_dotplot.txt)
gawk -v OFS='\t' -v isbp=$isbp_arina '{if ($1==isbp) print "arianalrfor","7B","9583954",$0 }' /home/palasser/results/bwa/isbp_REFSEQV2_vs_RENAN_PSEUDOV2/ISBPS_REFSEQV2_RENAN_PSEUDOV2_input_dotplot.txt \
>> $OUTPUT/corresp_position_phys_genomes_vs_CS_RESULTS_220530.tsv

isbp_landmark=$(gawk '{if ($2~"1A" && $5>=101110372-500 && $5<=101110372+500) print $1 }' /home/palasser/results/bwa/isbp_RENAN_PSEUDOV2_vs_10GENOMES/ISBPS_RENAN_PSEUDOV2_vs_CDC_Landmark_input_dotplot.txt |tail -n1)
gawk -v OFS='\t' -v isbp=$isbp_landmark '{if ($1==isbp) print "landmark","1A","101110372",$0 }' /home/palasser/results/bwa/isbp_REFSEQV2_vs_RENAN_PSEUDOV2/ISBPS_REFSEQV2_RENAN_PSEUDOV2_input_dotplot.txt \
>> $OUTPUT/corresp_position_phys_genomes_vs_CS_RESULTS_220530.tsv

gawk -v OFS='\t' '{if ($2~"3B" && $5>=8000000-100000 && $5<=8000000+100000) print "RENAN","3B","8000000",$0 }' /home/palasser/results/bwa/isbp_REFSEQV2_vs_RENAN_PSEUDOV2/ISBPS_REFSEQV2_RENAN_PSEUDOV2_input_dotplot.txt \
>> $OUTPUT/corresp_position_phys_genomes_vs_CS_RESULTS_220530.tsv

gawk -v OFS='\t' '{if ($2~"2D" && $5>=22107109-120000 && $5<=22107109+120000) print "RENAN","2D","22107109",$0 }' /home/palasser/results/bwa/isbp_REFSEQV2_vs_RENAN_PSEUDOV2/ISBPS_REFSEQV2_RENAN_PSEUDOV2_input_dotplot.txt \
>> $OUTPUT/corresp_position_phys_genomes_vs_CS_RESULTS_220530.tsv