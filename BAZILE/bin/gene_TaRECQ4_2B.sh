### gene RECQ4_2B
#### orthologues de AtRECQ4: TraesCS2A02G304900, TraesCS2B02G321700 et TraesCS2D02G303500

## corresp Renan CS
#grep TraesCS2B02G321700 TaeRenan_refseqv2.1_CORRESPONDANCE_CSv1_CSv2.txt
TraesCS2B03G0837200     TraesRN2B01G00896500    TraesCS2B02G321700

### Renan v2.1
#grep TraesCS2B03G0837200 /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/current/genes/annot_v2.1/*HC.gff3 |grep $'\t''CDS' |head -n1
chr2B   WHEATOMICS_MAGATT_TRIANNOT      CDS     464987973       464988402       100     -       1       ID=TraesRN2B01G00896500.1.CDS1;Parent=TraesRN2B01G00896500.1;Name=TraesCS2B03G0837200.1.cds1;Target=TraesCS2B03G0837200.1 3532 3961 +


### Renan v2.0
#grep TraesRN2B0100862400 /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/current/genes/annot_v2.0/*HCgtf |grep $'\t''CDS' |head -n1
TaeRenan_refseq_v2.0_genesHC.gtf:chr2B  WHEATOMICS_MAGATT_TRIANNOT      CDS     464987973       464988402       .       -       1       transcript_id "TraesRN2B0100862400.1"; gene_id "TraesRN2B0100862400"; gene_name "TraesRN2B0100862400";


### Chinese spring
#grep TraesCS2B02G321700 /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_IDmaping.csv
TraesCS2B02G321700      TraesCS2B03G0837200