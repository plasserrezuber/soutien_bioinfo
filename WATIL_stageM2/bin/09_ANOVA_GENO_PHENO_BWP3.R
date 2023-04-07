library(data.table)
library(tidyverse)
library(FactoMineR)
library(broom)
library(ggplot2)
library(nlme)
library(lme4)

######################
## solution (intermédaire entre GWAS et anova/lm
# pour faire des test d'asso quand pas beaucoup d'indiv: 
# package sommmer, fonction mmer en mettant kinship comme part variable, snp effet fixes
# vignettes/tuto en ligne: https://github.com/covaruber/sommer
#############################################################

setwd('Y:/ANALYSES/data/storage_proteins_chr1/Exige_FSOV')
dirout='Y:/ANALYSES/results/storage_proteins_chr1/Exige_FSOV'


#################################################################################
## INPUT DATA
#################################################################################

#################################################################################
## DATA GENO
## INFORMATION SUR FORMAT VCF; FORMAT SNP
## avec msa2vcf: allele REF = allele MAJORITAIRE ; allele ALT = allele MINORITAIRE
## les positions sont celles de l'alignement multiple 

geno=read_tsv("GENO_PacBio_94indiv_18_glutenin_Exige.tsv", col_names = TRUE)
pheno_list=read_tsv("liste_PHENO_72_accessions_FSOV_Exige.txt", col_names = TRUE)
pheno_list$ERGE=as.character(pheno_list$ERGE)
colnames(pheno_list)=gsub("code_pacbio","sample", colnames(pheno_list))


#### NB: 71 indiv en commun, pheno ET geno
g=geno%>%filter(sample%in%pheno_list$sample)

#################################################################################
## filtres MAF et NA sur data geno

## TABLE calcul freq allelique
freq_allelic=g%>%select(-sample)%>%
  pivot_longer(everything(), names_to = "marker", values_to="allele")%>%
  group_by(marker,allele)%>%
  summarise(n=n())%>%
  mutate(freq_perc=n/nrow(g)*100)

## TABLE stat MAF et NA par marker
MAF_NA=left_join(freq_allelic%>%group_by(marker)%>%summarise(MAF=min(freq_perc)),
                 freq_allelic%>%filter(is.na(allele))%>%select(marker,freq_perc)%>%rename(NAperc=freq_perc), by="marker")

## verif nombre genes representes
nb_genes=as.numeric(MAF_NA%>%filter(MAF>=2.5 & (NAperc<=10 || is.na(NAperc)))%>%select(marker)%>%
  mutate(marker=gsub(pattern = "_[0-9]*$",replacement = "", marker))%>%distinct()%>%summarise(n()))

## LISTE markers OK
### une approche (from Renaud RINCENT) peut consister a ne pas filtrer, et regarder a posteriori la MAF du marker si asso significative
marker_ok=MAF_NA%>%filter(MAF>=2.5 & (NAperc<=10 || is.na(NAperc)))%>%select(marker)%>%distinct()%>%as_vector()

### SELECT COLONNES GENO
g=g[,c("sample", colnames(g)[colnames(g)%in%marker_ok])]


#################################################################################
### DATA STRUCTURE
#s=read_tsv("wp3_4403ind_8741hap_k8_rep5_structure.txt", col_names = TRUE)

ped=fread("exige_alixan_2020/2020_CSI_DRY/2020_CSI_DRY_filtered_pruned_recoded.ped", sep=" ", header=F)
colnames(ped)[1:2]=c("variete","ERGE")

ped=ped%>%select(-variete,-V3,-V4,-V5,-V6)
ped=data.frame(ped)
rownames(ped)=ped$ERGE
ped=ped[,-1]
ped=as.matrix(ped)


for (i in 1:ncol(ped)) {
  alleles=unique(ped[,i])
  
  for (j in 1:length(alleles)) {
    
    if (alleles[j]!=0)   { ped[ped[,i]==alleles[j],i]=j }
    
    # Replace NA by minor allele frequency computed as follow:
    if (alleles[j]==0)   { ped[ped[,i]==alleles[j],i]=NA }
  }
}

ped=data.frame(ped)

for (i in 1:ncol(ped)) {
  ped[,i]=as.numeric(ped[,i])
  ped[is.na(ped[,i]),i]=mean(x =ped[,i], na.rm=T)
}

write_tsv(ped, "2020_BWP3_GENO_filtered_pruned_NUMERIC_without_NA_for_PCA.ped")

acp<-CA(ped,graph = F)

acpdf=data.frame(acp$row$coord)
acpdf=cbind(rownames(acpdf),acpdf)
colnames(acpdf)[1]="ERGE"

acpdf=left_join(acpdf, pheno_list, by="ERGE")
acpdf[is.na(acpdf$origine),"origine"]="notEXIGE"
acpdf[is.na(acpdf$date_obt),"date_obt"]="notEXIGE"

pdf(paste(dirout,"AC_sur_BWP3_GENO_FSOV_Exige.pdf", sep="/"), width = 12, height = 8)

# eigen values barplot
barplot(acp$eig[1:20,1], ylab="Eigen values")
# percentage of explained variance barplot
barplot(acp$eig[1:20,2], ylab="% explianed variance")	

## Factorial plan graphs
# Genotypes graph

plot(acpdf[,"Dim.1"],	# Dim.1 is X axe
     acpdf[,"Dim.2"],	# Dim.2 is Y axe
     main= "Genetic diversity",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 2",
     col=c("red", "grey", "blue")[as.numeric(as.factor(acpdf$date_obt))],
     sub="red=ancien - grey=notEXIGE - blue=recent"
)
abline(h=0,v=0,lty=2)			# adding lines

plot(acpdf[,"Dim.1"],	# Dim.1 is X axe
     acpdf[,"Dim.2"],	# Dim.2 is Y axe
     main= "BWP3 Genetic diversity",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 2",
     col=c("red", "blue", "green", "black", "grey", "pink")[as.numeric(as.factor(acpdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - grey=notEXIGE - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines


plot(acpdf[,"Dim.1"],	
     acpdf[,"Dim.3"],	
     main= "BWP3 Genetic diversity",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 3",
     col=c("red", "blue", "green", "black", "grey", "pink")[as.numeric(as.factor(acpdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - grey=notEXIGE - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines

plot(acpdf[,"Dim.2"],
     acpdf[,"Dim.3"],
     main= "BWP3 Genetic diversity",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 2",
     ylab="Dim 3",
     col=c("red", "blue", "green", "black", "grey", "pink")[as.numeric(as.factor(acpdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - grey=notEXIGE - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines

dev.off()


s=acp$row$contrib[,c(1:3)]
s=cbind(rownames(s),s)
colnames(s)=gsub(" ","", colnames(s))
colnames(s)[1]="ERGE"
s=as_tibble(s)
s=s%>%filter(ERGE%in%pheno_list$ERGE)
s$Dim1=as.double(s$Dim1)
s$Dim3=as.double(s$Dim3)
s$Dim2=as.double(s$Dim2)

# sdim1=s%>%rowwise(ERGE)%>%mutate(Dim1=ifelse(Dim1==max(Dim1,Dim2,Dim3), 1, 0))%>%select(ERGE,Dim1)
# sdim2=s%>%rowwise(ERGE)%>%mutate(Dim2=ifelse(Dim2==max(Dim1,Dim2,Dim3), 1, 0))%>%select(ERGE,Dim2)
# sdim3=s%>%rowwise(ERGE)%>%mutate(Dim3=ifelse(Dim3==max(Dim1,Dim2,Dim3), 1, 0))%>%select(ERGE,Dim3)
# 
# structure=sdim1%>%left_join(sdim2)%>%left_join(sdim3)

#################################################################################
## STAT PHENO: moyenne pour indiv avec repet pour meme condition
#################################################################################

#################################################################################
## DATA PHENO
phenoA=fread("Alixan.txt", header = TRUE, sep="\t")
colnames(phenoA)=gsub("code_pacbio","sample", colnames(phenoA))

phenoM=fread("Mons.txt", header = TRUE, sep="\t")
colnames(phenoM)=gsub("code_pacbio","sample", colnames(phenoM))
colnames(phenoM)=gsub("mode","TTT", colnames(phenoM))

## declarer numero premiere col pheno
numPhenoA=12
numPhenoM=15

## vecteurs colonnes pheno
colphenoA=colnames(phenoA)[numPhenoA:ncol(phenoA)]
colphenoM=colnames(phenoM)[numPhenoM:ncol(phenoM)]


pdf(paste(dirout,"Analayses_PHENO_Correlations_FSOV_Exige.pdf", sep="/"), width = 12, height = 8)

#################################################################################
### Alixan ########################################################
phenoA_irr=cbind("IRR",phenoA%>%filter(TTT=="IRR")%>%group_by(sample)%>%summarise(across(c(ERGE,colphenoA), mean, na.rm=T)))
colnames(phenoA_irr)[1]="TTT"
row.names(phenoA_irr)=phenoA_irr$ERGE

phenoA_dry=cbind("DRY",phenoA%>%filter(TTT=="DRY")%>%group_by(sample)%>%summarise(across(c(ERGE,colphenoA), mean, na.rm=T)))
colnames(phenoA_dry)[1]="TTT"
row.names(phenoA_dry)=phenoA_dry$ERGE


### Correlations between phenotypic variables TTT=IRR
###############
cormat=round(cor(phenoA_irr[,colphenoA], method="pearson", use = "complete.obs"),2)

## replace upper triangle of the matrix by NAs
#cormat[upper.tri(cormat, diag=TRUE)]=NA

## add rownames of the data frame as new column named "var1"
var_names=colnames(data.frame(cormat))
df_cormat=cbind(var1=var_names, data.frame(cormat))

### reshape correlation matrix
r_cormat=tibble(df_cormat%>%gather(var_names,key="var2",value="pearson_coeff")) 
#r_cormat=tibble(df_cormat%>%gather(var_names,key="var2",value="pearson_coeff")%>%drop_na(pearson_coeff))

### plot correlation matrix
ggplot(r_cormat, aes(x=var1, y=var2, fill=pearson_coeff)) + 
  geom_tile()+
  #geom_text(aes(var1,var2,label=pearson_coeff), colour="white")+
  theme(axis.text.x = element_text(angle = 40, hjust=1))+
  coord_fixed()+
  ggtitle("correlations pearson Alixan TTT=IRR")


############################
### ACP QUANTITATIVE TRAITS Alixan TTT=IRR ########################################################
resPCA=PCA(phenoA_irr[,colphenoA], graph=F)
resPCAdf=data.frame(resPCA$ind$coord)
resPCAdf=cbind(rownames(resPCAdf),resPCAdf)
colnames(resPCAdf)[1]="ERGE"

resPCAdf=left_join(resPCAdf, pheno_list, by="ERGE")

# eigen values barplot
barplot(resPCA$eig[1:20,1], ylab="Eigen values")
# percentage of explained variance barplot
barplot(resPCA$eig[1:20,2], ylab="% explianed variance")	

## Factorial plan graphs
# Genotypes graph

plot(resPCAdf[,"Dim.1"],	# Dim.1 is X axe
     resPCAdf[,"Dim.2"],	# Dim.2 is Y axe
     main= "PCA PHENO Alixan 2020 IRR",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 2",
     col=c("red", "blue")[as.numeric(as.factor(resPCAdf$date_obt))],
     sub=" blue=recent - red=ancien"
)
abline(h=0,v=0,lty=2)			# adding lines

plot(resPCAdf[,"Dim.1"],	# Dim.1 is X axe
     resPCAdf[,"Dim.2"],	# Dim.2 is Y axe
     main= "PCA PHENO Alixan 2020 IRR",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 2",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines


plot(resPCAdf[,"Dim.1"],	
     resPCAdf[,"Dim.3"],	
     main= "PCA PHENO Alixan 2020 IRR",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 3",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines

plot(resPCAdf[,"Dim.2"],
     resPCAdf[,"Dim.3"],
     main= "PCA PHENO Alixan 2020 IRR",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 2",
     ylab="Dim 3",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines

plot.PCA(resPCA, choix="varcor", title="PCA PHENO Alixan 2020 IRR")

### Correlations between phenotypic variables TTT=DRY
###############
cormat=round(cor(phenoA_dry[,colphenoA], method="pearson", use = "complete.obs"),2)

## replace upper triangle of the matrix by NAs
#cormat[upper.tri(cormat, diag=TRUE)]=NA

## add rownames of the data frame as new column named "var1"
var_names=colnames(data.frame(cormat))
df_cormat=cbind(var1=var_names, data.frame(cormat))

### reshape correlation matrix
r_cormat=tibble(df_cormat%>%gather(var_names,key="var2",value="pearson_coeff")) 
#r_cormat=tibble(df_cormat%>%gather(var_names,key="var2",value="pearson_coeff")%>%drop_na(pearson_coeff))

### plot correlation matrix
ggplot(r_cormat, aes(x=var1, y=var2, fill=pearson_coeff)) + 
  geom_tile()+
  #geom_text(aes(var1,var2,label=pearson_coeff), colour="white")+
  theme(axis.text.x = element_text(angle = 40, hjust=1))+
  coord_fixed()+
  ggtitle("correlations pearson Alixan TTT= DRY")

### ACP QUANTITATIVE TRAITS Alixan TTT=DRY ########################################################
resPCA=PCA(phenoA_dry[,colphenoA], graph=F)
resPCAdf=data.frame(resPCA$ind$coord)
resPCAdf=cbind(rownames(resPCAdf),resPCAdf)
colnames(resPCAdf)[1]="ERGE"

resPCAdf=left_join(resPCAdf, pheno_list, by="ERGE")

# eigen values barplot
barplot(resPCA$eig[1:20,1], ylab="Eigen values")
# percentage of explained variance barplot
barplot(resPCA$eig[1:20,2], ylab="% explianed variance")	

## Factorial plan graphs
# Genotypes graph

plot(resPCAdf[,"Dim.1"],	# Dim.1 is X axe
     resPCAdf[,"Dim.2"],	# Dim.2 is Y axe
     main= "PCA PHENO Alixan 2020 DRY",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 2",
     col=c("red", "blue")[as.numeric(as.factor(resPCAdf$date_obt))],
     sub=" blue=recent - red=ancien"
)
abline(h=0,v=0,lty=2)			# adding lines

plot(resPCAdf[,"Dim.1"],	# Dim.1 is X axe
     resPCAdf[,"Dim.2"],	# Dim.2 is Y axe
     main= "PCA PHENO Alixan 2020 DRY",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 2",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines


plot(resPCAdf[,"Dim.1"],	
     resPCAdf[,"Dim.3"],	
     main= "PCA PHENO Alixan 2020 DRY",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 3",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines

plot(resPCAdf[,"Dim.2"],
     resPCAdf[,"Dim.3"],
     main= "PCA PHENO Alixan 2020 DRY",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 2",
     ylab="Dim 3",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines

plot.PCA(resPCA, choix="varcor", title="PCA PHENO Alixan 2020 DRY")

#################################################################################
## Mons #########################################################
phenoM_N=cbind("N",phenoM%>%filter(TTT=="N")%>%group_by(sample)%>%summarise(across(c(ERGE,colphenoM), mean, na.rm=T)))
colnames(phenoM_N)[1]="TTT"
row.names(phenoM_N)=phenoM_N$ERGE

phenoM_N0=cbind("N0",phenoM%>%filter(TTT=="N0")%>%group_by(sample)%>%summarise(across(c(ERGE,colphenoM), mean, na.rm=T)))
colnames(phenoM_N0)[1]="TTT"
row.names(phenoM_N0)=phenoM_N0$ERGE


### Correlations between phenotypic variables TTT=N
###############
cormat=round(cor(phenoM_N[,colphenoM], method="spearman", use = "complete.obs"),2)

## add rownames of the data frame as new column named "var1"
var_names=colnames(data.frame(cormat))
df_cormat=cbind(var1=var_names, data.frame(cormat))

### reshape correlation matrix
r_cormat=tibble(df_cormat%>%gather(var_names,key="var2",value="pearson_coeff")) 

### plot correlation matrix
ggplot(r_cormat, aes(x=var1, y=var2, fill=pearson_coeff)) + 
  geom_tile()+
  #geom_text(aes(var1,var2,label=pearson_coeff), colour="white")+
  theme(axis.text.x = element_text(angle = 40, hjust=1))+
  coord_fixed()+
  ggtitle("correlations pearson Mons TTT=N")

### ACP QUANTITATIVE TRAITS Mons TTT=N ########################################################
resPCA=PCA(phenoM_N[,colphenoM], graph=F)
resPCAdf=data.frame(resPCA$ind$coord)
resPCAdf=cbind(rownames(resPCAdf),resPCAdf)
colnames(resPCAdf)[1]="ERGE"

resPCAdf=left_join(resPCAdf, pheno_list, by="ERGE")

# eigen values barplot
barplot(resPCA$eig[1:20,1], ylab="Eigen values")
# percentage of explained variance barplot
barplot(resPCA$eig[1:20,2], ylab="% explianed variance")	

## Factorial plan graphs
# Genotypes graph

plot(resPCAdf[,"Dim.1"],	# Dim.1 is X axe
     resPCAdf[,"Dim.2"],	# Dim.2 is Y axe
     main= "PCA PHENO Mons 2020 N",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 2",
     col=c("red", "blue")[as.numeric(as.factor(resPCAdf$date_obt))],
     sub=" blue=recent - red=ancien"
)
abline(h=0,v=0,lty=2)			# adding lines

plot(resPCAdf[,"Dim.1"],	# Dim.1 is X axe
     resPCAdf[,"Dim.2"],	# Dim.2 is Y axe
     main= "PCA PHENO Mons 2020 N",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 2",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines


plot(resPCAdf[,"Dim.1"],	
     resPCAdf[,"Dim.3"],	
     main= "PCA PHENO Mons 2020 N",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 3",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines

plot(resPCAdf[,"Dim.2"],
     resPCAdf[,"Dim.3"],
     main= "PCA PHENO Mons 2020 N",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 2",
     ylab="Dim 3",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines

plot.PCA(resPCA, choix="varcor", title="PCA PHENO Mons 2020 N")

### Correlations between phenotypic variables TTT=N0
###############
cormat=round(cor(phenoM_N0[,colphenoM], method="pearson", use = "complete.obs"),2)

## add rownames of the data frame as new column named "var1"
var_names=colnames(data.frame(cormat))
df_cormat=cbind(var1=var_names, data.frame(cormat))

### reshape correlation matrix
r_cormat=tibble(df_cormat%>%gather(var_names,key="var2",value="pearson_coeff")) 

### plot correlation matrix
ggplot(r_cormat, aes(x=var1, y=var2, fill=pearson_coeff)) + 
  geom_tile()+
  #geom_text(aes(var1,var2,label=pearson_coeff), colour="white")+
  theme(axis.text.x = element_text(angle = 40, hjust=1))+
  coord_fixed()+
  ggtitle("correlations pearson Mons TTT=N0")

### ACP QUANTITATIVE TRAITS Mons TTT=N ########################################################
resPCA=PCA(phenoM_N0[,colphenoM], graph=F)
resPCAdf=data.frame(resPCA$ind$coord)
resPCAdf=cbind(rownames(resPCAdf),resPCAdf)
colnames(resPCAdf)[1]="ERGE"

resPCAdf=left_join(resPCAdf, pheno_list, by="ERGE")

# eigen values barplot
barplot(resPCA$eig[1:20,1], ylab="Eigen values")
# percentage of explained variance barplot
barplot(resPCA$eig[1:20,2], ylab="% explianed variance")	

## Factorial plan graphs
# Genotypes graph

plot(resPCAdf[,"Dim.1"],	# Dim.1 is X axe
     resPCAdf[,"Dim.2"],	# Dim.2 is Y axe
     main= "PCA PHENO Mons 2020 N0",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 2",
     col=c("red", "blue")[as.numeric(as.factor(resPCAdf$date_obt))],
     sub=" blue=recent - red=ancien"
)
abline(h=0,v=0,lty=2)			# adding lines

plot(resPCAdf[,"Dim.1"],	# Dim.1 is X axe
     resPCAdf[,"Dim.2"],	# Dim.2 is Y axe
     main= "PCA PHENO Mons 2020 N0",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 2",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines


plot(resPCAdf[,"Dim.1"],	
     resPCAdf[,"Dim.3"],	
     main= "PCA PHENO Mons 2020 N0",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 1",
     ylab="Dim 3",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines

plot(resPCAdf[,"Dim.2"],
     resPCAdf[,"Dim.3"],
     main= "PCA PHENO Mons 2020 N0",	# title
     pch=16,				# symbol circle
     cex=1,				# half size symbol
     asp=1,       # orthonormal basis
     xlab="Dim 2",
     ylab="Dim 3",
     col=c("red", "blue", "green", "black", "pink")[as.numeric(as.factor(resPCAdf$origine))],
     sub="red=AmNord - blue=AmSud - green=Asie - black=Europe - pink=Océanie"
)
abline(h=0,v=0,lty=2)			# adding lines

plot.PCA(resPCA, choix="varcor", title="PCA PHENO Mons 2020 N0")

dev.off()

#################################################################################
### Histogrammes Alixan

list=left_join(phenoA_irr[,c(1,2)], phenoA_dry[,c(1,2)], by="sample")[,2]
phenoA_mean=rbind(phenoA_irr%>%filter(sample%in%list),phenoA_dry%>%filter(sample%in%list))

pdf(paste(dirout,"histo_Alixan.pdf", sep="/"), width = 12, height = 8)

starts=seq(from=1, to=51, by=5)
ends=seq(from=5, to=55, by=5)

for (i in 1:length(starts)) {

dA_gathered=phenoA_mean%>%gather(colphenoA,key="trait", value="measure")%>%filter(trait%in%colphenoA[starts[i]:ends[i]])
h=ggplot(dA_gathered, aes(measure))+
  geom_histogram(stat="bin")+
  facet_grid(rows=vars(TTT), cols=vars(trait), scales="free")
plot(h)
}
dev.off()


### Correlation PHENO 2 a 2 Alixan IRR vs DRY
pdf(paste(dirout,"corr_PHENO_IRR_vs_DRY_Alixan.pdf", sep="/"), width = 12, height = 8)
par(mfrow=c(2,3))

for (t in 1:length(colphenoA)) {
  
  plot(phenoA_irr[phenoA_irr$sample%in%list,colphenoA[t]], phenoA_dry[phenoA_dry$sample%in%list,colphenoA[t]],
       main=paste("ALIXAN",colphenoA[t]))
}
dev.off()

### Histogrammes Mons

phenoM_N[,2]==phenoM_N0[,2]
phenoM_mean=rbind(phenoM_N,phenoM_N0)

pdf(paste(dirout,"histo_Mons.pdf", sep="/"), width = 12, height = 8)

starts=seq(from=1, to=51, by=5)
ends=seq(from=5, to=55, by=5)

for (i in 1:length(starts)) {
  
  dM_gathered=phenoM_mean%>%gather(colphenoM,key="trait", value="measure")%>%filter(trait%in%colphenoM[starts[i]:ends[i]])
  h=ggplot(dM_gathered, aes(measure))+
    geom_histogram(stat="bin")+
    facet_grid(rows=vars(TTT), cols=vars(trait), scales="free")
  plot(h)
}
dev.off()



### Correlation PHENO 2 a 2 Mons IRR vs DRY
pdf(paste(dirout,"corr_PHENO_N_vs_N0_Mons.pdf", sep="/"), width = 12, height = 8)
par(mfrow=c(2,3))

for (t in 1:length(colphenoM)) {
  
  plot(phenoM_N[,colphenoM[t]], phenoM_N0[,colphenoM[t]],
       main=paste("MONS",colphenoM[t]))
}
dev.off()

#################################################################################
## join data pheno avec structure

phenoA_mean$ERGE=as.character(phenoA_mean$ERGE)
psA=phenoA_mean%>%left_join(s, by="ERGE")%>%
  relocate(starts_with("Dim", ignore.case=F), .after=ERGE)

phenoM_mean$ERGE=as.character(phenoM_mean$ERGE)
psM=phenoM_mean%>%left_join(s, by="ERGE")%>%
  relocate(starts_with("Dim", ignore.case=F), .after=ERGE)

#################################################################################
## test anova ALIXAN
print(length(colphenoA))
#############################################
## TEST sur "TKW15_g","HMW1Dx_2260" ALIXAN

#############################################

for (t in 1:length(colphenoA)) {
  
  #t=2  #test
  trait=colphenoA[t]
  print(t)
  
  if (sum(is.na(psA[,t+6]))<=14) {
    
    ## join geno pheno 
    dA=left_join(psA[,c(1:6,t+6)], geno, by="sample")
    
    #test=dA[!is.na(dA$HMW1Dx_2260),c("TTT","sample","ERGE","Dim1","Dim2","Dim3","TKW15_g","HMW1Dx_2260")]
    #colnames(test)[7]="TRAIT"
    #colnames(test)[8]="allele"
    #test$allele=factor(test$allele)
    #test$TTT=factor(test$TTT)
    #summary(glm(TRAIT ~ Dim1 + Dim2 + Dim3 + TTT + allele, data=test))$coefficients
  
    dAmodel=dA%>%rename(TRAIT=colphenoA[t])%>%select(TTT, sample, TRAIT, starts_with(c("Dim","HMW","LMW"), ignore.case=F))%>%
      pivot_longer(starts_with(c("HMW","LMW")), names_to = "marker", values_to="allele")%>%
      mutate(marker = as.factor(marker))
  
  
    vec_polymorph=dAmodel%>%group_by(marker)%>%
      summarise(niv=n_distinct(allele, na.rm=T))%>%
      filter(niv>=2)%>%
      select(marker)%>%
      as_vector()
    vec_polymorph=droplevels(vec_polymorph)
  
    #correction de Bonferroni:
    seuil_signif=0.05/length(vec_polymorph)
  
    ## filtre marqueurs polymorphes
    dAmodel=dAmodel%>%filter(marker%in%vec_polymorph)   # dAmodel=dAmodel[dAmodel$marker%in%vec_polymorph,]
  
    ## declarer les variables facteurs fixes ou non fixes
    dAmodel$allele=factor(dAmodel$allele)
    dAmodel$TTT=factor(dAmodel$TTT)
    # dAmodel$Dim1=factor(dAmodel$Dim1)
    # dAmodel$Dim2=factor(dAmodel$Dim2)
    # dAmodel$Dim3=factor(dAmodel$Dim3)
  
  
    ## MODEL LINEAIRE
    #test[[1]]$coefficients
    res_lm=dAmodel%>%nest(data=-marker)%>%
      mutate(test=map(data, ~ summary(lm(TRAIT ~ Dim1 + Dim2 + Dim3 + TTT + allele, data=.x))), table=map(test, tidy))%>%
      unnest(table)%>%select(marker, term, estimate, std.error, statistic, p.value)%>%filter(term=="allele1" & p.value<=seuil_signif)
  
  
    if (dim(res_lm)[1]!=0) { res_lm=res_lm%>%cbind(trait)} else { 
      res_lm[1,]=rep(NA,7)
      res_lm=cbind(res_lm, trait) }
  
    #rbind sur la boucle
  
    if (t==1) { RES=res_lm } else {
      RES=rbind(RES,res_lm) }
  
  }

}

write_tsv(RES, paste(dirout,"Tendances_asso_Exige_Alixan.tab", sep="/"))

#################################################################################
## test anova MONS
print(length(colphenoM))

for (t in 1:length(colphenoM)) {
  trait=colphenoM[t]
  
  print(t)
  
  if (sum(is.na(psM[,t+6]))<=14) {
    
    ## join geno pheno 
    dM=left_join(psM[,c(1:6,t+6)], geno, by="sample")
    
    dMmodel=dM%>%rename(TRAIT=colphenoM[t])%>%select(TTT, sample, TRAIT, starts_with(c("Dim","HMW","LMW"), ignore.case=F))%>%
      pivot_longer(starts_with(c("HMW","LMW")), names_to = "marker", values_to="allele")%>%
      mutate(marker = as.factor(marker))
    
    
    vec_polymorph=dMmodel%>%group_by(marker)%>%
      summarise(niv=n_distinct(allele, na.rm=T))%>%
      filter(niv>=2)%>%
      select(marker)%>%
      as_vector()
    vec_polymorph=droplevels(vec_polymorph)
    
    #correction de Bonferroni:
    seuil_signif=0.05/length(vec_polymorph)
    
    ## filtre marqueurs polymorphes
    dMmodel=dMmodel%>%filter(marker%in%vec_polymorph)   # dMmodel=dMmodel[dMmodel$marker%in%vec_polymorph,]
    
    ## declarer les variables facteurs fixes ou non fixes
    dMmodel$allele=factor(dMmodel$allele)
    dMmodel$TTT=factor(dMmodel$TTT)
    
    
    ## MODEL LINEAIRE
    res_lm=dMmodel%>%nest(data=-marker)%>%
      mutate(test=map(data, ~ summary(lm(TRAIT ~ Dim1 + Dim2 + Dim3 + TTT + allele, data=.x))), table=map(test, tidy))%>%
      unnest(table)%>%select(marker, term, estimate, std.error, statistic, p.value)%>%filter(term=="allele1" & p.value<=seuil_signif)
    
    
    if (dim(res_lm)[1]!=0) { res_lm=res_lm%>%cbind(trait)} else { 
      res_lm[1,]=rep(NA,7)
      res_lm=cbind(res_lm, trait) }
    
    #rbind sur la boucle
    
    if (t==1) { RESM=res_lm } else {
      RESM=rbind(RESM,res_lm) }
    
  }
  
}

write_tsv(RESM, paste(dirout,"Tendances_asso_Exige_Mons.tab", sep="/"))

########################################################################
## codage

code_markers=t(as.matrix(geno[c(1,2),]))
code_markers=cbind(row.names(code_markers), code_markers)
code_markers=data.frame(code_markers[-1,])
colnames(code_markers)=c("marker", "REF=MAJ code 0", "ALT=MIN code 1, 2, etc")
write_tsv(code_markers, paste(dirout,"code_markers_18_genes_glu_FSOV_Exige.tab", sep="/"), col_names=T)
