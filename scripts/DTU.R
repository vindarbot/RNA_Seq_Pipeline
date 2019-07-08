#BiocManager::install("docopt")

library(docopt)

if (!require("rats")) devtools::install_github("bartongroup/rats", ref="master")
library(rats)

library(tidyverse)
library(wasabi)

"Analyse des gènes différentiellement transcrits (RATs)
Usage:
  DTU.R --cutoff=cutoff --padj=<padj> --gtf=<gtf> --cond_controle=<cond_controle> --cond_traitment=<cond_traitment> --outprefix=<outprefix> 
  Options:
    --cutoff=cutoff                  Log 2 Fold change à appliquer entre les deux contions
    --padj=<padj>                    Seuil de pvalue ajustée à appliquer 
    --gtf=<gtf>                      GTF file
    --cond_controle=<cond_controle>  Condition controle
    --cond_traitment=<cond_traitment> Condition traitement
    --outprefix=<outprefix>          Dossier avec les résultats
" -> doc

opts <- docopt(doc)

cutoff <- opts[['cutoff']]
cutoff <- as.numeric(cutoff)  
padj <- opts[['padj']]
padj <- as.numeric(padj)
outprefix <- opts[['outprefix']]
gtf <- opts[['gtf']]
cond_controle <- opts[['cond_controle']]
cond_traitment <- opts[['cond_traitment']]

ident <- annot2ids(gtf)

# Dossiers parents des quantifications, des echantillons controles et parents
samples_A <- list.files(path='Salmon/quants', pattern=cond_controle, full.names = TRUE)
samples_B <- list.files(path='Salmon/quants', pattern=cond_traitment, full.names = TRUE)

samples_A <- as.vector(samples_A)
samples_B <- as.vector(samples_B)

# 2. Calculate length- & library-normalised abundances. 
#    Scale them to 1M reads for TPM values.
mydata_2 <- fish4rodents(A_paths= samples_A, B_paths= samples_B, 
                         annot=ident, scaleto=100000000)

mydtu_2 <- call_DTU(annot= ident, boot_data_A= mydata_2$boot_data_A, 
                    boot_data_B= mydata_2$boot_data_B, verbose= FALSE,abund_thresh = 5, p_thresh = padj,dprop_thresh = cutoff, threads = 2)


# 4. Plot significance VS effect size:
#plot_overview(mydtu_2)


# 5a. Get all gene and transcript identifiers per category 
# (significant DTU, no DTU, Not Applicable):
myids_2 <- get_dtu_ids(mydtu_2)

#dtu_summary(mydtu_2)
# 5b. Get all gene and transcript identifiers implicated in isoform switching:


DTU_genes <- myids_2$`DTU genes (gene test)`
write_tsv(as.data.frame(DTU_genes),"DTU/DTU.txt")

#dir = setwd("~/Desktop/Data/FLUX-SIMU/AS_genes_detected/RATs/")
write_tsv(mydtu_2$Abundances[[1]], "DTU/abondance_control.txt")
write_tsv(mydtu_2$Abundances[[2]], "DTU/abondance_traitment.txt")


#print(get_plurality_ids(mydtu_2))

#plot_gene(mydtu_2,'AT1G74420')
