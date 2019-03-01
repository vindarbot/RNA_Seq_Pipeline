library(BiocManager)
BiocManager::install("DRIMSeq")
source("http://bioconductor.org/biocLite.R")
biocLite("devtools")    # only if devtools not yet installed
biocLite("COMBINE-lab/wasabi")
biocLite("pachterlab/sleuth")
library(wasabi)
library(sleuth)
library(biomaRt)
library(dplyr)
library(tximport)
library(DESeq2)
library(DRIMSeq)

source("https://bioconductor.org/biocLite.R")
pkgs <- rownames(installed.packages())
biocLite(pkgs, type="source")

dir = setwd("~/Desktop/Data")

list.files(system.file("extdata", package = "tximportData"))

seqlevels(TxDb.Athaliana.BioMart.plantsmart25)

sfdirs <- file.path("DTE", c("Col_1","Col_2","hon4_1","hon4_2","hon4_3"))

prepare_fish_for_sleuth(sfdirs)

sample_id <- dir(file.path("DTE"))

condition=c(rep("Control",2),rep("hon4",3))

sfdata=data.frame(sample=list.files("DTE"), path=sfdirs, condition=c(rep("Control",2),rep("hon4",3)),  stringsAsFactors = F)



plantMart <- useMart(biomart="plants_mart",host="plants.ensembl.org",dataset = "athaliana_eg_gene")
t2g <- getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id","description","external_gene_name"),mart = plantMart)
ttg <- dplyr::rename(t2g, target_id= ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

so <- sleuth_prep(sfdata, ~condition)

so <- sleuth_fit(so, ~condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')

so <- sleuth_lrt(so, 'reduced', 'full')

results <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)


sign <- results[which(results$qval<0.1),]

write.table(as.data.frame(sign),"~/Desktop/Data/DTE/results_AtRTD2.txt")






TxDb <- makeTxDbFromGFF("reference.gtf", format="gtf", dataSource="TAIR", organism="Arabidopsis thaliana")


write_tsv(as.data.frame(sign$target_id),"~/Desktop/ANALYSE_COLD_HS/acquired.txt")



so_gene <- sleuth_prep(sfdata, ~condition,target_mapping = ttg,aggregation_column = 'ens_gene')

















Col_1 <- read_tsv("DTE/Col_1/quant.sf")
Col_2 <- read_tsv("DTE/Col_2/quant.sf")
hon4_1 <- read_tsv("DTE/hon4_1/quant.sf")
hon4_2 <- read_tsv("DTE/hon4_2/quant.sf")
hon4_3 <- read_tsv("DTE/hon4_3/quant.sf")


### script python pour généré transcript_to_gene dans Desktop/SNAKEMAKE/Salmon/scripts
tx2gene <- read.table("~/Desktop/SNAKEMAKE_DEG/Salmon/transcript_to_gene.tsv")

files <- file.path("DTE", c("Col_1", "Col_2", "hon4_1","hon4_2","hon4_3"), "quant.sf")
names(files) <- c("Col_1", "Col_2", "hon4_1","hon4_2","hon4_3")

txi <- tximport(files,type="salmon",tx2gene = tx2gene, txOut = TRUE)





sampleTable <- data.frame(condition = factor(c(rep("Col",2), rep("hon4",3))))
rownames(sampleTable) <- colnames(txi$counts)



dds <- DESeqDataSetFromTximport(txi,
                                sampleTable,
                                design= ~ condition)

dds <- DESeq(dds)

rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup="condition",ntop = 1000)
res <- DESeq2::results(dds)

genes_up = res[ which(res$padj < 0.10 & res$log2FoldChange > 1), ]

genes_up <- genes_up[order(genes_up$padj, decreasing = F),]

## Gènes réprimés lorsque le gène testé est muté
genes_down = res[ which(res$padj < 0.10 & res$log2FoldChange < -1), ]

genes_down <- genes_down[order(genes_down$padj, decreasing = F),]

genes_down <- as.data.frame(genes_down)
genes_up <- as.data.frame(genes_up)

identifiants_genes_up <- rownames(genes_up)
identifiants_genes_down <- rownames(genes_down)

identifiants_genes_up <- gsub("[._][A-Za-z0-9]+", "",identifiants_genes_up)

identifiants_genes_down <- gsub("[._][A-Za-z0-9]+", "",identifiants_genes_down)

write.table(as.data.frame(unique(identifiants_genes_up),"~/Desktop/Data/DTE/transcripts_up.txt"))




identifiants <- c(identifiants_genes_down,identifiants_genes_up)

genes_up$ID <- identifiants_genes_up
genes_down$ID <- identifiants_genes_down

## Ensemble des gènes différentiellement exprimés.
genes_signif = rbind(genes_up, genes_down)






gene_description <- read.delim("~/Desktop/RNA_SEQ/gene_description_20131231.txt",
                               header=FALSE, quote="")

genenames = gsub("[.][1234567890]", "",
                 gene_description[,1])

gene_description[,1]=genenames

## Ici on récupère la description des gènes différentiellement exprimés
genes_match_rows <- match(rownames(genes_signif), gene_description[,1])

genes_match_up <- match(rownames(genes_up), gene_description[,1])

Code_to_description <- gene_description[genes_match_rows,c(1,3)]

Code_to_description_up <- gene_description[genes_match_up,c(1,3)]

colnames(Code_to_description) <- c("Code","Description")

colnames(Code_to_description_up) <- c("Code","Description")

genes_match_down <- match(rownames(genes_down), gene_description[,1])

Code_to_description_down <- gene_description[genes_match_down,c(1,3)]

colnames(Code_to_description_down) <- c("Code","Description")


genes_up$Description <- Code_to_description_up$Description

genes_down$Description <- Code_to_description_down$Description


# Pour récupérer les ID symbols à partir des TAIR ID.

ensembl <- useMart(biomart="plants_mart",host="plants.ensembl.org",dataset = "athaliana_eg_gene")
TAIRID <- genes_down$ID
symbol <- getBM(attributes=c("tair_symbol","tair_locus"), mart=ensembl)

genes_up <- (merge(x=symbol,y=genes_up,by.x="tair_locus",by.y="ID"))
genes_up <- genes_up[c(1,2,4,8,9)]
genes_up <- genes_up[order(genes_up$padj,decreasing = F),]

genes_down <- (merge(x=symbol,y=genes_down,by.x="tair_locus",by.y="ID"))
genes_down <- genes_down[c(1,2,4,8,9)]
genes_down <- genes_down[order(genes_down$padj,decreasing = F),]


write_tsv(as.data.frame(genes_up),"~/Desktop/salmon_up.txt")
write_tsv(as.data.frame(genes_down),"~/Desktop/salmon_down.txt")
write_tsv(as.data.frame(identifiants),"~/Desktop/salmonid.txt")
