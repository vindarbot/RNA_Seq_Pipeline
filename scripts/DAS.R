BiocManager::install("IsoformSwitchAnalyzeR")

library(IsoformSwitchAnalyzeR)
library(tximport)


dir = setwd("~/Desktop/Data")


salmonQuant <- importIsoformExpression(
  parentDir = ("DTE/"),
  addIsofomIdAsColumn = TRUE
)


Design <- data.frame(
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
)

importGTF(
  
)

aSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = Design,
  isoformExonAnnoation = c("Reference/reference_AtRTDv2.gtf"),
  isoformNtFasta = c("Reference/transcriptome.fasta"),
  addAnnotatedORFs=TRUE,
  showProgress = FALSE
)
?importRdata()
# Removal of single isoform genes is the default setting in preFilter() since these genes, per definition, cannot have changes in isoform usage

preFilter(aSwitchList, geneExpressionCutoff = 5)

exampleSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchList,
  reduceToSwitchingGenes=TRUE
)

exampleSwitchListAnalyzed <- analyzeAlternativeSplicing(exampleSwitchListAnalyzed, quiet=TRUE)

consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')

exampleSwitchListAnalyzed <- analyzeSwitchConsequences(
  exampleSwitchListAnalyzed,
  consequencesToAnalyze = consequencesOfInterest, 
  dIFcutoff = 0.4, # very high cutoff for fast runtimes
  showProgress=FALSE
)

table(exampleSwitchListAnalyzed$AlternativeSplicingAnalysis$IR)



analyzeORF(
  exampleSwitchListAnalyzed,
  genomeObject = "Reference/reference_AtRTDv2",
  minORFlength=100,
  orfMethod = "longest",
  cds = NULL,
  PTCDistance = 50,
  startCodons="ATG",
  stopCodons=c("TAA", "TAG", "TGA"),
  showProgress=TRUE,
  quiet=FALSE
)

exampleSwitchListAnalyzed$isoformRepExpression

