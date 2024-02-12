# load all libraries used in this tutorial except mlr3
library("TCGAbiolinks")
library("SummarizedExperiment")
#library("DESeq2")
library("dplyr")
library("ggplot2")
library("survival")
#library("survminer")
#library("M3C")
#library("glmnet")
library("plotmo")
library("grpreg")
library("SGL")
library("psbcGroup")
library("psbcSpeedUp")
library("GGally")
library("BhGLM")
library("risksetROC")
library("riskRegression")
library("peperr")
library("c060")
library("rms")
library("survAUC")



# download TCGA breast cancer (BRCA) mRNA-Seq data using GDC api method
query_exp <- TCGAbiolinks::GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq",
  sample.type = c("Primary Tumor")
)
#TCGAbiolinks::GDCdownload(query = query_exp, method = "api")
dat_gene <- TCGAbiolinks::GDCprepare(query = query_exp)

SummarizedExperiment::assays(dat)$unstranded[1:5, 1:2]



##TCGAbiolinks::getProjectSummary("TCGA-BRCA")

#use DESeq2 normalization


meta <- colData(dat_gene)[, c("project_id", "submitter_id", "age_at_diagnosis", "ethnicity", "gender", "days_to_death", "days_to_last_follow_up", "vital_status", "paper_BRCA_Subtype_PAM50", "treatments")]
meta$treatments <- unlist(lapply(meta$treatments, function(xx) {
  any(xx$treatment_or_therapy == "yes")
}))
dds <- DESeq2::DESeqDataSetFromMatrix(assays(dat_gene)$unstranded, colData = meta, design = ~1)
dds2 <- DESeq2::estimateSizeFactors(dds)
RNA_count <- DESeq2::counts(dds2, normalized = TRUE) # save
RNA_count[1:5, 1:2]

save(dat_gene, file="dat_gene.RData")
save(dds2, file="dds2.RData")




query_mut <- TCGAbiolinks::GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  #workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  #experimental.strategy = "RNA-Seq",
  #sample.type = c("Primary Tumor")
)
#TCGAbiolinks::GDCdownload(query = query_mut, method = "api")
dat_mut <- TCGAbiolinks::GDCprepare(query = query_mut)

SummarizedExperiment::assays(dat_mut2)$unstranded[1:5, 1:2]
dat_mut = dat_mut2

save(dat_mut, file="dat_mut2.RData")

#dat_gene$sample_id seems to correspond to dat_mut$case_id
