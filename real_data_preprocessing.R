library("DESeq2")
library("M3C")
library("SummarizedExperiment")
library("tidyr")
library(reshape2)
library("KEGGgraph")
library(KEGGREST)
library(biomaRt)

load("./dat_gene.RData")
load("./dat_mut.RData")



# remove patients who died less than a month after diagnosis, and non-females
dat_gene = dat_gene[,-which(dat_gene$days_to_death <= 30)] # removes 4 
dat_gene = dat_gene[,which(dat_gene$gender == "female") ] # removes 12 male and 1 NA
#have 148 events

#remove so that each patient only has one sample
dat_gene = dat_gene[,!duplicated(dat_gene$bcr_patient_barcode)] # 5 removed

#check
length(dat_gene$bcr_patient_barcode) == length(unique(dat_gene$bcr_patient_barcode))


#in mutation data: keep only sample ids that are present in gene expression data
dim(dat_mut[dat_mut$Tumor_Sample_UUID %in% dat_gene$bcr_patient_barcode])


# find patients that are common in both mutation data and gene expression
subject_mut = dat_mut$Tumor_Sample_Barcode

for (i in 1:length(subject_mut)) {
  subject_mut[i] = paste( strsplit(dat_mut$Tumor_Sample_Barcode[i], "-")[[1]][1:4], collapse="-")
  
}

dat_mut$bcr_patient_barcode = subject_mut
dat_mut = dat_mut[which(subject_mut %in% dat_gene$bcr_patient_barcode),]

#remove multiple sample from same patients
duplicates = which(duplicated(dat_mut[c("Hugo_Symbol", "bcr_patient_barcode")]))
dat_mut = dat_mut[-duplicates,]
#for future: check which has most information, then keep that one? 
# but in those values I checked I saw no difference, also same tumor barcode. 


#frequency filtering of the mutation data
mut_value_counts = table( dat_mut$Hugo_Symbol )
num_patients = length(dat_gene$bcr_patient_barcode)

prev_limit = 0.02 #prevalence limit
dat_mut_filtered = dat_mut[which( dat_mut$Hugo_Symbol %in% names(which(mut_value_counts > prev_limit*num_patients)) ),]
dat_mut_filtered = dat_mut_filtered[which( dat_mut_filtered$Hugo_Symbol %in% names(which(mut_value_counts < (1-prev_limit)*num_patients)) ),]

#prev_limit 0.02 gives 208, 0.03 gives 60



#rna normalisation
meta <- colData(dat_gene)[, c("project_id", "submitter_id", "age_at_diagnosis", "ethnicity", "gender", "days_to_death", "days_to_last_follow_up", "vital_status", "paper_BRCA_Subtype_PAM50", "treatments")]
meta$treatments <- unlist(lapply(meta$treatments, function(xx) {
  any(xx$treatment_or_therapy == "yes")
}))
dds <- DESeq2::DESeqDataSetFromMatrix(assays(dat_gene)$unstranded, colData = meta, design = ~1)
dds2 <- DESeq2::estimateSizeFactors(dds)

RNA_count <- DESeq2::counts(dds2, normalized = TRUE)

#only use protein coded genes
filtered_rna <- RNA_count[rowData(dat_gene)$gene_type == "protein_coding", ]

#variance filtering of gene expression data
RNA_log2count <- log2(filtered_rna + 1)
filtered <- M3C::featurefilter(RNA_log2count, percentile = 100, method = "var", topN = 5)
filtered_rna1 <- filtered$filtered_data

#dat_gene_filtered = dat_gene[dat_gene$barcode %in% filtered_rna1,]

#cumulative variance filtering. Keeping features that explain 50% of total variance.
cumsum_var <- cumsum(filtered$statistics$var)
cumsum_cutoff <- cumsum_var[length(cumsum_var)] * 0.5
filtered_names <- filtered$statistics$feature[cumsum_var < cumsum_cutoff]



gexp_filtered = filtered_rna1[filtered_names,]

gexp_colnames = colnames(gexp_filtered)
for (i in 1:length(gexp_colnames)) {
  gexp_colnames[i] = paste( strsplit(gexp_colnames[i], "-")[[1]][1:4], collapse="-")
  
}

colnames(gexp_filtered) = gexp_colnames


#after filtering mutation data, need to transform data set from long format to wide:
mut_wide = data.frame(matrix(0, ncol=dim(gexp_filtered)[2], nrow=length(unique(dat_mut_filtered$Hugo_Symbol))))
colnames(mut_wide) = colnames(gexp_filtered)
rownames(mut_wide) = unique(dat_mut_filtered$Hugo_Symbol)

for (mut in rownames(mut_wide)) {
  patients_w_mut = dat_mut_filtered[which(dat_mut_filtered$Hugo_Symbol == mut), "bcr_patient_barcode"]
  mut_wide[mut, unique(patients_w_mut)] = 1
}



X = rbind(gexp_filtered, mut_wide)


# now need to get the graph for MRF prior of the covariates

#graph of proteins for breast cancer, from KEGG Pathway website
kegg_graph <- parseKGML2Graph(file="hsa05224.xml", genesOnly=FALSE)
kegg_df = parseKGML2DataFrame(file="hsa05224.xml")

# graph of proteins --> graph of genes


pathway_info = keggGet("hsa05224")



#create mapping between Ensembl IDs of gene expression data to Entrez Gene IDs
ensembl <- useMart("ensembl")

dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

ensembl_ids = first_elements <- sapply(strsplit(rownames(gexp_filtered), "[.]"), `[`, 1)

# Retrieve Entrez Gene IDs for the Ensembl IDs of gene expression data
entrez_ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                    filters = "ensembl_gene_id",
                    values = ensembl_ids,
                    mart = dataset)


#define mapping from Hugo Symbol to Entrez Gene ID
entrez_mut = data.frame(Hugo_Symbol=rownames(mut_wide), Entrez_Gene_Id= rownames(mut_wide))

for (hugo in rownames(mut_wide)) {
  entrez_mut[entrez_mut$Hugo_Symbol == hugo,]$Entrez_Gene_Id = dat_mut_filtered[dat_mut_filtered$Hugo_Symbol == hugo, "Entrez_Gene_Id"][1]
}



#need mapping of KEGG IDs to Entrez Gene IDs, so all data uses Entrez IDs
kegg_gene_ids = nodes(kegg_graph)

# mappings to Entrez Gene IDs for the KEGG gene IDs
kegg_entrez_mappings = keggConv("ncbi-geneid", kegg_gene_ids)

edge_to = keggConv("ncbi-geneid", kegg_df$to)
edge_from = keggConv("ncbi-geneid", kegg_df$from)

edge_to = as.integer( sapply(strsplit(edge_to, "[:]"), `[`, 2) )
edge_from = as.integer( sapply(strsplit(edge_from, "[:]"), `[`, 2) )

entrez_nodes = rep(0, length(kegg_entrez_mappings))

for (k in 1:length(kegg_entrez_mappings)) {
  entrez_nodes[k] = strsplit(kegg_entrez_mappings[k], "[.]")
}


entrez_nodes = as.integer( sapply(strsplit(kegg_entrez_mappings, "[:]"), `[`, 2) )




#some checks: how many mutation features where genes also appear in RNA seq
sum(entrez_mut$Entrez_Gene_Id %in% entrez_ids$entrezgene_id)

#number of genes in KEGGgraph of breast cancer, that appears in X
sum(entrez_nodes %in% entrez_mut$Entrez_Gene_Id) + sum(entrez_nodes %in% entrez_ids$entrezgene_id)

nodes_in_X = c(entrez_nodes[which(entrez_nodes %in% entrez_mut$Entrez_Gene_Id)], entrez_nodes[which(entrez_nodes %in% entrez_ids$entrezgene_id)])



# construct the graph of X
G = matrix(0, ncol=num_features, num_features)
index_shift = dim(gexp_filtered)[1]

for (i in 1:length(kegg_df$from)) {
  #if both to and from appear among the features, we add the edge. 
  if (edge_to[i] %in% entrez_ids$entrezgene_id & edge_from[i] %in% entrez_ids$entrezgene_id) {
    index1 = which(entrez_ids$entrezgene_id == edge_to[i])
    index2 = which(entrez_ids$entrezgene_id == edge_from[i])
    
    G[index1,index2] = G[index2,index1] = 1
  }
  
  #same but for mutation data
  if (edge_to[i] %in% entrez_mut$Entrez_Gene_Id & edge_from[i] %in% entrez_mut$Entrez_Gene_Id) {
    index1 = index_shift + which(entrez_mut$Entrez_Gene_Id == edge_to[i])
    index2 = index_shift + which(entrez_mut$Entrez_Gene_Id == edge_from[i])
    
    
    G[index1,index2] = G[index2,index1] = 1
  }
}


#if we want to link mutation features and gene expression features corresponding to the same gene
which_mut = index_shift + which(entrez_mut$Entrez_Gene_Id %in% entrez_ids$entrezgene_id)

which_gexp = which(entrez_ids$entrezgene_id %in% entrez_mut$Entrez_Gene_Id)

G[cbind(which_mut, which_gexp)] = 1
G[cbind(which_gexp, which_mut)] = 1



#remember split of data set into training and test
