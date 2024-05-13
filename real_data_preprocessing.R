library("DESeq2")
library("M3C")
library("SummarizedExperiment")
library("tidyr")
library(reshape2)
library("KEGGgraph")
library(KEGGREST)
library(biomaRt)
library(graph) #for graphNEL 

load("./dat_gene.RData")
load("./dat_mut.RData")

# remove patients who died less than a month after diagnosis, and non-females
dat_gene = dat_gene[,-which(dat_gene$days_to_death <= 30)] # removes 4 
dat_gene = dat_gene[,which(dat_gene$gender == "female") ] # removes 12 male and 1 NA
dat_gene = dat_gene[,-which(dat_gene$days_to_last_follow_up <= 0)]
#have 148 events

#only use subtype basal 
dat_gene = dat_gene[,which(dat_gene$paper_BRCA_Subtype_PAM50 == "Basal")]

#remove so that each patient only has one sample
#dat_gene = dat_gene[,!duplicated(dat_gene$bcr_patient_barcode)] # 5 removed
dat_gene = dat_gene[,!duplicated(substr(dat_gene$bcr_patient_barcode, 1, 12))] # 16 removed

#check
length(dat_gene$bcr_patient_barcode) == length(unique(dat_gene$bcr_patient_barcode))


#in mutation data: keep only sample ids that are present in gene expression data
dim(dat_mut[dat_mut$Tumor_Sample_UUID %in% dat_gene$bcr_patient_barcode])


# find patients that are common in both mutation data and gene expression
# 
 subject_mut = dat_mut$Tumor_Sample_Barcode
# 
 for (i in 1:length(subject_mut)) {
   subject_mut[i] = paste( strsplit(dat_mut$Tumor_Sample_Barcode[i], "-")[[1]][1:4], collapse="-")
   
 }
# 
dat_mut$bcr_patient_barcode = subject_mut
dat_mut = dat_mut[which(subject_mut %in% dat_gene$bcr_patient_barcode),]

#remove multiple sample from same patients
#duplicates = which(duplicated(dat_mut[c("Hugo_Symbol", "bcr_patient_barcode")]))
#dat_mut = dat_mut[-duplicates,]
#for future: check which has most information, then keep that one? 
# but in those values I checked I saw no difference, also same tumor barcode. 


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

#cumulative variance filtering. Keeping features that explain 70% of total variance.
cumsum_var <- cumsum(filtered$statistics$var)
cumsum_cutoff <- cumsum_var[length(cumsum_var)] * 0.1
filtered_names <- filtered$statistics$feature[cumsum_var < cumsum_cutoff]



gexp_filtered = filtered_rna1[filtered_names,]

gexp_colnames = colnames(gexp_filtered)
for (i in 1:length(gexp_colnames)) {
  gexp_colnames[i] = paste( strsplit(gexp_colnames[i], "-")[[1]][1:4], collapse="-")
  
}

colnames(gexp_filtered) = gexp_colnames




#graph of proteins for breast cancer, from KEGG Pathway website
kegg_graph <- parseKGML2Graph(file="hsa05224.xml", genesOnly=FALSE)
kegg_df = parseKGML2DataFrame(file="hsa05224.xml")

# graph of proteins --> graph of genes


#pathway_info = keggGet("hsa05224")
#keggGet("N00049")

#kegg_df = parseKGML2DataFrame(file="hsa05224.xml")





#need mapping of KEGG IDs to Entrez Gene IDs, so all data uses Entrez IDs
kegg_gene_ids = nodes(kegg_graph)

# mappings to Entrez Gene IDs for the KEGG gene IDs
kegg_entrez_mappings = keggConv("ncbi-geneid", kegg_gene_ids)
#kegg_entrez_mappings = keggConv("ncbi-proteinid", kegg_gene_ids)

edge_to = keggConv("ncbi-geneid", kegg_df$to)
edge_from = keggConv("ncbi-geneid", kegg_df$from)

edge_to = as.integer( sapply(strsplit(edge_to, "[:]"), `[`, 2) )
edge_from = as.integer( sapply(strsplit(edge_from, "[:]"), `[`, 2) )

entrez_nodes = rep(0, length(kegg_entrez_mappings))

for (k in 1:length(kegg_entrez_mappings)) {
  entrez_nodes[k] = strsplit(kegg_entrez_mappings[k], "[.]")
}


entrez_nodes = as.integer( sapply(strsplit(kegg_entrez_mappings, "[:]"), `[`, 2) )



#frequency filtering of the mutation data
mut_value_counts = table( dat_mut$Hugo_Symbol )
num_patients = ncol(dat_gene)

prev_limit = 0.05 #prevalence limit
dat_mut_filtered = dat_mut[which( dat_mut$Hugo_Symbol %in% names(which(mut_value_counts > prev_limit*num_patients & mut_value_counts < (1-prev_limit)*num_patients)) ),]
#dat_mut_filtered = dat_mut_filtered[which( dat_mut_filtered$Hugo_Symbol %in% names(which(mut_value_counts < (1-prev_limit)*num_patients)) ),]
#above: need to change to one line: if first lines removes then second prev limit will be skewed

#add mutation features w prev limit > 0.02, which are in the graph:
prev_limit_nodes = 0.02
nodes_mut = dat_mut[which(dat_mut$Hugo_Symbol %in% names(which(mut_value_counts > prev_limit_nodes*num_patients & mut_value_counts < prev_limit*num_patients))),]
#nodes_mut = nodes_mut[which( nodes_mut[which(dat_mut$Entrez_Gene_Id %in% entrez_nodes),]
nodes_mut = nodes_mut[which(nodes_mut$Entrez_Gene_Id %in% entrez_nodes),]
dat_mut_filtered = rbind(dat_mut_filtered, nodes_mut)

already_added_nodes = dat_mut_filtered[which(dat_mut_filtered$Entrez_Gene_Id %in% entrez_nodes),]

#prev_limit 0.02 gives 205, 0.03 gives 59

#force pam50 genes to be included in 
pam50_names = c("UBE2T", "BIRC5", "NUF2", "CDC6", "CCNB1", "TYMS", "MYBL2", "CEP55",
  "MELK", "NDC80", "RRM2", "UBE2C", "CENPF", "PTTG1", "EXO1", "ORC6",
  "ANLN", "CCNE1", "CDC20", "MKI67", "KIF2C", "ACTR3B", "MYC", "EGFR",
  "KRT5", "PHGDH", "CDH3", "MIA", "KRT17", "FOXC1", "SFRP1", "KRT14",
  "ESR1", "SLC39A6", "BAG1", "MAPT", "PGR", "CXXC5", "MLPH", "BCL2",
  "MDM2", "NAT1", "FOXA1", "BLVRA", "MMP11", "GPR160", "FGFR4", "GRB7", "TMEM45B", "ERBB2")

#idx <- which(dat_mut$Hugo_Symbol %in% pam50_names)

#TCGA_PAM50 <- dat_mut[idx,]


# PAM50_already_added = which(rownames(TCGA_PAM50$Hugo_Symbol) %in% dat_mut_filtered$Hugo_Symbol)
# if (identical(PAM50_already_added, integer(0))) {
#   dat_mut_filtered = rbind(dat_mut_filtered, TCGA_PAM50)
#   
# } else {
#   dat_mut_filtered = rbind(dat_mut_filtered, TCGA_PAM50[-PAM50_already_added,])
#   
# }





#force the PAM50 genes to be included in the gene expression data
# idx <- which(rowData(dat_gene)$gene_name %in% pam50_names)
#                
# 
# # extract the PAM50 genes of TCGA-BRCA patients
# TCGA_PAM50 <- RNA_count[idx, ]
# 
# PAM50_already_added = which(rownames(TCGA_PAM50) %in% rownames(gexp_filtered))
# if (identical(PAM50_already_added, integer(0))) {
#   gexp_filtered = rbind(gexp_filtered, TCGA_PAM50)
# } else {
#   gexp_filtered = rbind(gexp_filtered, TCGA_PAM50[-PAM50_already_added,])
# }



#after filtering mutation data, need to transform data set from long format to wide:
mut_wide = data.frame(matrix(0, ncol=dim(gexp_filtered)[2], nrow=length(unique(dat_mut_filtered$Hugo_Symbol))))
colnames(mut_wide) = colnames(gexp_filtered)
rownames(mut_wide) = unique(dat_mut_filtered$Hugo_Symbol)

for (mut in rownames(mut_wide)) {
  patients_w_mut = dat_mut_filtered[which(dat_mut_filtered$Hugo_Symbol == mut), "bcr_patient_barcode"]
  mut_wide[mut, unique(patients_w_mut)] = 1
}






# now need to get the graph for MRF prior of the covariates

#

#force the genes present in the graph to be selected as covariates
ensembl_ids = sapply(strsplit(rownames(filtered_rna1), "[.]"), `[`, 1)
ensembl <- useMart("ensembl")

dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)


dat_gene_id_mapping <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                    filters = "ensembl_gene_id",
                    values = ensembl_ids,
                    mart = dataset)

rna_entrezID = rep(NA, nrow(gexp_filtered)) #dat_gene_id_mapping$entrezgene_id
which_corr = which( ensembl_ids %in%  dat_gene_id_mapping$ensembl_gene_id )
rna_entrezID[which( ensembl_ids %in%  dat_gene_id_mapping$ensembl_gene_id )] = dat_gene_id_mapping$entrezgene_id[which( ensembl_ids %in%  dat_gene_id_mapping$ensembl_gene_id )]
rna_ensembl = rep(NA, nrow(gexp_filtered))
rna_ensembl[which_corr] = dat_gene_id_mapping$ensembl_gene_id[which_corr]

right_order = match(ensembl_ids, rna_ensembl)
rna_entrezID = rna_entrezID[right_order]

idx <- which(rna_entrezID %in% entrez_nodes)
initial_gexp = gexp_filtered

nodes_rna <- filtered_rna1[idx, ]

nodes_already_added = which(rownames(nodes_rna) %in% rownames(gexp_filtered))
if (identical(nodes_already_added, integer(0))) {
  gexp_filtered = rbind(gexp_filtered, nodes_rna)
} else {
  gexp_filtered = rbind(gexp_filtered, nodes_rna[-nodes_already_added,])
}

vars_in_graph = c( nodes_already_added, (nrow(initial_gexp)+1):(nrow(initial_gexp) + nrow(nodes_rna)) )



#create mapping between Ensembl IDs of gene expression data to Entrez Gene IDs
ensembl <- useMart("ensembl")

mart = useEnsembl(biomart = "ensembl",
                  dataset = "hsapiens_gene_ensembl")
dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

ensembl_ids = sapply(strsplit(rownames(gexp_filtered), "[.]"), `[`, 1)


# Retrieve Entrez Gene IDs for the Ensembl IDs of gene expression data
entrez_ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                    filters = "ensembl_gene_id",
                    values = ensembl_ids,
                    mart = dataset)#,
                    #uniqueRows = F)

#entrez_ids$ensembl_gene_id = unique(entrez_ids$ensembl_gene_id)
#entrez_ids$entrezgene_id = unique(entrez_ids$entrezgene_id)

gexp_entrezID = rep(NA, nrow(gexp_filtered)) #dat_gene_id_mapping$entrezgene_id
gexp_entrezID[which( ensembl_ids %in%  entrez_ids$ensembl_gene_id )] = entrez_ids$entrezgene_id[which( ensembl_ids %in%  entrez_ids$ensembl_gene_id )]

which_corr = which( ensembl_ids %in%  entrez_ids$ensembl_gene_id )
rna_entrezID[which( ensembl_ids %in%  entrez_ids$ensembl_gene_id )] = entrez_ids$entrezgene_id[which( ensembl_ids %in%  entrez_ids$ensembl_gene_id )]
rna_ensembl = rep(NA, nrow(gexp_filtered))
rna_ensembl[which_corr] = entrez_ids$ensembl_gene_id[which_corr]

right_order = match(ensembl_ids, rna_ensembl)
gexp_entrezID = gexp_entrezID[right_order]



#define mapping from Hugo Symbol to Entrez Gene ID
entrez_mut = data.frame(Hugo_Symbol=rownames(mut_wide), Entrez_Gene_Id= rownames(mut_wide))

for (hugo in rownames(mut_wide)) {
  entrez_mut[entrez_mut$Hugo_Symbol == hugo,]$Entrez_Gene_Id = dat_mut_filtered[dat_mut_filtered$Hugo_Symbol == hugo, "Entrez_Gene_Id"][1]
}



X = rbind(gexp_filtered, mut_wide)

#vars_in_graph shows the col-index of X of what variables are in the KEGG breast cancer graph
vars_in_graph = c(vars_in_graph, which(colnames(X) %in% unique(nodes_mut$Hugo_Symbol)),
                                 which(colnames(X) %in% unique(already_added_nodes$Hugo_Symbol)) ) #add indices for mutation features in graph

sum( duplicated( sapply(strsplit(colnames(X)[vars_in_graph], "[.]"), `[`, 1) ) )

entrez_gexp = entrez_ids[which(!duplicated(entrez_ids$ensembl_gene_id)),]$entrezgene_id


#some checks: how many mutation features where genes also appear in RNA seq
sum(entrez_mut$Entrez_Gene_Id %in% entrez_ids$entrezgene_id)

#number of genes in KEGGgraph of breast cancer, that appears in X
#could have duplicates
sum(entrez_nodes %in% entrez_mut$Entrez_Gene_Id) + sum(entrez_nodes %in% entrez_ids$entrezgene_id)

#no duplicate version
sum(entrez_nodes %in% unique(c(entrez_mut$Entrez_Gene_Id, entrez_ids$entrezgene_id)))
sum(entrez_nodes %in% unique(c(entrez_mut$Entrez_Gene_Id, gexp_entrezID)))

nodes_in_X = c(entrez_nodes[which(entrez_nodes %in% entrez_mut$Entrez_Gene_Id)], entrez_nodes[which(entrez_nodes %in% entrez_ids$entrezgene_id)])


all_features_entrez = c(gexp_entrezID, entrez_mut$Entrez_Gene_Id)

# construct the graph of X
num_features = nrow(X)
G = matrix(0, ncol=num_features, num_features)
index_shift = nrow(gexp_filtered)

for (i in 1:length(kegg_df$from)) {
  # if both to and from appear among the features, we add the edge
  if (edge_to[i] %in% all_features_entrez & edge_from[i] %in% all_features_entrez) {
    index1 = which(all_features_entrez == edge_to[i])
    index2 = which(all_features_entrez == edge_from[i])
    G[index1, index2] = G[index2, index1] = 1
  }
}

# 
# for (i in 1:length(kegg_df$from)) {
#   #if both to and from appear among the features, we add the edge. 
#   if (edge_to[i] %in% entrez_ids$entrezgene_id & edge_from[i] %in% entrez_ids$entrezgene_id) {
#     index1 = which(entrez_ids$entrezgene_id == edge_to[i])
#     index2 = which(entrez_ids$entrezgene_id == edge_from[i])
#     
#     G[index1,index2] = G[index2,index1] = 1
#   }
#   
#   #same but for mutation data
#   if (edge_to[i] %in% entrez_mut$Entrez_Gene_Id & edge_from[i] %in% entrez_mut$Entrez_Gene_Id) {
#     index1 = index_shift + which(entrez_mut$Entrez_Gene_Id == edge_to[i])
#     index2 = index_shift + which(entrez_mut$Entrez_Gene_Id == edge_from[i])
#     
#     
#     G[index1,index2] = G[index2,index1] = 1
#   }
# }
# 

#if we want to link mutation features and gene expression features corresponding to the same gene
which_mut = index_shift + which(entrez_mut$Entrez_Gene_Id %in% gexp_entrezID)

which_gexp = which(gexp_entrezID %in% entrez_mut$Entrez_Gene_Id)


G[cbind(which_mut, which_gexp)] = 1
G[cbind(which_gexp, which_mut)] = 1


#how many of PAM50 genes appear in graph? Must transform gene names to Entrez Gene ID
ens = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("external_gene_name", "entrezgene_id"),
                 filters = "external_gene_name",
                 values = pam50_names,
                 mart = ens)

pam50_entrezID = rep(NA, length(pam50_names)) #dat_gene_id_mapping$entrezgene_id
pam50_entrezID[which( pam50_names %in%  gene_map$external_gene_name )] = gene_map$external_gene_name[which( pam50_names %in%  gene_map$external_gene_name )]

which_corr = which( pam50_names %in%  gene_map$external_gene_name )
pam50_entrezID[which_corr] = gene_map$external_gene_name[which_corr]
rna_gene_name = rep(NA, length(pam50_names))
rna_gene_name[which_corr] = gene_map$external_gene_name[which_corr]

right_order = match(pam50_names, rna_gene_name)
pam50_entrezID = pam50_entrezID[right_order]


sum(entrez_nodes %in% gene_map$entrezgene_id) #only 5

#occurrences <- table(factor(entrez_nodes, levels = unique(entrez_nodes)), factor(all_features_entrez, levels = unique(all_features_entrez)))

pam50_entrezID = rep(NA, length(pam50_names)) #dat_gene_id_mapping$entrezgene_id
which_corr = which( pam50_names %in%  gene_map$external_gene_name )
pam50_entrezID[which_corr] = gene_map$entrezgene_id[which_corr]
rna_ensembl = rep(NA, length(pam50_names))
rna_ensembl[which_corr] = gene_map$external_gene_name[which_corr]

right_order = match(pam50_names, rna_ensembl)
pam50_entrezID = pam50_entrezID[right_order]




counts <- numeric(length(entrez_nodes))

for (i in seq_along(entrez_nodes)) {
  counts[i] <- sum(entrez_nodes[i] == all_features_entrez)
}

# Print the result
print(counts)

X = t(X)


alive_subjects = which(meta$vital_status=="Alive")
dead_subjects = which(meta$vital_status == "Dead")
meta$time = rep(0,nrow(meta))
meta[alive_subjects,]$time = meta[alive_subjects,]$days_to_last_follow_up
meta[dead_subjects,]$time = meta[dead_subjects,]$days_to_death

#
if (any(is.na(meta$time))) {
  which_na = which(is.na(meta$time))
  meta = meta[-which_na,]
  X = X[-which_na,]
  
}

#add two clinical vars to X and G:
clin = data.frame(age = meta$age_at_diagnosis, treatment=as.integer(meta$treatments) )
X = cbind(X, clin)
clin_graph_row = data.frame(age=rep(0,ncol(G)), treatment = rep(0, ncol(G)))

G = rbind(G, t(clin_graph_row))
clin_graph_col = data.frame(age=rep(0,nrow(G)), treatment = rep(0, nrow(G)))
X = as.matrix(X)

G = cbind(G, clin_graph_col)

G["age", "treatment"] = G["treatment", "age"] = 2

#save(X, meta, file="./RealData/dataset10p.RData")
#save(G, file="./RealData/graph10p.RData")

#save(X, meta, file="./RealData/dataset10p_basal.RData")
#save(G, file="./RealData/graph10p_basal.RData")

#save(G, file="./RealData/graph10p_basal_empty.RData")

#save(X, meta, file="./RealData/dataset_basal.RData")
#save(G, file="./RealData/graph_basal_kegg.RData")

#save(G, file="./RealData/graph_basal_empty.RData")

#save(G, file="./RealData/graph_basal_full.RData" )
G = as.matrix(G)
save(G, file="./RealData/graph_w_clin_full.RData")
load("./RealData/graph_w_clin_full.RData")
save(X, meta, vars_in_graph, file="./RealData/march20long/dataset_w_clin.RData")
#G = matrix(0, nrow=nrow(G), ncol=ncol(G))
#save(G, file="./RealData/graph_w_clin_empty.RData")
