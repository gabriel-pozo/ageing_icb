#Import metadata
metadata <- read.csv("D:/Scanpy/pseudobulk_ICB/metadata.csv")
rownames(metadata) <- metadata$barcode.1
metadata$barcode.1 <- NULL

#Load barcodes, features and counts matrix
counts <- readMM("D:/Scanpy/pseudobulk_ICB/matrix.mtx")
genes <- read_tsv("D:/Scanpy/pseudobulk_ICB/features.tsv", col_names = FALSE)
cell_ids <- read_tsv("D:/Scanpy/pseudobulk_ICB/barcodes.tsv", col_names = FALSE)

#Create seurat obj
rownames(counts) <- genes$X1
colnames(counts) <- cell_ids$X1
sce <- CreateSeuratObject(counts = counts, meta.data = metadata)

#Load SenMayo geneset
features <- read.delim(file = "D:Scanpy/senescence.txt", header = F, sep = "\t")
features <- features$V1

#Subset matrix to retain only senescence-associated genes
sce <- subset(sce, features = features)
sce$cell_type2 <- "T cells"

#Prepare and create singlecellexperiment obj
sce <- as.SingleCellExperiment(sce3)
sce <- prepSCE(sce, 
               kid = "cell_type2", # subpopulation IDs (e.g., cell types)
               gid = "status",  # condition (e.g., ctrl/AD)
               sid = "patient_id",   # sample IDs (e.g., patient ID)
               drop = TRUE)
               
# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))
sids

# Total number of samples 
ns <- length(sids)
ns

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$sample_id)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reorder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")
ei

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum")
                       
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb)

# Explore the different components of list
str(pb)

# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$sample_id)

# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 


metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id) 

metadata$cluster_id <- factor(metadata$cluster_id)

head(metadata, n = 10)

# Generate vector of cluster IDs
clusters <- levels(metadata$cluster_id)
clusters

clusters[1]

#DEA between T cells from Elderly and Non-elderly samples
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id
cluster_metadata

# Retain all counts
counts <- pb[[clusters[1]]]

cluster_counts <- as.data.frame(as.matrix(counts[, which(colnames(counts) %in% rownames(cluster_metadata))]))

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))    

# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group_id)
                              
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
p1 <- DESeq2::plotPCA(rld, intgroup = "group_id")

# Extract counts matrix and run PCA manually
test <- assay(rld)
pca_data=prcomp(t(test))
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
df_pca_data=data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = colnames(test), condition=rld@colData$group_id)

# Compare the plots
p2 <- ggplot(df_pca_data, aes(PC1,PC2, color = condition))+
  geom_point(size=5)+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))
p1+p2

# Get genes contribution to each PC
pca_weights <- as.data.frame(pca_data$rotation)
write.csv(pca_weights, file='D:Scanpy/artigo_ICB/pca.csv', quote=F, row.names=T)
