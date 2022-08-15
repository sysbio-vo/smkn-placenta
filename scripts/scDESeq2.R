# Adapted from https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

gc()

library(tidyverse)
library(Seurat)
library(DESeq2)
library(SingleCellExperiment)
library(scater)
library(Matrix.utils)
library(magrittr)
library(pheatmap)

data <- readRDS("../data/SP082_SeuratIntegrated_DoubletsAmbientRNAcorrected_filtered_20220808.rds")
original <- data[[]]
sample_ids <- levels(as.factor(original$sample_id))

# Reading data
dirs <- list.dirs("../data/Outs_mtxs", recursive = FALSE)

samples <- c()

for (dir in dirs) {
  sample <- Read10X(data.dir = dir)
  sample <- CreateSeuratObject(counts = sample, project = "SP082")
  sample@meta.data$sampleID <- substr(sapply(strsplit(dir, "/"), `[`, 4), 1, 3)
  samples <- c(sample, samples)
}

merged <- merge(samples[[1]], y = samples[2:length(samples)],
                project = "SP082")


# Matching cell IDs with original data object
merged@meta.data$cellID <- paste(sapply(strsplit(rownames(merged@meta.data),"-"), `[`, 1), merged@meta.data$sampleID, sep="_")
original$sample_id <- gsub(" ", "", original$sample_id)
original$cellID <- paste(sapply(strsplit(rownames(original),"-"), `[`, 1), original$sample_id, sep="_")

# Filtering out cells
merged.filtered <- subset(merged, subset = cellID %in% original$cellID)

# Adding cell type with matching indices
merged.filtered@meta.data$cellType <- original$cell_typev2[match(merged.filtered@meta.data$cellID, original$cellID)]
merged.filtered@meta.data$condition <- original$condition[match(merged.filtered@meta.data$cellID, original$cellID)]


# DE analysis

counts <- merged.filtered@assays$RNA@counts 
metadata <- merged.filtered@meta.data
metadata$cellType <- factor(metadata$cellType)
metadata$sampleID <- factor(metadata$sampleID)
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
groups <- colData(sce)[, c("cellType", "sampleID")]

cells <- purrr::set_names(levels(sce$cellType))
nc <- length(cells)
samples <- purrr::set_names(levels(sce$sampleID))
ns <- length(samples)


# Generate sample level metadata

## Determine the number of cells per sample
table(sce$sampleID)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sampleID))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(samples, sce$sampleID)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], n_cells, row.names = NULL) %>% select(-"cellType")

ei

# Perform QC if not already performed
dim(sce)

# Calculate quality control (QC) metrics
sce <- addPerCellQCMetrics(sce)

# Get cells w/ few/many detected genes
sce$is_outlier <- isOutlier(metric = sce$total, nmads = 2, type = "both", log = TRUE)

# Remove outlier cells
sce <- sce[, !sce$is_outlier]
dim(sce)

## Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

dim(sce)

# Aggregate across cluster-sample groups
groups <- colData(sce)[, c("cellType", "sampleID")]
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)

dim(pb)

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
table(sce$cellType, sce$sampleID)

# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(cells), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(cells), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(cells), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cellType = de_cluster_ids,
                    sampleID = de_samples)

gg_df <- left_join(gg_df, ei[, c("sampleID", "cellID", "condition")]) 


metadata <- gg_df %>%
  dplyr::select(cellType, sampleID, condition) 

metadata 

# Generate vector of cluster IDs
clusters <- levels(as.factor(metadata$cellType))
clusters

# Subset the metadata to only the one cell type
type = "STB"
cluster_metadata <- metadata[which(metadata$cellType == type), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sampleID
head(cluster_metadata)

# Subset the counts to only the B cells
counts <- pb[[type]]

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))], check.names = F)

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ condition)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA

DESeq2::plotPCA(rld, intgroup = "condition")

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("condition"), drop=F])

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
# Plot dispersion estimates
plotDispEsts(dds)


# Output results of Wald test for contrast for stim vs ctrl
contrast <- c("condition", levels(as.factor(cluster_metadata$condition))[2], levels(as.factor(cluster_metadata$condition))[1])

# resultsNames(dds)
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

resultsNames(dds)

res <- lfcShrink(dds, 
                 coef =  2,
                 res=res)

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl

# Write all results to file
write.table(res_tbl,
            paste0("../results/", type, "_", "all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE,
            sep = "\t")