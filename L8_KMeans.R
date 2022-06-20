setwd("C:/Users/eleon/Desktop/data_analysis")

# install the mbkmeans Bioconductor package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mbkmeans")

# apply UMAP on the single cell dataset provided: first load the dataset with the counts in a variable
saver_RNA <- read.table("L8saver_RNA-5c.csv", sep = ",", header = T, row.names = 1)

# Create the function counts2cpm to perform the step from counts to cpm that is required after imputation (see L6)
counts2cpm <- function(file, sep = ","){
  tmp <- read.table(file, sep=sep, header=T, row.names=1)
  col.sum <- apply(tmp, 2, sum)
  tmp1 <- t(tmp)/col.sum
  tmp1 <- t(tmp1)
  tmp1 <- tmp1 * 1000000
  write.table(tmp1, "L8cpm.csv", sep=",",
              col.names=NA)
  write.table(log2(tmp1 + 1), "L8log2cpm.csv",
              sep=",", col.names=NA) 
}

# recall the function and run it on your file
counts2cpm(file="L8saver_RNA-5c.csv", sep=",")
file.rename(from="L8log2cpm.csv", to="L8_RNA_5c_log2cpm.csv")

# load libraries to run UMAP
library(umap)
library(ggplot2)
# load the log2cpm 
log2cpm_RNA <- read.table("L8_RNA_5c_log2cpm.csv", sep=",", header=T, row.names=1)
# create a vector with the names of the cell lines and factor it
celllines_label <- strsplit(colnames(log2cpm_RNA), "\\.")  # . as separation
matrix_cellnames <- matrix(unlist(celllines_label), ncol = 2, byrow = T)
cell_line <- as.factor(matrix_cellnames[,2])
# This can also be done this way:
# celllines_labels <- sapply(strsplit(names(umap_data), '\\.'), function(x)x[2]) # . as separation
# cell_line <- as.factor(umap.labels)
umap_results <- umap(t(log2cpm_RNA), random_state=111, n_epochs = 1000)
f=data.frame(x=as.numeric(umap_results$layout[,1]),y=as.numeric(umap_results$layout[,2]), cell_line = cell_line)
rownames(f) <- colnames(log2cpm_RNA)

# plotting UMAP (aggiorna i files)
sp <- ggplot(f, aes(x=x,y=y, color = cell_line)) + geom_point(pch=19, cex=0.3)
pdf("L8UMAP.pdf")
print(sp)
dev.off()

#run Tsne
library(Rtsne)
library(ggplot2)

# we run it on the same umap_data variable we did umap on before
#umap_data <- read.table("L8_RNA_5c_log2cpm.csv", sep=",", header=T, row.names=1)
#umap.labels <- sapply(strsplit(names(umap_data), '\\.'), function(x)x[2])
#cell_line <- as.factor(umap.labels)

set.seed(111)
# We do not want any initial step of PCA before doing tSNE, so pca = F
tsne_results <- Rtsne(as.matrix(t(umap_data)), pca=FALSE, perplexity=30,theta=0.0)
t=data.frame(x = as.numeric(tsne_results$Y[,1]),y = tsne_results$Y[,2])
# Before plotting tSNE you want to have something inside the variable t
# which is the one you are going to plot to differentiate the cell lines: put cell_lines as a column of t
t$cell_line <- cell_line
t$cell_line <- as.factor(t$cell_line)
rownames(t) <- colnames(log2cpm_RNA)

# plotting tSNE
tsneplot <- ggplot(t, aes(x=x,y=y, color = cell_line)) + geom_point(pch=19, cex=0.3)
pdf("L8tSNE.pdf")
print(tsneplot)
dev.off()

# now  that we have the low dimensional data we can run kmeans
# I want to run kmeans on the UMAP output, the umap_result object, on the column layout
library(mbkmeans)
kmeans_umap <- kmeans(umap_results$layout, 5)
# to plot you want to plot the umap with the k means clustering color code (must be numeric!):
plot_mbkmeans <- ggplot(f, aes(x=x, y=y)) + geom_point(pch = 19, cex = 0.3, color = as.numeric(kmeans_umap$cluster))
pdf("L8mbkmeans_result")
print(plot_mbkmeans)
dev.off()

# do the same with clusterR
library(ClusterR)
kmeans_clusterR <- KMeans_arma(umap_results$layout, clusters = 5, n_iter = 10, seed_mode = "random_subset", verbose = T, CENTROIDS = NULL)
# The next variable is used to color
pr <- predict_KMeans(umap_results$layout, kmeans_clusterR)
plot_clusterR <- ggplot(f, aes(x=x, y=y)) + geom_point(pch = 19, cex = 0.3, color=as.numeric(pr))
pdf("L8kmeans_ClusterR_result")
print(plot_clusterR)
dev.off()
