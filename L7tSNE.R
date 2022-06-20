library(Rtsne)
library(ggplot2)

# We do tSNE on ctrl samples first
tmp <- read.table("L6ctrl_log2CPM.csv", sep=",", header=T, row.names=1)
tmp.labels <- sapply(strsplit(names(tmp), '\\.'), function(x)x[2])
cell_line <- as.factor(tmp.labels)
set.seed(111)
# We do not want any initial step of PCA before doing tSNE 
tsne_out_ctrl <- Rtsne(as.matrix(t(tmp)), pca=FALSE, perplexity=30,
                  theta=0.0)
f=data.frame(x = as.numeric(tsne_out_ctrl$Y[,1]),y = tsne_out_ctrl$Y[,2])
#plotting tSNE
tsneplot <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3)
pdf("L7ctrl_noPCA.pdf")
print(tsneplot)
dev.off()

# Then we do the same for Osimertinib treatment
tmp <- read.table("L6osi_log2CPM.csv", sep=",", header=T, row.names=1)
tmp.labels <- sapply(strsplit(names(tmp), '\\.'), function(x)x[2])
cell_line <- as.factor(tmp.labels)
set.seed(111)

# Again no PCA
tsne_out_osi <- Rtsne(as.matrix(t(tmp)), pca=FALSE, perplexity=30,
                  theta=0.0)
f=data.frame(x = as.numeric(tsne_out_osi$Y[,1]),y = tsne_out_osi$Y[,2])
# plotting tSNE
tsneplot2 <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3)
pdf("L7osi_noPCA.pdf")
print(tsneplot2)
dev.off()