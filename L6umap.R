library(umap)
library(ggplot2)
#CTRL samples
ctrl <- read.table("L6ctrl_log2cpm.csv", sep=",", header=T, row.names=1)
ctrl.labels <- sapply(strsplit(names(ctrl), '\\.'), function(x)x[2])
cell_line <- as.factor(ctrl.labels)
ctrl.umap <- umap(t(ctrl), random_state=111, n_epochs = 1000)
f=data.frame(x=as.numeric(ctrl.umap$layout[,1]),y=as.numeric(ctrl.umap$
                                                               layout[,2]))
#plotting UMAP
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3, aes(colour =
                                                                  cell_line))
pdf("ctrlUMAP.pdf")
print(sp)
dev.off()

#OSI samples
osi <- read.table("L6osi_log2cpm.csv", sep=",", header=T, row.names=1)
osi.labels <- sapply(strsplit(names(osi), '\\.'), function(x)x[2])
cell_line <- as.factor(osi.labels)
osi.umap <- umap(t(osi), random_state=111, n_epochs = 1000)
f=data.frame(x=as.numeric(osi.umap$layout[,1]),y=as.numeric(osi.umap$
                                                              layout[,2]))
#plotting UMAP
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3, aes(colour =
                                                                  cell_line))
pdf("osiUMAP.pdf")
print(sp)
dev.off()