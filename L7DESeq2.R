setwd("C:/Users/eleon/Desktop/data_analysis")

# put the bulkRNA seq file with the counts without normalization (L6) into a variable
cpm <- read.table("L6bulk_pc9_wo_w_osi.txt", header=TRUE, sep = "\t", row.names=1)
# with header=T it considers the first line as a header
# with row.names=1 it means that the 1st column of the table becomes the rowname of the dataframe
# these 2 options are conceptually the same but one for columns and the other for rows

# now check that the file is correct
head(cpm)

#use strsplit to create a large list containing the symbol of the gene and the ensemble ID
symbol_ID <- strsplit(row.names(cpm), ":") #y

#convert the list into a matrix
symbol_ID_matrix <- matrix(unlist(symbol_ID), ncol = 2, byrow = T) #z

# use the paste function to create a large character vector inverting the order: first ensemble ID and then symbol
ID_symbol <- paste(symbol_ID_matrix[,2], symbol_ID_matrix[,1], sep=":")

#put the new row names in the original dataframe (cpm) where we had all the data 
row.names(cpm) <- ID_symbol

# make the 2 matrices with the ctrl data (columns 1,2,3 of the initial dataframe) and acute (78,8,9) or chronic (4,5,6)
acute <-cpm[,c(1,2,3,7,8,9)]
chronic <-cpm[,c(1,2,3,4,5,6)]

#save the variable as .txt file. It is needed to do the row.names and col.names
write.table(acute, file="L7acute_matrix.txt", sep="\t", row.names = TRUE, col.names = TRUE)
write.table(chronic, file="L7chronic_matrix.txt", sep="\t", row.names = TRUE, col.names = TRUE)

# Now we proceed with deSEQ. It wants 3 parameters: 
  # 1. countData (the matrix with the data we want to analyze. We have created it before)
  # 2. colData (our dataframe containing the samples and the conditions)
  # 3. design = the parameter for which we have to perform differential expression. 
  #    Here we put our condition vector with the 2 conditions to be compared for diff expr analysis

#read the new file and set rownames and colnames
acute_df <- read.table("L7acute_matrix.txt",header = TRUE, sep="\t", row.names = 1)
chronic_df <- read.table("L7chronic_matrix.txt",header = TRUE, sep="\t", row.names = 1)

# transform into a matrix
acute_matrix <- as.matrix(acute_df, rownames.force = NA)
chronic_matrix <- as.matrix(chronic_df, rownames.force = NA)

#create the vectors containing the samples and the conditions
samples_acute <- colnames(acute_matrix)
samples_chronic <- colnames(chronic_matrix)
conditions<-c(rep("ctrl",3), rep("treat",3))

#create the dataframe combining the 2 vectors. The names of the vectors will be considered automatically as header
coldata_acute <- data.frame(samples_acute, conditions)
coldata_chronic <- data.frame(samples_chronic, conditions)

#put the rownames taking them from the vectors we created before
rownames(coldata_acute)<-samples_acute
rownames(coldata_chronic)<-samples_chronic

#factor the conditions column
coldata_acute$conditions<-factor(coldata_acute$conditions)
coldata_chronic$conditions<-factor(coldata_chronic$conditions)

# control that effectively it is a factor
class(coldata_acute$conditions)
class(coldata_chronic$conditions)

# now we are ready to do DESeq on acute
library("DESeq2")
dds_acute <- DESeqDataSetFromMatrix(countData = acute_matrix, colData = coldata_acute, design = ~ conditions)
dds_chronic <- DESeqDataSetFromMatrix(countData = chronic_matrix, colData = coldata_chronic, design = ~ conditions)

# Do the study on this output
dds_acute_study <- DESeq(dds_acute)
dds_chronic_study <- DESeq(dds_chronic)

# Save the results
res_acute <- results(dds_acute_study)
res_chronic <- results(dds_chronic_study)

# Save into a file
write.table(res_acute, file="L7DESeq_results_acute.txt", sep="\t", row.names=TRUE,col.names=TRUE)
write.table(res_chronic, file="L7DESeq_results_chronic.txt", sep="\t", row.names=TRUE,col.names=TRUE)

# Filter the data for adjp <0.05. 
res_acute_filter <- res_acute [which(res_acute$padj <= 0.05),]
res_chronic_filter <- res_chronic [which(res_chronic$padj <= 0.05),]
# This command creates a new variable taking all the rows from the res_* dataframe that have a value in the $padj column <0.05. 
# All the columns of the res_* dataframe are selected and maintained. 
# You can see that which is followed by [..., ] the "..." are the rows, the space means all columns.

#save into a file
write.table(res_acute_filter, file="L7DESeq_results_acute_filter.txt", sep="\t", row.names=TRUE,col.names=TRUE)
write.table(res_chronic_filter, file="L7DESeq_results_chronic_filter.txt", sep="\t", row.names=TRUE,col.names=TRUE)

#compute the absolute value of Log2FC column only
res_acute_filter$log2FoldChange <- abs(res_acute_filter$log2FoldChange)
res_chronic_filter$log2FoldChange <- abs(res_chronic_filter$log2FoldChange)

# Filter the data already filtered for p-value also for the log2FC
res_acute_filter <- res_acute_filter [which(res_acute_filter$log2FoldChange >= 1),]
res_chronic_filter <- res_chronic_filter [which(res_chronic_filter$log2FoldChange >= 1),]

# Save the results in a file
write.table(res_acute_filter, file="L7DESeq_results_acute_filter.txt", sep="\t", row.names=TRUE,col.names=TRUE)
write.table(res_chronic_filter, file="L7DESeq_results_chronic_filter.txt", sep="\t", row.names=TRUE,col.names=TRUE)

# put the file osi_log2cpm.csv from lesson 6 in a variable
# it contains data from single cell experiments to intersect with the bulk data
osi <- read.table("L6osi_log2CPM.csv", sep=",", header=T, row.names=1)

#perform the filtering of the single cell data using the DE genes of RNAseq acute treatment
intersection_osi_acute <- intersect(row.names(res_acute_filter), row.names(osi))

# do the same on the chronic
intersection_osi_chronic <- intersect(row.names(res_chronic_filter), row.names(osi))


# select the row of the data frame (osi) of the single cell experiment contaning the DiffExpr 
# genes in common between the single cell data and the DE genes acute or chronic treatment of bulkRNAseq
osi_acute <- osi[which(row.names(osi)%in%intersection_osi_acute), ]
osi_chronic <- osi[which(row.names(osi)%in%intersection_osi_chronic), ]

# save the data in a csv file needed to run the UMAP function
# actually this is not necessary
#write.csv(osi_acute, "L7osi_acute_log2cpm.csv")
#write.csv(osi_chronic, "L7osi_chronic_log2cpm.csv")
# with write.csv you do not have to add the sep parameter

#load the library necessary to run UMAP
library(umap)
library(ggplot2)

#run and plot UMAP on both files
# osi_acute <- read.table("L7osi_acute_log2cpm.csv", sep=",", header=T, row.names=1)
osi.labels_acute <- sapply(strsplit(names(osi_acute), '\\.'), function(x)x[2])
cell_line_acute <- as.factor(osi.labels_acute)
osiUMAP_acute <- umap(t(osi_acute), random_state=111, n_epochs = 1000)
f = data.frame(x = as.numeric(osiUMAP_acute$layout[,1]), y = as.numeric(osiUMAP_acute$
                                                              layout[,2]))
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3)
pdf("L7osi_acute_UMAP.pdf")
print(sp)
dev.off()

# osi_chronic <- read.table("L7osi_chronic_log2cpm.csv", sep=",", header=T, row.names=1)
osi.labels_chronic <- sapply(strsplit(names(osi_chronic), '\\.'), function(x)x[2])
cell_line_chronic <- as.factor(osi.labels_chronic)
osiUMAP_chronic <- umap(t(osi_chronic), random_state=111, n_epochs = 1000)
f = data.frame(x=as.numeric(osiUMAP_chronic$layout[,1]),y=as.numeric(osiUMAP_chronic$
                                                              layout[,2]))

sp_chronic <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3)
pdf("L7osi_chronic_UMAP.pdf")
print(sp_chronic)
dev.off()

# Now we want to filter the scRNAseq dataset with the genes that are differentially expressed both in the acute and in the chronic
# First perform the union on the bulk data combining acute and chronic DE genes.
union_acute_chronic <- union(row.names(res_chronic_filter),row.names(res_acute_filter))

# filter the scRNAseq osi dataset with the diff expr genes from chronic and acute united
osi_union <- osi[which(row.names(osi)%in%union_acute_chronic), ]

# save the dataset into a csv table to run UMAP
# not necessary
#write.csv(osiUMAP5, "L7osiUMAP_union_log2cpm.csv", sep = ",")

#load the library necessary to run UMAP
library(umap)
library(ggplot2)

# run UMAP function 
#osi_union <- read.table("L7osiUMAP5log2cpm.csv", sep=",", header=T, row.names=1)
osiUMAP_union <- umap(t(osi_union), random_state=111, n_epochs = 1000)
f = data.frame(x=as.numeric(osiUMAP_union$layout[,1]), y=as.numeric(osiUMAP_union$layout[,2]))

#plotting UMAP
sp <- ggplot(f, aes(x=x,y=y)) + geom_point(pch=19, cex=0.3)
pdf("L7osiUMAP_union.pdf")
print(sp)
dev.off()