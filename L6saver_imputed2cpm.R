setwd("C:/Users/eleon/Desktop/data_analysis")

#SAVER imputation from raw data (already done by the professor -> commented)
#install.packages("SAVER")
#library(SAVER)

#INPUT = "RNA-5c.csv"
#SEPARATOR = ","
#THREADS = 1 # the number of cores dedicated to run such analysis
#setwd(WDraw.data <- read.table(INPUT, sep=SEPARATOR,
#                       header = T, row.names=1)
#dataset<- as.matrix(raw.data)
#dataset.saver <- saver(dataset, ncores = THREADS,
#                       estimates.only = TRUE)
#write.table(dataset.saver, paste("saver", INPUT,
#                                 sep="_"), sep=SEPARATOR, col.names=NA)

# If you have 2 files (we have ctrl + osi) you have to copy this code and run it on the other file.


# Start from a SAVER imputed .csv file 
# Create the function counts2cpm
counts2cpm <- function(file, sep = ","){
  tmp <- read.table(file, sep=sep, header=T, row.names=1)
  col.sum <- apply(tmp, 2, sum)
  tmp1 <- t(tmp)/col.sum
  tmp1 <- t(tmp1)
  tmp1 <- tmp1 * 1000000
  write.table(tmp1, "L6cpm.csv", sep=",",
              col.names=NA)
  write.table(log2(tmp1 + 1), "log2cpm.csv",
              sep=",", col.names=NA)
}


# Recall the function and run it on your ctrl file
counts2cpm(file="L6saver_ctrl.csv", sep=",")
file.rename(from="log2cpm.csv", to="L6ctrl_log2cpm.csv")


# Recall again the function and run it on your osi file
counts2cpm(file="L6saver_osi.csv", sep=",")
file.rename(from="log2cpm.csv", to="L6osi_log2cpm.csv")

#now you can run UMAP or tSNE