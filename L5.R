#set correct working directory
setwd("C:/Users/eleon/Desktop/data_analysis")
#put the tab delimited file in an R object to be used with R - a variable
x <- read.table("L5dataset.txt", sep = "\t", header = TRUE, row.names = 1)

#x is now a dataframe containing the data from the tab delimited file
#convert it into a matrix
y <- data.matrix(x, rownames.force = NA)

#retrieve the cpm of the new matrix
cpm (y)

#save the changes in a file
write.table(y, file = "L5cpm.txt", append = FALSE, quote = TRUE, sep = "\t", col.names = NA)