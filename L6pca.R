setwd("C:/Users/eleon/Desktop/data_analysis")
# in this case you don't give the name to the dots because the sampleNames = FALSE -->left image
library(docker4seq)
pca(
  experiment.table = "./L6bulk_pc9_wo_w_osi.txt",
  type = c("counts"),
  covariatesInNames = TRUE,
  samplesName = FALSE,
  principal.components = c(1, 2),
  legend.position = c("bottomleft"),
  pdf = TRUE,
  output.folder = getwd()
)
file.rename("pca.pdf", "L6pca_without_names.pdf")

#in this case you give a name to the dots because the sampleNames=TRUE-->right image
library(docker4seq)
pca(
  experiment.table = "./L6bulk_pc9_wo_w_osi.txt",
  type = c("counts"),
  covariatesInNames = FALSE,
  samplesName = TRUE,
  principal.components = c(1, 2),
  legend.position = c( "bottomleft"),
  pdf = TRUE,
  output.folder = getwd()
)
file.rename("pca.pdf", "L6pca_with_names.pdf")