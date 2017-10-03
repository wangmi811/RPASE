## ---- eval=FALSE---------------------------------------------------------
#  library("vcfR")
#  the.vcf <- read.vcfR("file.vcf")
#  phased.variants <- cbind(the.vcf@fix, the.vcf@gt)

## ---- eval=FALSE---------------------------------------------------------
#  library(GenomicFeatures)
#  txdb <- makeTxDbFromGFF("file.gff", format="gff")
#  annotation <- as.data.frame(genes(txdb))[,c(1,2,3,6)]
#  colnames(annotation) <- c("chromosome","gene.start","gene.end","gene.name")

## ------------------------------------------------------------------------
library(RPASE)
data(example_phased_vcf)
example_phased_vcf
data(example_annotation)
example_annotation

## ------------------------------------------------------------------------
library(RPASE)
data(example_phased_vcf)
data(example_annotation)
phased_block_list <- createPhasedBlockList(example_phased_vcf, example_annotation)
phased_block_list

## ---- warning=FALSE------------------------------------------------------
library(RPASE)
data(example_phased_block_list)
phi.estimate(example_phased_block_list)

## ------------------------------------------------------------------------
library(RPASE)
data(example_phased_block_list)
runRPASE(example_phased_block_list, phi=0.3, quiet=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  data(example_phased_block_list)
#  my.phi <- 0.3
#  library(parallel)
#  my.cluster <- makeCluster(4)
#  res <- unlist(parLapply(my.cluster, example_phased_block_list, RPASE, phi=my.phi))
#  stopCluster(my.cluster)

