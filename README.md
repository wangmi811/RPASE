# How to use RPASE
## Introduction

The package RPASE (Read-backed Phasing-based ASE detection) implements an algorithm for identifying genes which show allele-specific expression (ASE) on a per-individual and per-gene basis, using RNA-seq read count data. It takes a VCF containing variants phased using read-backed phasing from [Genome Analysis Toolkit][GATK] and applies an extra-binomial exact test which takes the interdependency of local read counts into consideration.  RPASE leverages information across SNPs and controls for overdispersion in read counts.  The details of the statistical approach implemented by RPASE can be found in the accompanying paper.

Note that RPASE can be applied to any study unit (genes, exons, etc.) as specified by the annotation supplied.

##	Caveat Emptor

* RPASE operates on read counts from phased heterozygous SNPs produced by [GATK][GATK]. It is up to the user to call, phase and filter SNPs using GATK tools. Also, any known artifacts that could result in spurious ASE calls should be removed prior to analysis.

* RPASE provides a [function to estimate overdispersion in read counts](#estimating-overdispersion) in user data.  Alternatively, users may provide their own estimate of overdispersion.  It is up to the user to decide whether the overdispersion estimate provided by `phi.estimate` makes biological sense for their data.

* The p-values returned by RPASE are not adjusted for multiple test correction.

* Finally, large variation in coverage of alleles between SNPs within a phased block (e.g., 5|1000 for one SNP and 1000|10 for another) can result in large amounts of memory being used while calculating the exact p-value.  Such variation is rare in our experience (less than 0.05% of phased blocks) and may occur due to mapping error, repetitive regions, or inaccurate position annotation. among other causes.  The data may be filtered to remove such cases if they occur in the data and memory is limiting.

## Input data              

RPASE requires information on (1) genomic variants phased using read-backed phasing in [GATK][GATK], and (2) an annotation file describing study units.  Genomic variants would typically be provided in a single-sample VCF file, which can be read into R using a variety of packages.  For example, with [vcfR][vcfR]:

```{r, eval=FALSE}
library("vcfR")
the.vcf <- read.vcfR("file.vcf")
phased.variants <- cbind(the.vcf@fix, the.vcf@gt)
```

The annotation information is contained in a four-column dataframe.  This can be read directly from a [BED][BED] file which includes interval names, or could be modified from a [GFF/GTF][GFFGTF] file read using, for example, the [GenomicFeatures][GenomicFeatures] package, which is part of [BioConductor][BioConductor]:

```{r, eval=FALSE}
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("file.gff", format="gff")
annotation <- as.data.frame(genes(txdb))[,c(1,2,3,6)]
colnames(annotation) <- c("chromosome","gene.start","gene.end","gene.name")
```

We include some example variant and annotation dataframes as part of this package.

```{r}
library(RPASE)
data(example_phased_vcf)
example_phased_vcf
data(example_annotation)
example_annotation
```

### Generating a phased block list

The phased block list is a list of dataframes containing information for each study unit.  RPASE provides the `createPhasedBlockList` function for creating a phased block list using variant and annotation information in the formats described above.

```{r}
library(RPASE)
data(example_phased_vcf)
data(example_annotation)
phased_block_list <- createPhasedBlockList(example_phased_vcf, example_annotation)
phased_block_list
```

The `createPhasedBlockList` function provides additional options for filtering phased blocks on various criteria.  The `min.phased.sites` option specified the minimum number of phased SNPs that a phased block must contain.  The `min.coverage` option specifies the minimum number of reads a phased SNP must be covered by to be included in the phased block.  The `heter` option specified the minimum number of reads that must be present for each allele for a SNP to be considered heterozygous.  For defaults, see the function documentation.

The annotation information lists the study units requested by the user, but the variants may not contain sufficient phasing information to cover a complete study unit.  In such cases, a phased block will cover the extent possible.  If the phasing of variants within a study unit contains breaks, then the study unit will be covered by multiple phased blocks.  In `phased_block_list` created above, two phased blocks were created for the single annotation entry named `GeneA`.

Users could also create their own phased block list.  Note that the dataframe for each phased block contains three additional columns describing phase orientation and allele coverage.
The `index` column contains 1 indicating a phase of 0|1, and -1 for a phase of 1|0; within phased blocks, these values may differ from those present in the input variant information so that each phased block has a consistent phase index. The `AD1` and `AD2` columns contain read counts for alleles one and two.

## Estimating overdispersion

Overdispersion in read counts can be estimated directly from a phased block list using the `phi.estimate` function.  This is the function used to provide the default overdispersion estimate for the `runRPASE` function.

```{r, warning=FALSE}
library(RPASE)
data(example_phased_block_list)
phi.estimate(example_phased_block_list)
```

There are only five phased blocks in `example_phased_block_list`, so the overdispersion estimate is quite biased, hence a warning about sample size will be generated if these specific commands are run directly.  We provide a more realistic estimate for this dataset (0.3) in the example below.

## Running RPASE

The `runRPASE` function tests each phased block sequentially, using a default overdispersion estimate calculated from the input phased block list.  The value returned by `runRPASE` is a vector of named p-values, each representing the result of a RPASE test for a single phased block.  The null hypothesis for each test is that there is no allele-specific expression within this phased block.

```{r}
library(RPASE)
data(example_phased_block_list)
runRPASE(example_phased_block_list, phi=0.3, quiet=TRUE)
```

## Running RPASE in parallel

For users who want to use multiple tests in parallel, the `RPASE` function is available, which tests a single phased block.  For `RPASE`, the overdispersion estimate has no default value and must be precomputed using `phi.estimate` or another method.

```{r, eval=FALSE}
data(example_phased_block_list)
my.phi <- 0.3
library(parallel)
my.cluster <- makeCluster(4)
res <- unlist(parLapply(my.cluster, example_phased_block_list, RPASE, phi=my.phi))
stopCluster(my.cluster)
```


[GATK]: https://software.broadinstitute.org/gatk/
[vcfR]: https://cran.r-project.org/package=vcfR
[GenomicFeatures]: https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
[BioConductor]: https://bioconductor.org
[BED]: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
[GFFGTF]: http://www.ensembl.org/info/website/upload/gff.html
