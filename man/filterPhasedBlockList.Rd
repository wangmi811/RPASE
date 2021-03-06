% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RPASE.r
\name{filterPhasedBlockList}
\alias{filterPhasedBlockList}
\title{Filter an RPASE phased block list}
\usage{
filterPhasedBlockList(unfiltered_phased_block_list, min_phased_sites,
  min_coverage, heter, phased_vcf_data)
}
\arguments{
\item{unfiltered_phased_block_list}{A list of phased blocks from createUnfilteredPhasedBlockList}

\item{min_phased_sites}{minimum number of phased SNPs to include within a phased block}

\item{min_coverage}{minimum coverage for each phased SNP included within a phased block}

\item{heter}{minimum minor-allele coverage for a phased SNP to be included within a phased block}

\item{phased_vcf_data}{A data.frame with columns ordered according to the VCF standard containing information about phased SNPs}
}
\value{
A filtered RPASE phased block list
}
\description{
Filter an RPASE phased block list
}
\details{
For further details, see  \code{\link{createPhasedBlockList}}
}
