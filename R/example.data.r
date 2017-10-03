#' @name example_annotation
#' @docType data
#' @aliases example_annotation
#' @title Example input annotation data
#' @description For further details, see runMakeGeneList_and_Filter
#' @usage example_annotation
#' @format A data.frame containing gene annotations ordered as a four-column BED format
#' @keywords RPASE example data
NULL

#' @name example_phased_vcf
#' @docType data
#' @aliases example_phased_vcf
#' @title Example input vcf data.
#' @description For further details, see runMakeGeneList_and_Filter
#' @usage example_phased_vcf
#' @format A data.frame with columns ordered according to the VCF standard containing information about phased SNPs
#' @keywords RPASE example data
NULL

#' @name example_phased_block_list
#' @docType data
#' @aliases example_phased_block_list
#' @title Example input of a list of phased block
#' @description Instead of generating a list of phased block using functions runMakeGeneList_and_Filter, one can also alternatively provide a such list with the right formart as shown here.
#' @usage example_phased_block_list
#' @format A list of data.frame of each containing a phased block.In each data.frame, it should contain 6 colunms: 1)chromosome, 2)position, 3)index(of phasing), 4)gene, 5)AD1(counts from one allele), and 6)AD2(counts from another allele).
#' @keywords RPASE example data
NULL
