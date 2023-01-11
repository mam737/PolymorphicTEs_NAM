#### process_initial_EDTA_bed.R
#### Manisha Munasinghe - Last Updated - 01/11/23
#### Take Raw panEDTA TE Annotations in BED Format
#### Filter out nonTE Annotations
#### And store relevant attributes as unique columns
.libPaths('/path/to/Rlibs/')

library(dplyr)
library(stringr)

# Relevant Session Info When Running Script
#> sessionInfo()
#R version 4.0.4 (2021-02-15)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: CentOS Linux 7 (Core)#

#Matrix products: default
#BLAS:   /panfs/roc/msisoft/R/4.0.4/lib64/R/lib/libRblas.so
#LAPACK: /panfs/roc/msisoft/R/4.0.4/lib64/R/lib/libRlapack.so#

#locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C#

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base#

#other attached packages:
#[1] stringr_1.4.0 dplyr_1.0.9#

#loaded via a namespace (and not attached):
# [1] fansi_1.0.3      assertthat_0.2.1 utf8_1.2.2       crayon_1.5.1
# [5] R6_2.5.1         DBI_1.1.2        lifecycle_1.0.1  magrittr_2.0.3
# [9] pillar_1.7.0     stringi_1.7.6    rlang_1.0.2      cli_3.3.0
#[13] vctrs_0.4.1      generics_0.1.2   ellipsis_0.3.2   tools_4.0.4
#[17] glue_1.6.2       purrr_0.3.4      compiler_4.0.4   pkgconfig_2.0.3
#[21] tidyselect_1.1.2 tibble_3.1.6

## From the attribute column, extract the classification as is
extract_raw_classification <- function(attribute_string) {
  classification <- substring(attribute_string, regexpr("*Classification=", attribute_string),regexpr(";Se", attribute_string)-1)
  classification <- unlist(str_split(classification,pattern='='))[2]
  return(classification)
}

## From the attribute column, extract the classification as is
##     depending on value, reclassify into 10 larger groups to
##     reduce the number of classifications + obtain broader classifications
extract_classification <- function(attribute_string) {
  classification <- substring(attribute_string, regexpr("*Classification=", attribute_string),regexpr(";Se", attribute_string)-1)
  classification <- unlist(str_split(classification,pattern='='))[2]
  if (grepl('LINE',classification)) {
    upd_classification <- 'LINE'
  } else if (grepl('LTR',classification)) {
    upd_classification <- 'LTR' 
  } else if (grepl('Helitron',classification)) {
    upd_classification <- 'Helitron' 
  } else if (grepl('MITE',classification)) {
    upd_classification <- 'MITE' 
  } else if (grepl('DNA',classification)) {
    upd_classification <- 'TIR' 
  } else if (grepl('knob',classification)) {
    upd_classification <- 'knob' 
  } else if (grepl('Simple',classification)) {
    upd_classification <- 'Simple Repeat' 
  } else if (grepl('Cent',classification)) {
    upd_classification <- 'Cent/CentC' 
  } else if (grepl('rDNA',classification)) {
    upd_classification <- 'rDNA/spacer' 
  } else if (grepl('Low',classification)) {
    upd_classification <- 'Low Complexity' 
  } else (
    upd_classification <- "subtelomere/4-12-1"
  )
  return(upd_classification)
}

## From the attribute column, extract the family as is
extract_family <- function(attribute_string) {
  family <- substring(attribute_string, regexpr("*Name=", attribute_string),regexpr(";C", attribute_string)-1)
  family <- unlist(str_split(family,pattern='='))[2]
  return(family)
}

## From the attribute column, extract the method as is
##    to see if the call was made via structure vs homology
extract_method <- function(attribute_string) {
  method <- substring(attribute_string, regexpr("*Method=", attribute_string))
  method <- unlist(str_split(method,pattern='='))[2]
  
  if (nchar(method) > 10) {
    new_method <- unlist(str_split(method,pattern=';'))[1]
  } else {
    new_method <- method
  }
  return(new_method)
}

## From the attribute column, extract the identity as is
extract_identity <- function(attribute_string) {
  identity <- substring(attribute_string, regexpr("*dentity=", attribute_string),regexpr(";M", attribute_string)-1)
  identity <- unlist(str_split(identity,pattern='='))[2]
  return(identity)
}

args <- commandArgs(trailingOnly = TRUE)
NAM_line <- unlist(str_split(args,pattern='/'))[8]

bed_data <- read.delim(args, header=FALSE,stringsAsFactors=F)
colnames(bed_data) <- c('chr','start','end','name','score','strand','source','type','attributes')

## Get Raw Classification/Superfamily String
bed_data$raw_superfamily <- unlist(lapply(bed_data$attributes,extract_raw_classification))

## Extract Raw Classification/Superfamily String
##   Rename Depending on Raw to Reduce Number of Classifications
##   See function for details
bed_data$upd_superfamily <- unlist(lapply(bed_data$attributes,extract_classification))

bed_data$raw_family <- unlist(lapply(bed_data$attributes,extract_family))

## Get Method with which the TE called (structure v homology)
bed_data$method <- unlist(lapply(bed_data$attributes,extract_method))

## Get Identity
bed_data$identity <- unlist(lapply(bed_data$attributes,extract_identity))

bed_data$raw_superfamily[bed_data$raw_superfamily == 'LTR/Gypsy'] <- "LTR/Ty3"

## Remove rows w/ non-TE annotations & Helitrons
remove_classification_types <- c('knob/knob180',"subtelomere/4-12-1","knob/TR-1",
                                 "Simple_repeat", "Cent/CentC", "rDNA/spacer",
                                 "Low_complexity","DNA/Helitron")
retained_bed_data <- bed_data[!(bed_data$raw_superfamily %in% remove_classification_types),]


## For all rows called using homology, remove any where the type is repeat region
homology_reads <- retained_bed_data[retained_bed_data$method=='homology',]
retained_homology_reads <- homology_reads[!(homology_reads$type %in% c('repeat_region')),]

## For all rows called using structure, remove duplicates for LTRs
structure_reads <- retained_bed_data[retained_bed_data$method=='structural',]
retained_structure_reads <- structure_reads[!grepl('lTSD_*',unique(structure_reads$name)),]
retained_structure_reads <- retained_structure_reads[!grepl('lLTR_*',unique(retained_structure_reads$name)),] 
retained_structure_reads <- retained_structure_reads[!grepl('LTRRT_*',unique(retained_structure_reads$name)),]
retained_structure_reads <- retained_structure_reads[!grepl('rLTR_*',unique(retained_structure_reads$name)),]
retained_structure_reads <- retained_structure_reads[!grepl('rTSD_*',unique(retained_structure_reads$name)),]


finalized_list <- rbind(retained_homology_reads,retained_structure_reads)

#Adjust Classification Names
finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'LTR/CRM'] <- "LTR/Gypsy"
finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'LTR/Gypsy'] <- "LTR/Ty3"

finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'MITE/DTH'] <- "DTH/MITE"  
finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'MITE/DTA'] <- "DTA/MITE"
finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'MITE/DTT'] <- "DTT/MITE"
finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'MITE/DTM'] <- "DTM/MITE"
finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'MITE/DTC'] <- "DTC/MITE"

finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'DNA/DTH'] <- "DTH/DNA"
finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'DNA/DTA'] <- "DTA/DNA"
finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'DNA/DTT'] <- "DTT/DNA"
finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'DNA/DTM'] <- "DTM/DNA"
finalized_list$raw_superfamily[finalized_list$raw_superfamily == 'DNA/DTC'] <- "DTC/DNA"

output_screened_fname <- paste('/path/to/screened_TE_annotations/updated_screened_',NAM_line,sep='')

write.table(finalized_list,output_screened_fname,col.names=TRUE,quote=FALSE,sep='\t',row.names=F)
