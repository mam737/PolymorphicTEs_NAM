#### classify_polymorphic_TEAnnotations.R
#### Manisha Munasinghe - Last Updated - 01/11/23
#### Intersect TE Annotations with summarised AnchorWave outputs
#### to classify TE features as either shared, polymorphic, or ambiguous
#### in a pairwise comparison between B73 + NAM genome

.libPaths('/path/to/Rlibs/')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(ggplot2)
library(rlist)
# export PATH=/path/to/bedtools2/bin/:$PATH 
# Using bedtools v2.30.0
library(bedtoolsr) 

# Relevant Session Info
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
#[1] bedtoolsr_2.30.0-4 rlist_0.4.6.2      ggplot2_3.3.5      tidyr_1.2.0
#[5] data.table_1.14.2  stringr_1.4.0      dplyr_1.0.9#

#loaded via a namespace (and not attached):
# [1] magrittr_2.0.3   munsell_0.5.0    tidyselect_1.1.2 colorspace_2.0-3
# [5] R6_2.5.1         rlang_1.0.2      fansi_1.0.3      tools_4.0.4
# [9] grid_4.0.4       gtable_0.3.0     utf8_1.2.2       cli_3.3.0
#[13] DBI_1.1.2        withr_2.5.0      ellipsis_0.3.2   assertthat_0.2.1
#[17] tibble_3.1.6     lifecycle_1.0.1  crayon_1.5.1     purrr_0.3.4
#[21] vctrs_0.4.1      glue_1.6.2       stringi_1.7.6    compiler_4.0.4
#[25] pillar_1.7.0     generics_0.1.2   scales_1.2.0     pkgconfig_2.0.3

# For a given feature, extract which AW Blocks it overlaps
pull_AW_blocks <- function(id_name_val,bed_df,comp) {
	if (comp=='id') {
		AW_Blocks <- bed_df %>% subset(id_name==id_name_val) %>% 
		pull(AW_name) %>% unique()
	}
	if (comp=='asm') {
		AW_Blocks <- bed_df %>% subset(asm_name==id_name_val) %>% 
		pull(AW_name) %>% unique()
	}
	AW_Blocks <- paste(AW_Blocks,collapse=',')
	return(AW_Blocks)
}

# There are technically 5 possibly groupings
# alignable, structural variation in B73, structural variation in NAM, unalignable, or missing data
# ensure column for each of these groupings for consistency across all files
check_missing_cols <- function(final_df,ID_struct_tag,ASM_struct_tag) {	
	if (!('alignable_region' %in% colnames(final_df))) {
		final_df$alignable_region <- 0
	}
	if (!(ID_struct_tag %in% colnames(final_df))) {
		final_df[[ID_struct_tag]] <- 0
	}
	if (!(ASM_struct_tag %in% colnames(final_df))) {
		final_df[[ASM_struct_tag]] <- 0
	}
	if(!('unalignable' %in% colnames(final_df))) {
		final_df$unalignable <- 0
	}
	if(!('Missing_Data' %in% colnames(final_df))) {
		final_df$Missing_Data <- 0
	}
	return(final_df)
}

# Pass in directory containing summarised AnchorWave blocks for
# a pairwise comparison between B73 and another NAM genome
# stored as directory containing a file for each chromosome
args <- commandArgs(trailingOnly = TRUE)
AW_dir_name <- args[1]

# Store all summarised AnchorWave blocks across different chromosomes
# in a single list
summarised_AW <- lapply(list.files(AW_dir_name,full.names=T),fread)

# Extract which ID (B73) and ASM (NAM) genotype we're comparing
dir_head <-  paste(unlist(str_split(AW_dir_name,pattern='/'))[1:6],collapse='/')
comp_pattern <- paste(unlist(str_split(unlist(str_split(AW_dir_name,pattern='/'))[8],pattern='_'))[1:2],collapse='_')
id_lineage <- unlist(str_split(unlist(str_split(AW_dir_name,pattern='/'))[8],pattern='_'))[1]
asm_lineage <- unlist(str_split(unlist(str_split(AW_dir_name,pattern='/'))[8],pattern='_'))[2]

# Make Structural Insertion Tags
id_struct_tag <- paste('structural_insertion_in',id_lineage,sep='')
asm_struct_tag <- paste('structural_insertion_in',asm_lineage,sep='')

# Pull ID (B73) Filtered TE Annotations
id_te_dir_name <- paste(dir_head,'filtered_TE_anno',id_lineage,sep='/')
id_te <- lapply(list.files(id_te_dir_name,full.names=T),fread)

# Pull ASM (NAM) Filtered TE Annotations
asm_te_dir_name <- paste(dir_head,'filtered_TE_anno',asm_lineage,sep='/')
asm_te <- lapply(list.files(asm_te_dir_name,full.names=T),fread)

id_te_list <- list()
asm_te_list <- list()

# Go Chromosome by Chromosome through summarised AnchorWave + TE Annotations
for (i in 1:length(summarised_AW)) {

	### Extract all relevant information chromosome by chromosome
	# Extract summarised AW for a specific chromosome
	chr_summarised_AW <- summarised_AW[[i]]

	# Format summarised AW for bedtools
	# Extract ID (B73) coordinates for summarised AW
	chr_id_AW <- chr_summarised_AW %>% 
	select(contains(id_lineage),AW_Block,Reformat_Type)
	colnames(chr_id_AW) <- c('chrom','start','end','name','attributes')
	chr_id_AW$chrom <- paste('chr',chr_id_AW$chrom,sep='')

	# Extract asm coordinates for summarised AW
	chr_asm_AW <- chr_summarised_AW %>% 
	select(contains(asm_lineage),AW_Block,Reformat_Type)
	colnames(chr_asm_AW) <- c('chrom','start','end','name','attributes')
	chr_asm_AW$chrom <- paste('chr',chr_asm_AW$chrom,sep='')

	#Extract id TE bed
	chr_id_te <- as.data.table(id_te[[i]])
	chr_id_te <- chr_id_te %>% 
	select(chr,start,end,name,type,raw_superfamily,upd_superfamily,raw_family,method,class,condense_superfamily)
	chr_id <- unique(chr_id_te$chr)

	#Extract asm TE bed
	chr_asm_te <- as.data.table(asm_te[[i]])
	chr_asm_te <- chr_asm_te %>% 
	select(chr,start,end,name,type,raw_superfamily,upd_superfamily,raw_family,method,class,condense_superfamily)
	chr_asm <- unique(chr_asm_te$chr)

	### bedtools intersect on ID te bed
	id_te_AW_overlap <- bedtoolsr::bt.intersect(chr_id_AW,chr_id_te,wo=T)
	colnames(id_te_AW_overlap) <- c('AW_chr','AW_start','AW_end','AW_name','AW_attributes','id_chr','id_start','id_end','id_name','id_type', 'id_raw_superfamily','id_upd_superfamily','id_raw_family','id_method','id_class','id_condense_superfamily','bp_overlap')
	
	#Add columns for TE_size + prop_overlap
	id_te_AW_overlap$te_size <- id_te_AW_overlap$id_end - id_te_AW_overlap$id_start 
	id_te_AW_overlap$prop_overlap  <- id_te_AW_overlap$bp_overlap/id_te_AW_overlap$te_size
	
	#Pivot data format wider
	id_te_final <- id_te_AW_overlap %>% group_by(id_name,id_chr,id_start,id_end,id_type,id_raw_superfamily,id_upd_superfamily,id_raw_family,id_method,id_class,id_condense_superfamily,AW_attributes) %>% summarise(total_prop_overlap=sum(prop_overlap)) %>% tidyr::pivot_wider(.,id_cols=c('id_name','id_chr','id_start','id_end','id_method','id_type','id_class','id_raw_superfamily','id_upd_superfamily','id_condense_superfamily','id_raw_family',),names_from='AW_attributes',values_from='total_prop_overlap',values_fill=0) %>% arrange(id_start) %>% data.frame()

	#Add additional information on each gene lost during the transformation
	id_te_final$AW_Blocks <- unlist(lapply(id_te_final$id_name,pull_AW_blocks,bed_df=id_te_AW_overlap,comp='id'))

	#Add check in case one of the attributes was not present in this chromosome
	id_te_final <- check_missing_cols(id_te_final,id_struct_tag,asm_struct_tag)

	#Reorder
	id_te_final <- id_te_final[,c('id_name','id_chr','id_start','id_end','id_method', 'id_type', 'id_class', 'id_raw_superfamily', 'id_upd_superfamily', 'id_condense_superfamily','id_raw_family' ,'alignable_region',id_struct_tag,asm_struct_tag,'unalignable','Missing_Data','AW_Blocks')]
	colnames(id_te_final) <- c('TE_name','chr','start','end','method', 'type','class', 'raw_superfamily', 'upd_superfamily','condense_superfamily','raw_family' ,'alignable_region',id_struct_tag,asm_struct_tag,'unalignable','Missing_Data','AW_Blocks')

	#Obtain final classification for each gene based off of proportion
	id_te_final <- id_te_final %>% dplyr::mutate(classification=case_when(
		.[[12]] >= 0.95 ~ 'shared', #.[[12]] is column for alignable proportion overlap
		.[[13]] >= 0.95 ~ 'polymorphic', #[[13]] is column for SV in ID (B73) proportion overlap
		TRUE ~ 'ambiguous'
		))

	id_te_final <- id_te_final %>% data.frame() %>% arrange(start)

	#Add to full list across all chromosomes
	id_te_list <- list.append(id_te_list,id_te_final)

	### bedtools intersect on ID te bed
	asm_te_AW_overlap <- bedtoolsr::bt.intersect(chr_asm_AW,chr_asm_te,wo=T)
	colnames(asm_te_AW_overlap) <- c('AW_chr','AW_start','AW_end','AW_name','AW_attributes','asm_chr','asm_start','asm_end','asm_name','asm_type', 'asm_raw_superfamily','asm_upd_superfamily','asm_raw_family','asm_method','asm_class','asm_condense_superfamily','bp_overlap')
	
	#Add columns for TE_size + prop_overlap
	asm_te_AW_overlap$te_size <- asm_te_AW_overlap$asm_end - asm_te_AW_overlap$asm_start 
	asm_te_AW_overlap$prop_overlap  <- asm_te_AW_overlap$bp_overlap/asm_te_AW_overlap$te_size
	
	#Pivot data format wider
	asm_te_final <- asm_te_AW_overlap %>% group_by(asm_name,asm_chr,asm_start,asm_end,asm_type,asm_raw_superfamily,asm_upd_superfamily,asm_raw_family,asm_method,asm_class,asm_condense_superfamily, AW_attributes) %>% summarise(total_prop_overlap=sum(prop_overlap)) %>% tidyr::pivot_wider(.,id_cols=c('asm_name','asm_chr','asm_start','asm_end','asm_method','asm_type', 'asm_class','asm_raw_superfamily','asm_upd_superfamily','asm_condense_superfamily','asm_raw_family'),names_from='AW_attributes',values_from='total_prop_overlap',values_fill=0) %>% arrange(asm_start) %>% data.frame()
	
	#Add additional information lost during the transformation
	#Add additional information on each gene lost during the transformation
	asm_te_final$AW_Blocks <- unlist(lapply(asm_te_final$asm_name,pull_AW_blocks,bed_df=asm_te_AW_overlap,comp='asm'))

	#Add check in case one of the attributes was not present in this chromosome
	asm_te_final <- check_missing_cols(asm_te_final,id_struct_tag,asm_struct_tag)

	#Reorder
	asm_te_final <- asm_te_final[,c('asm_name','asm_chr','asm_start','asm_end','asm_method', 'asm_type', 'asm_class','asm_raw_superfamily', 'asm_upd_superfamily', 'asm_condense_superfamily','asm_raw_family' ,'alignable_region',id_struct_tag,asm_struct_tag,'unalignable','Missing_Data','AW_Blocks')]
	colnames(asm_te_final) <- c('TE_name','chr','start','end','method', 'type','class', 'raw_superfamily', 'upd_superfamily','condense_superfamily','raw_family' ,'alignable_region',id_struct_tag,asm_struct_tag,'unalignable','Missing_Data','AW_Blocks')

	#Obtain final classification for each gene based off of proportion
	asm_te_final <- asm_te_final %>% dplyr::mutate(classification=case_when(
		.[[12]] >= 0.95 ~ 'shared',#.[[12]] is column for alignable proportion overlap
		.[[14]] >= 0.95 ~ 'polymorphic',#[[14]] is column for SV in ASM (NAM) proportion overlap
		TRUE ~ 'ambiguous'
		))

	asm_te_final <- asm_te_final %>% data.frame() %>% arrange(start)

	#Add to full list across all chromosomes
	asm_te_list <- list.append(asm_te_list,asm_te_final)

}

complete_id_te <- do.call('rbind',id_te_list)
id_te_path <- paste('/path/to/polymorphic_TE_calls/B73/TE/',comp_pattern,'_TE_classification_by_n.tsv',sep='')
write.table(complete_id_te,id_te_path,quote=F,col.names=T,row.names=F)

complete_asm_te <- do.call('rbind',asm_te_list)
asm_te_path <- paste('/path/to/polymorphic_TE_calls/NAM/TE/',comp_pattern,'_TE_classification_by_n.tsv',sep='')
write.table(complete_asm_te,asm_te_path,quote=F,col.names=T,row.names=F)

