#### classify_polymorphic_geneAnnotations.R
#### Manisha Munasinghe - Last Updated - 01/11/23
#### Intersect gene Annotations with summarised AnchorWave outputs
#### to classify TE features as either shared, polymorphic, or ambiguous
#### in a pairwise comparison between B73 + NAM genome
#### doing both an exon-only + full-length approach

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

source('/path/to/function_scripts/classify_feature_functions.R')

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

# Pull ID (B73) Gene Annotations
# From Exon by Exon Annotation Generate Full Length Annotation
id_gff_dir_name <- paste(dir_head,'gff_gene_info',id_lineage,sep='/')
id_gff_exon <- lapply(list.files(id_gff_dir_name,full.names=T),fread)
id_gff_full <- lapply(id_gff_exon,generate_full_length)

# Pull ASM (NAM) Gene Annotations
# From Exon by Exon Annotation Generate Full Length Annotation
asm_gff_dir_name <- paste(dir_head,'gff_gene_info',asm_lineage,sep='/')
asm_gff_exon <- lapply(list.files(asm_gff_dir_name,full.names=T),fread)
asm_gff_full <- lapply(asm_gff_exon,generate_full_length)

id_exon_list <- list()
id_full_list <- list()

asm_exon_list <- list()
asm_full_list <- list()

# Go Chromosome by Chromosome through summarised AnchorWave + Gene Annotations
for (i in 1:length(summarised_AW)) {

	### Extract all relevant information chromosome by chromosome
	# Extract summarised AW for a specific chromosome
	chr_summarised_AW <- summarised_AW[[i]]

	# Format summarised AW for bedtools
	# Extract id coordinates for summarised AW
	chr_id_AW <- chr_summarised_AW %>% 
	select(contains(id_lineage),AW_Block,Reformat_Type)
	colnames(chr_id_AW) <- c('chrom','start','end','name','attributes')
	chr_id_AW$chrom <- paste('chr',chr_id_AW$chrom,sep='')

	# Extract asm coordinates for summarised AW
	chr_asm_AW <- chr_summarised_AW %>% 
	select(contains(asm_lineage),AW_Block,Reformat_Type)
	colnames(chr_asm_AW) <- c('chrom','start','end','name','attributes')
	chr_asm_AW$chrom <- paste('chr',chr_asm_AW$chrom,sep='')

	# Extract id exon-separated gene gff for chromosome
	chr_id_gff_exon <- as.data.table(id_gff_exon[[i]])
	chr_id <-  unique(chr_id_gff_exon$Chr)
	chr_id_gff_exon$exon_id <- unlist(lapply(chr_id_gff_exon$Info,extract_exon_classification))
	chr_id_gff_exon <- chr_id_gff_exon[,c(1,2,3,6,7)] 
	# Extract id full length gene gff for chromosome
	chr_id_gff_full <- as.data.table(id_gff_full[[i]])

	# Extract asm exon-separated gene gff for chromosome
	chr_asm_gff_exon <- as.data.table(asm_gff_exon[[i]])
	chr_asm <- unique(chr_asm_gff_exon$Chr)
	chr_asm_gff_exon$exon_id <- unlist(lapply(chr_asm_gff_exon$Info,extract_exon_classification))
	chr_asm_gff_exon <- chr_asm_gff_exon[,c(1,2,3,6,7)] 
	# Extract asm full length gene gff for chromosome
	chr_asm_gff_full <- as.data.table(asm_gff_full[[i]])

	### bedtools intersect on ID exon gff
	id_gff_exon_AW_overlap <- bedtoolsr::bt.intersect(chr_id_AW,chr_id_gff_exon,wo=T)
	colnames(id_gff_exon_AW_overlap) <- c('AW_chr','AW_start','AW_end','AW_name','AW_attributes','id_chr','id_start','id_end','id_name','id_exon','bp_overlap')
	
	#Add columns for exon size + prop overlap for each exon
	id_gff_exon_AW_overlap$exon_size <- id_gff_exon_AW_overlap$id_end - id_gff_exon_AW_overlap$id_start 
	id_gff_exon_AW_overlap$prop_overlap  <- id_gff_exon_AW_overlap$bp_overlap/id_gff_exon_AW_overlap$exon_size
	
	#Add columns summarising exon data when grouped by AW attribute (alignable_region,structural_insertion_inX,etc.) and gene name
	#Add info on total exons and proportion overlap
	id_gff_exon_AW_overlap_combined <- id_gff_exon_AW_overlap %>% 
		group_by(AW_attributes,id_name) %>% 
		summarise(tot_exon_bp_overlap=sum(bp_overlap),tot_exon_prop_overalp=sum(prop_overlap)/n(),num_exons=n()) %>% 
		ungroup() %>% group_by(id_name) %>% 
		mutate(tot_exons=sum(num_exons),prop_overlap=tot_exon_bp_overlap/sum(tot_exon_bp_overlap))

	#Pivot data format wider
	id_gff_final <- id_gff_exon_AW_overlap_combined %>% tidyr::pivot_wider(.,id_cols=c("id_name"),names_from='AW_attributes',values_from='prop_overlap',values_fill=0)

	#Add additional information on each gene lost during the transformation
	id_gff_final$AW_Blocks <- unlist(lapply(id_gff_final$id_name,pull_AW_blocks,bed_df=id_gff_exon_AW_overlap,comp='id'))
	id_gff_final$chr <- chr_id
	id_gff_final$start <- unlist(lapply(id_gff_final$id_name,add_gene_start,gff_exon=chr_id_gff_exon))
	id_gff_final$end <- unlist(lapply(id_gff_final$id_name,add_gene_end,gff_exon=chr_id_gff_exon))
	id_gff_final$exon_num <- unlist(lapply(id_gff_final$id_name,add_gene_exon_number,gff_exon=chr_id_gff_exon))

	#Add check in case one of the attributes was not present in this chromosome
	id_gff_final <- check_missing_cols(id_gff_final,id_struct_tag,asm_struct_tag)

	#Reorder
	id_gff_final <- id_gff_final[,c('id_name','chr','start','end','exon_num','alignable_region',id_struct_tag,asm_struct_tag,'unalignable','Missing_Data','AW_Blocks')]
	
	#Obtain final classification for each gene based off of proportion
	id_gff_final <- id_gff_final %>% data.frame() %>% dplyr::mutate(classification=case_when(
		.[[6]] >= 0.95 ~ 'shared',#.[[6]] is column for alignable proportion overlap in exon df
		.[[7]] >= 0.95 ~ 'polymorphic',#[[7]] is column for SV in ID (B73) proportion overlap in exon df
		TRUE ~ 'ambiguous'
		))

	id_gff_final <- id_gff_final %>% data.frame() %>% arrange(start)

	#Add to full list across all chromosomes
	id_exon_list <- list.append(id_exon_list,id_gff_final)

	### bedtools intersect on ID full gff
	id_gff_full_AW_overlap <- bedtoolsr::bt.intersect(chr_id_AW,chr_id_gff_full,wo=T)
	colnames(id_gff_full_AW_overlap) <- c('AW_chr','AW_start','AW_end','AW_name','AW_attributes','id_chr','id_start','id_end','id_name','bp_overlap')
	
	#Add columns for gene size + prop overlap for each exon
	id_gff_full_AW_overlap$gene_size <- id_gff_full_AW_overlap$id_end - id_gff_full_AW_overlap$id_start 
	id_gff_full_AW_overlap$prop_overlap  <- id_gff_full_AW_overlap$bp_overlap/id_gff_full_AW_overlap$gene_size
	
	#Add columns summarising gene data when grouped by AW attribute (alignable_region,structural_insertion_inX,etc.) and gene name
	#Combine info for a gene if multiple of same attribute so one total value/attribute
	#Add info on total bp + prop overlap for each attribute
	id_gff_full_AW_overlap_combined <- id_gff_full_AW_overlap %>% group_by(id_name,AW_attributes) %>%
		summarise(tot_bp_overlap=sum(bp_overlap),tot_prop_overlap=sum(prop_overlap))

	#Pivot data format wider
	id_gff_full_final <- id_gff_full_AW_overlap_combined %>% select(AW_attributes,id_name,tot_prop_overlap) %>% tidyr::pivot_wider(.,id_cols=c("id_name"),names_from='AW_attributes',values_from='tot_prop_overlap',values_fill=0)

	#Add additional information on each gene lost during the transformation
	id_gff_full_final$AW_Blocks <- unlist(lapply(id_gff_full_final$id_name,pull_AW_blocks,bed_df=id_gff_full_AW_overlap,comp='id'))
	id_gff_full_final$chr <- chr_id
	id_gff_full_final$start <- unlist(lapply(id_gff_full_final$id_name,add_gene_start,gff_exon=chr_id_gff_exon))
	id_gff_full_final$end <- unlist(lapply(id_gff_full_final$id_name,add_gene_end,gff_exon=chr_id_gff_exon))
	
	#Add check in case one of the attributes was not present in this chromosome
	id_gff_full_final <- check_missing_cols(id_gff_full_final,id_struct_tag,asm_struct_tag)

	#Reorder
	id_gff_full_final <- id_gff_full_final[,c('id_name','chr','start','end','alignable_region',id_struct_tag,asm_struct_tag,'unalignable','Missing_Data',"AW_Blocks")]
	
	#Obtain final classification for each gene based off of proportion
	id_gff_full_final <- id_gff_full_final %>% data.frame() %>% dplyr::mutate(classification=case_when(
		.[[5]] >= 0.95 ~ 'shared',#.[[5]] is column for alignable proportion overlap in full df
		.[[6]] >= 0.95 ~ 'polymorphic',#[[6]] is column for SV in ID (B73) proportion overlap in full df
		TRUE ~ 'ambiguous'
		))	

	id_gff_full_final <- id_gff_full_final %>% data.frame() %>% arrange(start)

	#Add to full list across all chromosomes
	id_full_list <- list.append(id_full_list,id_gff_full_final)

	### bedtools intersect on asm exon gff
	asm_gff_exon_AW_overlap <- bedtoolsr::bt.intersect(chr_asm_AW,chr_asm_gff_exon,wo=T)
	colnames(asm_gff_exon_AW_overlap) <- c('AW_chr','AW_start','AW_end','AW_name','AW_attributes','asm_chr','asm_start','asm_end','asm_name','asm_exon','bp_overlap')

	#Add columns for exon size + prop overlap for each exon
	asm_gff_exon_AW_overlap$exon_size <- asm_gff_exon_AW_overlap$asm_end - asm_gff_exon_AW_overlap$asm_start 
	asm_gff_exon_AW_overlap$prop_overlap  <- asm_gff_exon_AW_overlap$bp_overlap/asm_gff_exon_AW_overlap$exon_size

	#Add columns summarising exon data when grouped by AW attribute (alignable_region,structural_insertion_inX,etc.) and gene name
	#Add info on total exons and proportion overlap
	asm_gff_exon_AW_overlap_combined <- asm_gff_exon_AW_overlap %>% 
		group_by(AW_attributes,asm_name) %>% 
		summarise(tot_exon_bp_overlap=sum(bp_overlap),tot_exon_prop_overalp=sum(prop_overlap)/n(),num_exons=n()) %>% 
		ungroup() %>% group_by(asm_name) %>% 
		mutate(tot_exons=sum(num_exons),prop_overlap=tot_exon_bp_overlap/sum(tot_exon_bp_overlap))

	#Pivot data format wider
	asm_gff_final <- asm_gff_exon_AW_overlap_combined %>% tidyr::pivot_wider(.,id_cols=c("asm_name"),names_from='AW_attributes',values_from='prop_overlap',values_fill=0)

	#Add additional information on each gene lost during the transformation
	asm_gff_final$AW_Blocks <- unlist(lapply(asm_gff_final$asm_name,pull_AW_blocks,bed_df=asm_gff_exon_AW_overlap,comp='asm'))
	asm_gff_final$chr <- chr_asm
	asm_gff_final$start <- unlist(lapply(asm_gff_final$asm_name,add_gene_start,gff_exon=chr_asm_gff_exon))
	asm_gff_final$end <- unlist(lapply(asm_gff_final$asm_name,add_gene_end,gff_exon=chr_asm_gff_exon))
	asm_gff_final$exon_num <- unlist(lapply(asm_gff_final$asm_name,add_gene_exon_number,gff_exon=chr_asm_gff_exon))

	#Add check in case one of the attributes was not present in this chromosome
	asm_gff_final <- check_missing_cols(asm_gff_final,asm_struct_tag,asm_struct_tag)

	#Reorder
	asm_gff_final <- asm_gff_final[,c('asm_name','chr','start','end','exon_num','alignable_region',id_struct_tag,asm_struct_tag,'unalignable','Missing_Data','AW_Blocks')]
	
	#Obtain final classification for each gene based off of proportion
	asm_gff_final <- asm_gff_final %>% data.frame() %>% dplyr::mutate(classification=case_when(
		.[[6]] >= 0.95 ~ 'shared',#.[[6]] is column for alignable proportion overlap in exon df
		.[[8]] >= 0.95 ~ 'polymorphic',#[[8]] is column for SV in ASM (NAM) proportion overlap in exon df
		TRUE ~ 'ambiguous'
		))

	asm_gff_final <- asm_gff_final %>% data.frame() %>% arrange(start,end)

	#Add to full list across all chromosomes
	asm_exon_list <- list.append(asm_exon_list,asm_gff_final)

	### bedtools intersect on asm full gff
	asm_gff_full_AW_overlap <- bedtoolsr::bt.intersect(chr_asm_AW,chr_asm_gff_full,wo=T)
	colnames(asm_gff_full_AW_overlap) <- c('AW_chr','AW_start','AW_end','AW_name','AW_attributes','asm_chr','asm_start','asm_end','asm_name','bp_overlap')

	#Add columns for gene size + prop overlap for each exon
	asm_gff_full_AW_overlap$gene_size <- asm_gff_full_AW_overlap$asm_end - asm_gff_full_AW_overlap$asm_start 
	asm_gff_full_AW_overlap$prop_overlap  <- asm_gff_full_AW_overlap$bp_overlap/asm_gff_full_AW_overlap$gene_size

	#Add columns summarising gene data when grouped by AW attribute (alignable_region,structural_insertion_inX,etc.) and gene name
	#Combine info for a gene if multiple of same attribute so one total value/attribute
	#Add info on total bp + prop overlap for each attribute
	asm_gff_full_AW_overlap_combined <- asm_gff_full_AW_overlap %>% group_by(asm_name,AW_attributes) %>%
		summarise(tot_bp_overlap=sum(bp_overlap),tot_prop_overlap=sum(prop_overlap))

	#Pivot data format wider
	asm_gff_full_final <- asm_gff_full_AW_overlap_combined %>% select(AW_attributes,asm_name,tot_prop_overlap) %>% tidyr::pivot_wider(.,id_cols=c("asm_name"),names_from='AW_attributes',values_from='tot_prop_overlap',values_fill=0)

	#Add additional information on each gene lost during the transformation
	asm_gff_full_final$AW_Blocks <- unlist(lapply(asm_gff_full_final$asm_name,pull_AW_blocks,bed_df=asm_gff_full_AW_overlap,comp='asm'))
	asm_gff_full_final$chr <- chr_asm
	asm_gff_full_final$start <- unlist(lapply(asm_gff_full_final$asm_name,add_gene_start,gff_exon=chr_asm_gff_exon)) #Fine to use gff_exon here, it'll still extract the first coordinate, so no need to write a new function
	asm_gff_full_final$end <- unlist(lapply(asm_gff_full_final$asm_name,add_gene_end,gff_exon=chr_asm_gff_exon))

	#Add check in case one of the attributes was not present in this chromosome
	asm_gff_full_final <- check_missing_cols(asm_gff_full_final,asm_struct_tag,asm_struct_tag)

	#Reorder
	asm_gff_full_final <- asm_gff_full_final[,c('asm_name','chr','start','end','alignable_region',id_struct_tag,asm_struct_tag,'unalignable','Missing_Data',"AW_Blocks")]
	
	#Obtain final classification for each gene based off of proportion
	asm_gff_full_final <- asm_gff_full_final %>% data.frame() %>% dplyr::mutate(classification=case_when(
		.[[5]] >= 0.95 ~ 'shared',#.[[5]] is column for alignable proportion overlap in full df
		.[[7]] >= 0.95 ~ 'polymorphic',#[[7]] is column for SV in ASM (NAM) proportion overlap in full df
		TRUE ~ 'ambiguous'
		))	

	asm_gff_full_final <- asm_gff_full_final %>% data.frame() %>% arrange(start)

	#Add to full list across all chromosomes
	asm_full_list <- list.append(asm_full_list,asm_gff_full_final)
}

complete_id_exon <- do.call('rbind',id_exon_list)
id_exon_path <- paste('/path/to/polymorphic_gene_calls/B73/exon_calls/',comp_pattern,'_gene_classification_by_exon.tsv',sep='')
write.table(complete_id_exon,id_exon_path,quote=F,col.names=T,row.names=F)

complete_id_full <- do.call('rbind',id_full_list)
id_full_path <- paste('/path/to/polymorphic_gene_calls/B73/full_calls/',comp_pattern,'_gene_classification_by_full.tsv',sep='')
write.table(complete_id_full,id_full_path,quote=F,col.names=T,row.names=F)


complete_asm_exon <- do.call('rbind',asm_exon_list)
asm_exon_path <- paste('/path/to/polymorphic_gene_calls/NAM/exon_calls/',comp_pattern,'_gene_classification_by_exon.tsv',sep='')
write.table(complete_asm_exon,asm_exon_path,quote=F,col.names=T,row.names=F)

complete_asm_full <- do.call('rbind',asm_full_list)
asm_full_path <- paste('/path/to/polymorphic_gene_calls/NAM/full_calls/',comp_pattern,'_gene_classification_by_full.tsv',sep='')
write.table(complete_asm_full,asm_full_path,quote=F,col.names=T,row.names=F)
