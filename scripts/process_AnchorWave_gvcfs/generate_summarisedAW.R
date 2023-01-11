#### generate_summarisedAW.R
#### Manisha Munasinghe - Last Updated - 01/11/23
#### Take Parsed AnchorWave Output and Bin Entries
#### as either alignable (which includes nonvariant, SNP, and InDel < 50bp)
#### as structural variant sequence (which are SVs with > 50bp of sequence in one genotype but maps back to a single bp in the other)
#### or unalignable sequence

.libPaths('/path/to/Rlibs/')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(ggplot2)
library(rlist)

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
#[1] rlist_0.4.6.2     ggplot2_3.3.5     tidyr_1.2.0       data.table_1.14.2
#[5] stringr_1.4.0     dplyr_1.0.9#

#loaded via a namespace (and not attached):
# [1] magrittr_2.0.3   munsell_0.5.0    tidyselect_1.1.2 colorspace_2.0-3
# [5] R6_2.5.1         rlang_1.0.2      fansi_1.0.3      tools_4.0.4
# [9] grid_4.0.4       gtable_0.3.0     utf8_1.2.2       cli_3.3.0
#[13] DBI_1.1.2        withr_2.5.0      ellipsis_0.3.2   assertthat_0.2.1
#[17] tibble_3.1.6     lifecycle_1.0.1  crayon_1.5.1     purrr_0.3.4
#[21] vctrs_0.4.1      glue_1.6.2       stringi_1.7.6    compiler_4.0.4
#[25] pillar_1.7.0     generics_0.1.2   scales_1.2.0     pkgconfig_2.0.3

# We restrict Structural Variants to SVs with size = 0 in one genotype
# and size > 50 bp in the other
reformat_SV_type <- function(row) {
	ID_lineage <- unlist(str_split(names(row)[1],pattern='_'))[1]
	ASM_lineage <- unlist(str_split(names(row)[4],pattern='_'))[1]

	supplied_row <- as.character(row)
	ID_size <- as.numeric(supplied_row[3]) - as.numeric(supplied_row[2])
	ASM_size <- as.numeric(supplied_row[6]) - as.numeric(supplied_row[5])

	if (ID_size > 0 & ASM_size > 0) {
		updated_type <- 'unalignable'
	} else if (ID_size > 0 & ASM_size==0) {
		updated_type <- paste('structural_insertion_in',ID_lineage,sep='')
	} else if (ID_size==0 & ASM_size >0) {
		updated_type <- paste('structural_insertion_in',ASM_lineage,sep='')
	}
	return(updated_type)
}

# There are gaps in the AnchorWave alignment that represent Missing Data
# We identify them so we are aware
identify_missing_blocks <- function(df) {
	missing_blocks <- list()
	for (i in 1:(nrow(df)-1) ) {
		if (df[[i,3]] + 1 != df[[i+1,2]]) {

			id_chr <- df[[i,1]]
			id_start_pos <- df[[i,3]] + 1
			id_end_pos <- df[[i+1,2]] - 1
			asm_chr <- df[[i,4]]
			asm_start_pos <- df[[i,6]] + 1
			asm_end_pos <- df[[i+1,5]]-1
			type <- 'Missing_Data'
			missing_row <- c(id_chr,id_start_pos,id_end_pos,asm_chr,asm_start_pos,asm_end_pos,type)
			missing_blocks <- list.append(missing_blocks,missing_row)
		}
	}

	full_missing.df <- do.call('rbind',missing_blocks) %>% data.frame()
	if (nrow(full_missing.df) != 0) {
		colnames(full_missing.df) <- colnames(df)
		full_missing.df[,1] <- as.integer(full_missing.df[,1])
		full_missing.df[,2] <- as.integer(full_missing.df[,2])
		full_missing.df[,3] <- as.integer(full_missing.df[,3])
		full_missing.df[,4] <- as.integer(full_missing.df[,4])
		full_missing.df[,5] <- as.integer(full_missing.df[,5])
		full_missing.df[,6] <- as.integer(full_missing.df[,6])		
	} 
	return(full_missing.df)
}

# Pass in Parsed AnchorWave alignment
# For speed we split this into 
args <- commandArgs(trailingOnly = TRUE)
AW_dir_name <- args[1]

AW_files <- list.files(AW_dir_name,full.names=T)

ID_lineage <- unlist(str_split(unlist(str_split(AW_dir_name,pattern='/'))[8],pattern='_'))[1]
ASM_lineage <- unlist(str_split(unlist(str_split(AW_dir_name,pattern='/'))[8],pattern='_'))[2]

# We use the chromosome endpoints to ensure we cover the entirety of both genomes
# If the AnchorWave pairwise alignment does not extend to the chromosome endpoints
# We will fill that gap and mark it as Missing Data
id_chr_filename <- paste('/home/brandvai/mmunasin/TE_Intron/store_data/NAM_chr_endpoints/Zm-',ID_lineage,'_chr_endpoints.txt',sep='')
id_chr_file <- fread(id_chr_filename)
colnames(id_chr_file) <- c('chr','start','end')

asm_chr_filename <- paste('/home/brandvai/mmunasin/TE_Intron/store_data/NAM_chr_endpoints/Zm-',ASM_lineage,'_chr_endpoints.txt',sep='')
asm_chr_file <- fread(asm_chr_filename)
colnames(asm_chr_file) <- c('chr','start','end')

# Loop through each chromosome pairwise alignment
for (AW_file_input in AW_files) {

	AW_chr_input <- fread(AW_file_input)

	raw_AW_chr <- as.character(AW_chr_input[1,1])
	AW_chr <- paste('chr',raw_AW_chr,sep='')

	ID_AW_startcoord <- as.numeric(AW_chr_input[1,2])
	ASM_AW_startcoord <- as.numeric(AW_chr_input[1,5])

	ID_AW_finalcoord <- as.numeric(AW_chr_input[nrow(AW_chr_input),3])
	ASM_AW_finalcoord <- as.numeric(AW_chr_input[nrow(AW_chr_input),6])

	# Retain true SVs (size = 0 in one genotype, size > 50 bp in the other)
	# All others classified as unalignable
	AW_SV <- AW_chr_input %>% subset(Type=='structural_variant')
	AW_SV$Reformat_Type <- unlist(apply(AW_SV,1,reformat_SV_type))

	# Check for any missing data in the pairwise alignment
	AW_missing <- identify_missing_blocks(df=AW_chr_input)

	if (nrow(AW_missing) !=0) {
		AW_missing$Reformat_Type <- 'Missing_Data'
		AW_SV <- rbind(AW_SV,AW_missing) %>% arrange(as.numeric(B73_StartPos))		
	}

	# We now fill in the remaining portions of the genome that are neither
	# true SVs or unalignable with alignable regions
	# Functionally combining nonvariant, SNP, and small InDel regions
	if (ID_AW_startcoord != as.numeric(AW_SV[1,2])) {
		id_chr <- as.numeric(AW_SV[1,1])
		id_row_start <- ID_AW_startcoord
		id_row_end <- as.numeric(AW_SV[1,2]) -1
		asm_chr <- as.numeric(AW_SV[1,4])
		asm_row_start <- ASM_AW_startcoord
		asm_row_end <- as.numeric(AW_SV[1,5]) -1	

		AW_SV <- AW_SV %>% add_row(!!!setNames(list(id_chr,id_row_start,id_row_end,asm_chr,asm_row_start,asm_row_end,'alignable_region','alignable_region'), names(.))) %>% arrange(.[[1]],.[[2]])		
	}

	if (ID_AW_finalcoord != as.numeric(AW_SV[nrow(AW_SV),3])) {
		id_chr <- as.numeric(AW_SV[nrow(AW_SV),1])
		id_row_start <- as.numeric(AW_SV[nrow(AW_SV),3]) + 1
		id_row_end <- ID_AW_finalcoord
		asm_chr <- as.numeric(AW_SV[nrow(AW_SV),4])
		asm_row_start <- as.numeric(AW_SV[nrow(AW_SV),6]) + 1
		asm_row_end <- ASM_AW_finalcoord	

		AW_SV <- AW_SV %>% add_row(!!!setNames(list(id_chr,id_row_start,id_row_end,asm_chr,asm_row_start,asm_row_end,'alignable_region','alignable_region'), names(.))) %>% arrange(.[[1]],.[[2]])
	}

	missing_rows <- data.frame(ID_Chr=as.numeric(),ID_StartPos=as.numeric(),
	ID_EndPos=as.numeric(),ASM_Chr=as.numeric(),ASM_StartPos=as.numeric(),
	ASM_EndPos=as.numeric(),Type=as.character(),Reformat_Type=as.character())	

	# We now move through our df, fillling in nonvariant regions between
	# our true SV and unalignable regions
	if (nrow(AW_SV) > 1) {
  		for (i in seq(1,nrow(AW_SV)-1)) {
    		junct_end <- as.numeric(AW_SV[i,3]) +1
    		junct_start <- as.numeric(AW_SV[i+1,2])
      
	    	if (junct_end!=junct_start) {
	      		id_chr <- as.numeric(AW_SV[nrow(AW_SV),1])
	      		id_row_start <- as.numeric(AW_SV[i,3]) + 1
	      		id_row_end <- as.numeric(AW_SV[i+1,2]) -1
	      		asm_chr <- as.numeric(AW_SV[nrow(AW_SV),4])
	      		asm_row_start <- as.numeric(AW_SV[i,6]) + 1
	      		asm_row_end <- as.numeric(AW_SV[i+1,5]) -1
	          
	      		missing_rows <- missing_rows %>% 
	        		add_row(!!!setNames(list(id_chr,id_row_start,id_row_end,asm_chr,asm_row_start,asm_row_end,'alignable_region','alignable_region'), names(.))) %>% 
	        		arrange(.[[1]],.[[2]])
	    	}
  		}
	}

	colnames(missing_rows) <- colnames(AW_SV)

	# Combine and give each block a unique identifier
	final_region <- rbind(AW_SV,missing_rows) %>% arrange(.[[1]],.[[2]]) %>% select(-c(Type))
	
	# Check Chromosome Start and End match
	id_official_start <- id_chr_file %>% subset(chr==AW_chr) %>% pull(start)
	id_official_end <- id_chr_file %>% subset(chr==AW_chr) %>% pull(end)

	asm_official_start <- asm_chr_file %>% subset(chr==AW_chr) %>% pull(start)
	asm_official_end <- asm_chr_file %>% subset(chr==AW_chr) %>% pull(end)

	if (id_official_start != final_region[[1,2]] | asm_official_start != final_region[[1,5]]) {
		row <- list(raw_AW_chr,id_official_start,final_region[[1,2]]-1,raw_AW_chr,asm_official_start,final_region[[1,5]]-1,'Missing_Data')
		final_region <- rbind(final_region,row)	%>% arrange(.[[1]],.[[2]])
	}

	if (id_official_end != final_region[[nrow(final_region),3]] | asm_official_end != final_region[[nrow(final_region),6]]) {
		row <- list(raw_AW_chr,final_region[[nrow(final_region),3]]+1,id_official_end,raw_AW_chr,final_region[[nrow(final_region),6]]+1,asm_official_end,'Missing_Data')
		final_region <- rbind(final_region,row)	%>% arrange(.[[1]],.[[2]])
	} 


	final_region$AW_Block <- paste(paste('chr',final_region$B73_Chr,sep=''),'AW_BlockID',seq(1,nrow(final_region)),sep='_')


	#### OUTPUT SUMMARISED AW REGIONS
	output_dir <- '/path/to/output_dir/summarised_AnchorWave_Regions/'
	output_subdir <- paste(ID_lineage,ASM_lineage,'regions',sep='_')

	final_output_dir <- paste(output_dir,output_subdir,sep='')
	if (dir.exists(final_output_dir)!=TRUE) {
		dir.create(final_output_dir)
	}
	chr_header <- paste('chr',as.character(unique(final_region[,1])),sep='')
	output_filename <- paste(chr_header,ID_lineage,ASM_lineage,'summarised_regions.tsv',sep='_')

	output_path <- paste(final_output_dir,output_filename,sep='/')
	write.table(final_region,output_path,quote=F,col.names=T,row.names=F)
}
