### extract_SNP_density.R
### Manisha Munasinghe - Last Updated -03/03/23
### Take Pairwise AnchorWave Alignments
### Obtain Sliding Window SNP density estimates

.libPaths('/path/to/Rlibs/')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(tibbletime)

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
#[1] tibbletime_0.1.6  tidyr_1.2.0       data.table_1.14.2 stringr_1.4.0
#[5] dplyr_1.0.9#

#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.8.3     fansi_1.0.3      assertthat_0.2.1 utf8_1.2.2
# [5] crayon_1.5.1     R6_2.5.1         DBI_1.1.2        lifecycle_1.0.1
# [9] magrittr_2.0.3   pillar_1.7.0     stringi_1.7.6    rlang_1.0.2
#[13] cli_3.3.0        vctrs_0.4.1      generics_0.1.2   ellipsis_0.3.2
#[17] tools_4.0.4      glue_1.6.2       purrr_0.3.4      compiler_4.0.4
#[21] pkgconfig_2.0.3  tidyselect_1.1.2 tibble_3.1.6

source('/path/to/function_scripts/SNP_density_functions.R')

args <- commandArgs(trailingOnly = TRUE)
AW_dir_name <- args[1]
# We start with the raw AnchorWave data which contains SNPs
AW_files <- list.files(AW_dir_name,full.names=T)

ID_lineage <- unlist(str_split(unlist(str_split(AW_dir_name,pattern='/'))[8],pattern='_'))[1]
ASM_lineage <- unlist(str_split(unlist(str_split(AW_dir_name,pattern='/'))[8],pattern='_'))[2]

# 10 files, one for each c'some
AW_data <- lapply(AW_files,fread)

# Establish Sliding Window Bin
sliding_window <- rollify(sum,window=4)
sliding_window_count_list <- list()

for (i in 1:length(AW_data)) {
	
	chr_AW <- AW_data[[i]]
	chr_id <- paste('chr', unique(chr_AW$B73_Chr),sep='')

	#Pull All SNP data and positions
	chr_AW_SNP <- chr_AW %>% subset(Type=='SNP') %>% select(B73_Chr,B73_StartPos)

	#Split the c'some into 250,000 (250kb) bins
	chr_start <- chr_AW[[1,2]]  #chr_start_end[[i,2]]
	chr_end <- chr_AW[[nrow(chr_AW),3]] #chr_start_end[[i,3]]
	kb_250_windows <- seq(chr_start,chr_end,250000)

	#Extract all 250kb window coordinates (window_id)
	#Get SNP counts for 250kb windows (n)
	#Then get SNP counts for 1Mb sliding in 250kb intervals (SNP_count)
	snp_count_sw <- chr_AW_SNP %>% 
	mutate(window_id=cut(B73_StartPos,kb_250_windows,dig.lab=5)) %>% 
	count(window_id,.drop=FALSE) %>% 
	mutate(SNP_count=sliding_window(n))

	#Make df of SNP count for each window
	final_window_count <- obtain_final_window_SNP_count(snp_count_sw,chr_id)
	
	#Calculate the amount of nonSV sequence in window
	all_alignable_seq <- c()

	chr_AW$row_index <- seq(1,nrow(chr_AW))

	for (j in 1:nrow(final_window_count)) {
		j_bin_start <- final_window_count[[j,2]]
		if (j_bin_start==0) {
			j_bin_start <- 1
		}
		j_bin_end <- final_window_count[[j,3]]

		# Pull out first AW row for the window
		start_row <- chr_AW %>% subset(between(j_bin_start, B73_StartPos,B73_EndPos))
		start_row_index <- start_row[[8]]

		# Pull out the final AW row for the window
		end_row <- chr_AW %>% subset(between(j_bin_end, B73_StartPos,B73_EndPos))
		end_row_index <- end_row[[8]]

		# If the start and and end row are the same, the AW region is just
		# a single entry
		if (all.equal(start_row,end_row) == TRUE) {
			AW_region <- start_row
		} else {
			# Pull out all entries in region
			# and only retain nonvariant and SNP entries
			temp <- chr_AW %>% 
			subset(between(row_index,start_row_index,end_row_index)) %>% 
			filter(Type %in% c('nonvariant_region','SNP'))

			#If there are actually nonvariant and SNP entries
			#get the size of each entry
			if (nrow(temp)!=0) {
				AW_region <- temp %>% 
				select(B73_Chr,B73_StartPos,B73_EndPos,Type) %>% 
				mutate(size=(B73_EndPos - B73_StartPos) + 1)
			}
		}

		#Count the amount of alignable sequence
		alignable_seq <- sum(AW_region$size)
 
 		#Adjust alignable seq if the first or last entry only partially
 		#overlaps the window
		if (alignable_seq != 0) {
			if (j_bin_start != AW_region[[1,2]] & j_bin_start > AW_region[[1,2]]) {
				alignable_seq <- alignable_seq - (j_bin_start - AW_region[[1,2]])
				print(alignable_seq)
			}
			if (j_bin_end!=AW_region[[nrow(AW_region),3]] & j_bin_end < AW_region[[nrow(AW_region),3]]) {
				alignable_seq <- alignable_seq -(AW_region[[nrow(AW_region),3]] - j_bin_end)
				print(alignable_seq)
			}			
		}

		all_alignable_seq <- c(all_alignable_seq,alignable_seq)
	}

	#Bind everything together!
	final_window_count$nonSV_bp <- all_alignable_seq
	final_window_count$norm_SNP_count <- final_window_count$SNP_count/final_window_count$nonSV_bp
	sliding_window_count_list[[i]] <- final_window_count

}

full_sliding_window <- do.call('rbind',sliding_window_count_list)

outfiles_name <- paste('/path/to/SNP_density_summaries/',ID_lineage,'_',ASM_lineage,'SNP_density_summary.tsv',sep='')
write.table(full_sliding_window,file=outfiles_name,quote=F,col.names=T,row.names=F)