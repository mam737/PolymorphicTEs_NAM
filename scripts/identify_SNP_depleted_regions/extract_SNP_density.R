.libPaths('/home/brandvai/mmunasin/Rlibs4')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(tibbletime)

args <- commandArgs(trailingOnly = TRUE)
AW_dir_name <- args[1]
# For testing purpose
#AW_dir_name <- '/home/brandvai/mmunasin/TE_Intron/store_data/AnchorWave_Regions/B73_B97_regions'
# We start with the raw AnchorWave data which contains SNPs
AW_files <- list.files(AW_dir_name,full.names=T)

ID_lineage <- unlist(str_split(unlist(str_split(AW_dir_name,pattern='/'))[8],pattern='_'))[1]
ASM_lineage <- unlist(str_split(unlist(str_split(AW_dir_name,pattern='/'))[8],pattern='_'))[2]

# 10 files, one for each c'some
AW_data <- lapply(AW_files,fread)

sliding_window <- rollify(sum,window=4)
sliding_window_count_list <- list()

for (i in 1:length(AW_data)) {
	
	chr_AW <- AW_data[[i]]

	#Pull All SNP data and positions
	chr_AW_SNP <- chr_AW %>% subset(Type=='SNP') %>% select(B73_Chr,B73_StartPos)

	#Split the c'some into 250,000 (250kb) bins
	chr_start <- chr_AW[[1,2]]  #chr_start_end[[i,2]]
	chr_end <- chr_AW[[nrow(chr_AW),3]] #chr_start_end[[i,3]]
	kb_250_windows <- seq(chr_start,chr_end,250000)


	snp_count_sw <- chr_AW_SNP %>% mutate(window_id=cut(B73_StartPos,kb_250_windows,dig.lab=5)) %>% count(window_id,.drop=FALSE) %>% mutate(SNP_count=sliding_window(n))

	start_bin_final <- unlist(str_split(as.character(snp_count_sw[[nrow(snp_count_sw)-4,1]]),pattern=','))[1]
	start_bin_final <- as.numeric(str_sub(start_bin_final,2,nchar(start_bin_final)))

	end_bin_final <- unlist(str_split(as.character(snp_count_sw[[nrow(snp_count_sw)-1,1]]),pattern=','))[2]
	end_bin_final <- as.numeric(str_sub(end_bin_final,1,nchar(end_bin_final)-1))

	start_bin_vec <- seq(0,start_bin_final,250000)
	end_bin_vec <- seq(1000000,end_bin_final,250000)
	bin_count <- snp_count_sw$SNP_count[seq(4,nrow(snp_count_sw)-1,by=1)]
	chr <- paste('chr', unique(chr_AW$B73_Chr),sep='')

	final_window_count <- data.table(chr=chr,BinStart=start_bin_vec,BinEnd=end_bin_vec,SNP_count=bin_count)
		
	all_alignable_seq <- c()

	chr_AW$row_index <- seq(1,nrow(chr_AW))
	for (j in 1:nrow(final_window_count)) {
		j_bin_start <- final_window_count[[j,2]]
		if (j_bin_start==0) {
			j_bin_start <- 1
		}
		j_bin_end <- final_window_count[[j,3]]


		start_row <- chr_AW %>% subset(between(j_bin_start, B73_StartPos,B73_EndPos))
		start_row_index <- start_row[[8]]

		end_row <- chr_AW %>% subset(between(j_bin_end, B73_StartPos,B73_EndPos))
		end_row_index <- end_row[[8]]

		if (all.equal(start_row,end_row) == TRUE) {
			AW_region <- start_row
		} else {
			temp <- chr_AW %>% subset(between(row_index,start_row_index,end_row_index)) %>% filter(Type %in% c('nonvariant_region','SNP'))

			if (nrow(temp)!=0) {
				AW_region <- temp %>% select(B73_Chr,B73_StartPos,B73_EndPos,Type) %>% mutate(size=(B73_EndPos - B73_StartPos) + 1)
			}
		}

		alignable_seq <- sum(AW_region$size)
 
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

	final_window_count$nonSV_bp <- all_alignable_seq
	final_window_count$norm_SNP_count <- final_window_count$SNP_count/final_window_count$nonSV_bp
	sliding_window_count_list[[i]] <- final_window_count

}

full_sliding_window <- do.call('rbind',sliding_window_count_list)

outfiles_name <- paste('/home/brandvai/mmunasin/TE_Intron/store_data/SNP_density_summaries/',ID_lineage,'_',ASM_lineage,'SNP_density_summary.tsv',sep='')
write.table(full_sliding_window,file=outfiles_name,quote=F,col.names=T,row.names=F)