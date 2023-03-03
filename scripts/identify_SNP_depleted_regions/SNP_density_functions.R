obtain_final_window_SNP_count <- function(snp_count_SW,chr_id) {
	start_bin_final <- unlist(str_split(as.character(snp_count_SW[[nrow(snp_count_SW)-4,1]]),pattern=','))[1]
	start_bin_final <- as.numeric(str_sub(start_bin_final,2,nchar(start_bin_final)))

	end_bin_final <- unlist(str_split(as.character(snp_count_SW[[nrow(snp_count_SW)-1,1]]),pattern=','))[2]
	end_bin_final <- as.numeric(str_sub(end_bin_final,1,nchar(end_bin_final)-1))

	start_bin_vec <- seq(0,start_bin_final,250000)
	end_bin_vec <- seq(1000000,end_bin_final,250000)
	bin_count <- snp_count_SW$SNP_count[seq(4,nrow(snp_count_SW)-1,by=1)]
	final_window_count <- data.table(chr=chr_id,BinStart=start_bin_vec,BinEnd=end_bin_vec,SNP_count=bin_count)

	return(final_window_count)
}

extract_lineage <- function(filename) {
	lineage <- unlist(str_split(unlist(str_split(filename,pattern='_'))[6],pattern='S'))[1]
	return(lineage)
}

obtain_SNP_densities <- function(SNP_dir) {
	files <- list.files(SNP_dir,full.names=T)
	file_lineages <- unlist(lapply(files,extract_lineage))
	SNP_data <- lapply(files,fread)
	SNP_data <- Map(cbind,SNP_data,lineage=file_lineages)
	full_SNP_df <- do.call('rbind',SNP_data)
	full_SNP_df$chr <- factor(full_SNP_df$chr,levels=chr_level_orders)
	return(full_SNP_df)
}

grab_chr <- function(df) {
	return(df[[1,1]])
}

grab_start <- function(df) {
	return(df[[1,2]])
}

grab_end <- function(df) {
	return(df[[nrow(df),3]])
}

grab_lineage <- function(df) {
	return(df[[1,7]])
}

find_SNPdepleted <- function(SNP_data) {
	# Subset it so nonSV base pair >= 900,000 + SNP_Count <= 100
	sub_SNP_data <- SNP_data %>% 
	subset(nonSV_bp >= 900000)  %>% 
	subset(SNP_count <= 100)

	# Split data into contiguous chunks
	split_SNP_data <- split(sub_SNP_data, cumsum(c(1,diff(sub_SNP_data$BinStart) != 250000)))

	# Keep data if at least 5 continuous blocks in a row
	size_split <- unlist(lapply(split_SNP_data,nrow))
	contig_split_SNP <- split_SNP_data[size_split >= 5]

	# Extract position data
	chr_vec <- unlist(lapply(contig_split_SNP,grab_chr))
	start_vec <- unlist(lapply(contig_split_SNP,grab_start))
	end_vec <- unlist(lapply(contig_split_SNP,grab_end))
	lineage_vec <- unlist(lapply(contig_split_SNP,grab_lineage))

	#Create data frame
	snp_depleted_regions <- data.frame(chr=chr_vec,start=start_vec,end=end_vec,lineage=lineage_vec)

	# Offset start and end coordinates by 100,000
	snp_depleted_regions$start_intersect <- snp_depleted_regions$start + 100000
	snp_depleted_regions$end_intersect <- snp_depleted_regions$end - 100000

	return(snp_depleted_regions)
}

extract_comp_intersect_coords <- function(sdf.df) {

	comp_start_intersect <- c()
	comp_end_intersect <- c()

	for (row in 1:nrow(sdf.df)) {
		chr <- as.character(sdf.df[[row,1]])
		start_coord <- sdf.df[[row,2]]
		end_coord <- sdf.df[[row,3]]
		lineage <- sdf.df[[row,4]]
		start_intersect <- sdf.df[[row,5]]
		end_intersect <- sdf.df[[row,6]]#

		dir_name <- paste('B73',lineage,'regions',sep='_')
		file_name <- paste(chr,'B73',lineage,'summarised_regions.tsv',sep='_')
		path_name <- paste('~/TE_Intron/store_data/summarised_AnchorWave_Regions/',dir_name,'/',file_name,sep='')#

		AW_file <- fread(path_name)

		#Extract AW row that start of SNP depleted region falls
		start_intersect_row <- AW_file %>% 
		subset( between(start_intersect,B73_StartPos, B73_EndPos))
	
		# How far(bp) into that region is the start of SNP depleted region 
		B73_start_offset <- start_intersect - start_intersect_row[[1,2]]
		# Obtain equally offset start coord in comparison genotype
		comp_start_intersect_coord <- start_intersect_row[[1,5]] + B73_start_offset
		comp_start_intersect <- c(comp_start_intersect,comp_start_intersect_coord)
	
		# Extract AW region where end of SNP depleted region falls	
		end_intersect_row <- AW_file %>% 
		subset( between(end_intersect,B73_StartPos, B73_EndPos))
	
		#How far(bp) into that region is the end of the SNP depleted region
		B73_end_offset <- end_intersect_row[[1,3]] - end_intersect
		#Obtain equally offset end coord in comparison genotype
		comp_end_intersect_coord <- end_intersect_row[[1,6]] - B73_end_offset
		comp_end_intersect <- c(comp_end_intersect,comp_end_intersect_coord)

		# We tested + confirmed that every start + end intersect falls into alignable region
	}
	sdf.df$comp_start_intersect <- comp_start_intersect
	sdf.df$comp_end_intersect <- comp_end_intersect
	sdf.df <- sdf.df %>% 
	dplyr::mutate(SNP_Depleted_ID=paste('SNP_Depleted_B73_',lineage,'_',1:n(),sep=''))
	return(sdf.df)
}

make_SNPdepleted_bed <- function(sdf.df,tag) {
	if (tag=='B73') {
		snp_depleted_bed <- snp_depleted_regions %>% select(1,5,6,4,9)
	} else {
		snp_depleted_bed <- snp_depleted_regions %>% select(1,7,8,4,9)
	}
	colnames(snp_depleted_bed) <- c('chr','start','end','ASM_Comp','SNP_Depleted_ID')
	return(snp_depleted_bed)
}

obtain_AnchorWave_files <- function(AW_path) {
	sub_dirs <- list.dirs(AW_path,recursive=F)
	all_AW_data <- list()
	for (dir_index in 1:length(sub_dirs)) {
		dir <- sub_dirs[dir_index]
		ID_lineage <- unlist(str_split(unlist(str_split(dir,pattern='/'))[8],pattern='_'))[1]
		ASM_lineage <- unlist(str_split(unlist(str_split(dir,pattern='/'))[8],pattern='_'))[2]
		AW_files <- list.files(dir,full.names=T)
		AW_data <- lapply(AW_files,fread)
		full_AW_data <- do.call('rbind',AW_data)
		full_AW_data$Lineage_Comp <- ASM_lineage
		colnames(full_AW_data) <- c('ID_chr','ID_start','ID_end','ASM_chr','ASM_start','ASM_end','Block_Type','AW_BlockID','Lineage_Comp')
		all_AW_data[[dir_index]] <- full_AW_data
	}

	all_AW_data <- do.call('rbind',all_AW_data)
	all_AW_data <- all_AW_data %>% mutate(ID_BlockSize = (ID_end-ID_start)+1,ASM_BlockSize=(ASM_end-ASM_start)+1 )
	all_AW_data <- all_AW_data %>% mutate(reduced_type = 
		case_when(Block_Type %in% c("Missing_Data","unalignable") ~ 'unalignable',
		Block_Type == "structural_insertion_inB73" ~ 'structural_insertion_inB73',
		Block_Type == 'alignable_region' ~ 'alignable_region',
		TRUE ~ 'structural_insertion_inNAM'))
	return(all_AW_data)
}

make_AW_bed <- function(AW_data,tag) {
	if (tag=='B73') {
		AW_bed <- AW_data %>% 
		select(ID_chr,ID_start,ID_end,AW_BlockID,Lineage_Comp,reduced_type) %>%
		dplyr::mutate(ID_chr=paste('chr',ID_chr,sep='')) %>% 
		dplyr::mutate(size=ID_end-ID_start) %>% 
		subset(reduced_type != 'structural_insertion_inNAM')
	}
	if (tag =='NAM') {
		AW_bed <- AW_data %>% 
		select(ASM_chr,ASM_start,ASM_end,AW_BlockID,Lineage_Comp,reduced_type) %>%
		dplyr::mutate(ASM_chr=paste('chr',ASM_chr,sep='')) %>% 
		dplyr::mutate(size=ASM_end-ASM_start) %>% 
		subset(reduced_type != 'structural_insertion_inB73')
	}
	colnames(AW_bed)[1:3] <- c("chr",'start','end')
	return(AW_bed)
}

pull_SNPDepleted_AW_overlaps <- function(SNPD.bed,AW.bed) {
	SNPD_AW_overlaps <- list()

	for (index in 1:length(unique(SNPD.bed$ASM_Comp))) {
		lineage_comp <- unique(SNPD.bed$ASM_Comp)[index]

		LC_SNPD.bed <- SNPD.bed %>% subset(ASM_Comp==lineage_comp)
		LC_AW.bed <- AW.bed %>% subset(Lineage_Comp==lineage_comp)
		LC_SNPD_AW_overlap <- bedtoolsr::bt.intersect(LC_SNPD.bed,LC_AW.bed,wo=T)

		if (nrow(LC_SNPD_AW_overlap) != 0) {
			colnames(LC_SNPD_AW_overlap) <- c('SNPD_chr','SNPD_start','SNPD_end','SNPD_Comp','SNPD_ID','AW_chr','AW_start','AW_end','AW_ID','AW_Comp','AW_type','AW_size','bp_overlap')
			LC_SNPD_AW_overlap <- LC_SNPD_AW_overlap %>% 
			dplyr::mutate(prop_AW_overlap=bp_overlap/AW_size)
			SNPD_AW_overlaps[[index]] <- LC_SNPD_AW_overlap
		}
	}
	SNPD_AW_overlaps <- do.call('rbind',SNPD_AW_overlaps)
	return(SNPD_AW_overlaps)
}

pull_SNPDepleted_TE_overlaps <- function(SNPD.bed,TE.bed) {
	SNPD_TE_overlaps <- list()

	for (index in 1:length(unique(SNPD.bed$ASM_Comp))) {
		lineage_comp <- unique(SNPD.bed$ASM_Comp)[index]

		LC_SNPD.bed <- SNPD.bed %>% subset(ASM_Comp==lineage_comp)
		LC_TE.bed <- TE.bed %>% subset(ASM_Comp==lineage_comp)
		LC_SNPD_TE_overlap <- bedtoolsr::bt.intersect(LC_SNPD.bed,LC_TE.bed,wo=T)

		if (nrow(LC_SNPD_TE_overlap) != 0) {
			colnames(LC_SNPD_TE_overlap) <- c('SNPD_chr','SNPD_start','SNPD_end','SNPD_Comp','SNPD_ID','TE_chr','TE_start','TE_end', 'TE_Lineage_Comp', 'TE_name','TE_class','TE_superfamily','TE_family', 'TE_method','TE_Classification','TE_size','bp_overlap')
			LC_SNPD_TE_overlap <- LC_SNPD_TE_overlap %>% 
			dplyr::mutate(prop_TE_overlap=bp_overlap/TE_size)
			SNPD_TE_overlaps[[index]] <- LC_SNPD_TE_overlap
		}
	}
	SNPD_TE_overlaps <- do.call('rbind',SNPD_TE_overlaps)
	return(SNPD_TE_overlaps)
}

process_TE_calls <- function(TE_dir){
	n_calls <- lapply(list.files(TE_dir,full.names=T),fread)
	n_file_order <- unlist(lapply(list.files(TE_dir,full.names=T),function(x) unlist(str_split(x,pattern='_'))[7]))
	n_calls <- mapply(cbind,n_calls,'ASM_Comp'=n_file_order,SIMPLIFY=F)
	n_colnames <- c("TE_name","chr","start","end","method","type","class","raw_superfamily","upd_superfamily","condense_superfamily","raw_family","alignable_region","structural_insertion_inID","structural_insertion_inASM","unalignable","Missing_Data","AW_Blocks","classification","ASM_Comp")
	n_calls  <- lapply(n_calls,setNames,n_colnames)

	all_n_calls <- do.call('rbind',n_calls)

	all_n_calls$ASM_Comp <- factor(all_n_calls$ASM_Comp,levels=NAM_level_orders)
	all_n_calls$chr <- factor(all_n_calls$chr,levels=chr_level_orders)
	all_n_calls$classification <- factor(all_n_calls$classification,levels=classification_level_orders)
	all_n_calls <- all_n_calls %>% subset(raw_family!='DTC_ZM00081_consensus')

	return(all_n_calls)
}

make_TE_Bed <- function(TE_calls) {
	TE_bed <- TE_calls %>% select(2,3,4,19,1,7,10,11,5,18) %>%
	dplyr::mutate(size=end-start)
	return(TE_bed)
}
