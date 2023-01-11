.libPaths('/home/brandvai/mmunasin/Rlibs4')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)

extract_lineage <- function(filename) {
	lineage <- unlist(str_split(unlist(str_split(filename,pattern='_'))[6],pattern='S'))[1]
	return(lineage)
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



NAM_level_orders <- c("B73","B97","Ky21","M162W","Ms71","Oh43","Oh7B","M37W","Mo18W","Tx303","HP301","P39","Il14H","CML52","CML69","CML103","CML228","CML247","CML277","CML322","CML333","Ki3","Ki11","NC350","NC358","Tzi8")
chr_level_orders <- c("chr1","chr2","chr3",'chr4','chr5','chr6','chr7','chr8','chr9','chr10')
classification_level_orders <- c("shared",'ambiguous','polymorphic')


#Load + Format Rolling Window SNP Density Across Lineages
SNP_density_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/SNP_density_summaries'
SNP_files <- list.files(SNP_density_dir, full.names=T)
filename_lineages <- unlist(lapply(SNP_files,extract_lineage))

SNP_data <- lapply(SNP_files,fread)
SNP_data <- Map(cbind,SNP_data,lineage=filename_lineages)
full_SNP_df  <- do.call('rbind',SNP_data)
chr_list <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
full_SNP_df$chr <- factor(full_SNP_df$chr,levels=chr_list)

# Subset it so nonSV base pair >= 900,000 + SNP_Count <= 100
sub_SNP_data <- full_SNP_df %>% subset(nonSV_bp >= 900000)  %>% subset(SNP_count <= 100)

# Split data into contiguous chunks
split_SNP_data <- split(sub_SNP_data, cumsum(c(1,diff(sub_SNP_data$BinStart) != 250000)))

# Keep data if at least 5 continuous blocks in a row
size_split <- unlist(lapply(split_SNP_data,nrow))

contig_split_SNP <- split_SNP_data[size_split >= 5]

chr_vec <- unlist(lapply(contig_split_SNP,grab_chr))
start_vec <- unlist(lapply(contig_split_SNP,grab_start))
end_vec <- unlist(lapply(contig_split_SNP,grab_end))
lineage_vec <- unlist(lapply(contig_split_SNP,grab_lineage))

snp_depleted_regions <- data.frame(chr=chr_vec,start=start_vec,end=end_vec,lineage=lineage_vec)
snp_depleted_regions$start_intersect <- snp_depleted_regions$start + 100000
snp_depleted_regions$end_intersect <- snp_depleted_regions$end - 100000

comp_start_intersect <- c()
comp_end_intersect <- c()

for (row in 1:nrow(snp_depleted_regions)) {
	chr <- as.character(snp_depleted_regions[[row,1]])
	start_coord <- snp_depleted_regions[[row,2]]
	end_coord <- snp_depleted_regions[[row,3]]
	lineage <- snp_depleted_regions[[row,4]]
	start_intersect <- snp_depleted_regions[[row,5]]
	end_intersect <- snp_depleted_regions[[row,6]]#

	dir_name <- paste('B73',lineage,'regions',sep='_')
	file_name <- paste(chr,'B73',lineage,'summarised_regions.tsv',sep='_')
	path_name <- paste('~/TE_Intron/store_data/summarised_AnchorWave_Regions/',dir_name,'/',file_name,sep='')#

	AW_file <- fread(path_name)

	start_intersect_row <- AW_file %>% subset( between(start_intersect,B73_StartPos, B73_EndPos))
	
	B73_start_offset <- start_intersect - start_intersect_row[[1,2]]
	comp_start_intersect_coord <- start_intersect_row[[1,5]] + B73_start_offset
	comp_start_intersect <- c(comp_start_intersect,comp_start_intersect_coord)
	
	end_intersect_row <- AW_file %>% subset( between(end_intersect,B73_StartPos, B73_EndPos))
	
	B73_end_offset <- end_intersect_row[[1,3]] - end_intersect
	comp_end_intersect_coord <- end_intersect_row[[1,6]] - B73_end_offset
	comp_end_intersect <- c(comp_end_intersect,comp_end_intersect_coord)

	# We tested + confirmed that every start + end intersect falls into alignable region
}

snp_depleted_regions$comp_start_intersect <- comp_start_intersect
snp_depleted_regions$comp_end_intersect <- comp_end_intersect
snp_depleted_regions <- snp_depleted_regions %>% 
dplyr::mutate(SNP_Depleted_ID=paste('SNP_Depleted_B73_',lineage,'_',1:n(),sep=''))

B73_snp_depleted_bed <- snp_depleted_regions %>% select(1,5,6,4,9)
colnames(B73_snp_depleted_bed) <- c('chr','start','end','ASM_Comp','SNP_Depleted_ID')

NAM_snp_depleted_bed <- snp_depleted_regions %>% select(1,7,8,4,9)
colnames(NAM_snp_depleted_bed) <- c('chr','start','end','ASM_Comp','SNP_Depleted_ID')

head_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/summarised_AnchorWave_Regions'
sub_dirs <- list.dirs(head_dir,recursive=F)
chr_order <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')

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
			TRUE ~ 'structural_insertion_inNAM'
	))

B73_AW_bed <- all_AW_data %>% 
select(ID_chr,ID_start,ID_end,AW_BlockID,Lineage_Comp,reduced_type) %>%
dplyr::mutate(ID_chr=paste('chr',ID_chr,sep='')) %>% 
dplyr::mutate(size=ID_end-ID_start) %>% 
subset(reduced_type != 'structural_insertion_inNAM')
colnames(B73_AW_bed)[1:3] <- c("chr",'start','end')

NAM_AW_bed <- all_AW_data %>%
select(ASM_chr,ASM_start,ASM_end,AW_BlockID,Lineage_Comp,reduced_type) %>%
dplyr::mutate(ASM_chr=paste('chr',ASM_chr,sep='')) %>% 
dplyr::mutate(size=ASM_end-ASM_start) %>% 
subset(reduced_type != 'structural_insertion_inB73')
colnames(NAM_AW_bed)[1:3] <- c("chr",'start','end')

#B73_AW_insertions <- all_AW_data %>% subset(reduced_type=='structural_insertion_inB73')
#B73_AW_insertions_bed <- B73_AW_insertions %>% 
#select(ID_chr,ID_start,ID_end,AW_BlockID,Lineage_Comp,reduced_type) %>%
#dplyr::mutate(ID_chr=paste('chr',ID_chr,sep='')) %>% 
#dplyr::mutate(size=ID_end-ID_start)
#colnames(B73_AW_insertions_bed)[1:3] <- c("chr",'start','end')

B73_SNPD_AW_overlaps <- list()
for (index in 1:length(unique(B73_snp_depleted_bed$ASM_Comp))) {
	lineage_comp <- unique(B73_snp_depleted_bed$ASM_Comp)[index]

	B73_LC_snp_depleted_bed <- B73_snp_depleted_bed %>% subset(ASM_Comp==lineage_comp)
	B73_LC_AW_bed <- B73_AW_bed %>% subset(Lineage_Comp==lineage_comp)
	B73_LC_snp_AW_overlap <- bedtoolsr::bt.intersect(B73_LC_snp_depleted_bed,B73_LC_AW_bed,wo=T)
	if (nrow(B73_LC_snp_AW_overlap) != 0) {
		colnames(B73_LC_snp_AW_overlap) <- c('SNPD_chr','SNPD_start','SNPD_end','SNPD_Comp','SNPD_ID','AW_chr','AW_start','AW_end','AW_ID','AW_Comp','AW_type','AW_size','bp_overlap')
		B73_LC_snp_AW_overlap <- B73_LC_snp_AW_overlap %>% 
		dplyr::mutate(prop_AW_overlap=bp_overlap/AW_size)
		B73_SNPD_AW_overlaps[[index]] <- B73_LC_snp_AW_overlap
	}
}

B73_SNPD_AW_overlaps <- do.call('rbind',B73_SNPD_AW_overlaps)

fwrite(B73_SNPD_AW_overlaps, '/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SNP_Depleted_Summaries/B73_SNP_Depleted_AW_Overlaps.csv')

NAM_SNPD_AW_overlaps <- list()
for (index in 1:length(unique(NAM_snp_depleted_bed$ASM_Comp))) {
	lineage_comp <- unique(NAM_snp_depleted_bed$ASM_Comp)[index]

	NAM_LC_snp_depleted_bed <- NAM_snp_depleted_bed %>% subset(ASM_Comp==lineage_comp)
	NAM_LC_AW_bed <- NAM_AW_bed %>% subset(Lineage_Comp==lineage_comp)
	NAM_LC_snp_AW_overlap <- bedtoolsr::bt.intersect(NAM_LC_snp_depleted_bed,NAM_LC_AW_bed,wo=T)
	if (nrow(NAM_LC_snp_AW_overlap) != 0) {
		colnames(NAM_LC_snp_AW_overlap) <- c('SNPD_chr','SNPD_start','SNPD_end','SNPD_Comp','SNPD_ID','AW_chr','AW_start','AW_end','AW_ID','AW_Comp','AW_type','AW_size','bp_overlap')
		NAM_LC_snp_AW_overlap <- NAM_LC_snp_AW_overlap %>% 
		dplyr::mutate(prop_AW_overlap=bp_overlap/AW_size)
		NAM_SNPD_AW_overlaps[[index]] <- NAM_LC_snp_AW_overlap
	}
}

NAM_SNPD_AW_overlaps <- do.call('rbind',NAM_SNPD_AW_overlaps)

fwrite(NAM_SNPD_AW_overlaps, '/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SNP_Depleted_Summaries/NAM_SNP_Depleted_AW_Overlaps.csv')

B73_TE_n_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/polymorphic_TE_calls/B73/TE_n'
NAM_TE_n_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/polymorphic_TE_calls/NAM/TE_n'

all_B73_TE_n_calls <- process_TE_calls(B73_TE_n_dir)
all_NAM_TE_n_calls <- process_TE_calls(NAM_TE_n_dir)


B73_TE_bed <- all_B73_TE_n_calls %>% select(2,3,4,19,1,7,10,11,5,18) %>%
dplyr::mutate(size=end-start)

NAM_TE_bed <- all_NAM_TE_n_calls %>% select(2,3,4,19,1,7,10,11,5,18)%>%
dplyr::mutate(size=end-start)

B73_SNPD_TE_overlaps <- list()
for (index in 1:length(unique(B73_snp_depleted_bed$ASM_Comp))) {
	lineage_comp <- unique(B73_snp_depleted_bed$ASM_Comp)[index]

	B73_LC_snp_depleted_bed <- B73_snp_depleted_bed %>% subset(ASM_Comp==lineage_comp)
	B73_LC_TE_bed <- B73_TE_bed %>% subset(ASM_Comp==lineage_comp)
	B73_LC_snp_TE_overlap <- bedtoolsr::bt.intersect(B73_LC_snp_depleted_bed,B73_LC_TE_bed,wo=T)
	if (nrow(B73_LC_snp_TE_overlap) != 0) {
		colnames(B73_LC_snp_TE_overlap) <- c('SNPD_chr','SNPD_start','SNPD_end','SNPD_Comp','SNPD_ID','TE_chr','TE_start','TE_end', 'TE_Lineage_Comp', 'TE_name','TE_class','TE_superfamily','TE_family', 'TE_method','TE_Classification','TE_size','bp_overlap')
		B73_LC_snp_TE_overlap <- B73_LC_snp_TE_overlap %>% 
		dplyr::mutate(prop_TE_overlap=bp_overlap/TE_size)
		B73_SNPD_TE_overlaps[[index]] <- B73_LC_snp_TE_overlap
	}
}

B73_SNPD_TE_overlaps <- do.call('rbind',B73_SNPD_TE_overlaps)

fwrite(B73_SNPD_TE_overlaps, '/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SNP_Depleted_Summaries/B73_SNP_Depleted_TE_Overlaps.csv')

NAM_SNPD_TE_overlaps <- list()
for (index in 1:length(unique(NAM_snp_depleted_bed$ASM_Comp))) {
	lineage_comp <- unique(NAM_snp_depleted_bed$ASM_Comp)[index]

	NAM_LC_snp_depleted_bed <- NAM_snp_depleted_bed %>% subset(ASM_Comp==lineage_comp)
	NAM_LC_TE_bed <- NAM_TE_bed %>% subset(ASM_Comp==lineage_comp)
	NAM_LC_snp_TE_overlap <- bedtoolsr::bt.intersect(NAM_LC_snp_depleted_bed,NAM_LC_TE_bed,wo=T)
	if (nrow(NAM_LC_snp_TE_overlap) != 0) {
		colnames(NAM_LC_snp_TE_overlap) <- c('SNPD_chr','SNPD_start','SNPD_end','SNPD_Comp','SNPD_ID','TE_chr','TE_start','TE_end', 'TE_Lineage_Comp', 'TE_name','TE_class','TE_superfamily','TE_family', 'TE_method','TE_Classification','TE_size','bp_overlap')
		NAM_LC_snp_TE_overlap <- NAM_LC_snp_TE_overlap %>% 
		dplyr::mutate(prop_TE_overlap=bp_overlap/TE_size)
		NAM_SNPD_TE_overlaps[[index]] <- NAM_LC_snp_TE_overlap
	}
}

NAM_SNPD_TE_overlaps <- do.call('rbind',NAM_SNPD_TE_overlaps)

fwrite(NAM_SNPD_TE_overlaps, '/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SNP_Depleted_Summaries/NAM_SNP_Depleted_TE_Overlaps.csv')

