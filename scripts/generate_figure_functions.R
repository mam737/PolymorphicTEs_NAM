load_filtered_TE_data <- function(TE_dir) {
	filtered_TE_anno <- list()

	for (dir_index in 2:length(list.dirs(TE_dir))) {
		subdir <- list.dirs(TE_dir)[dir_index]
		lineage <- unlist(str_split(subdir,pattern='/'))[8]
		lineage_TEs <- do.call('rbind',lapply(list.files(subdir,full.names=T),fread))
		filtered_TE_anno[[dir_index]] <- lineage_TEs
	}
	filtered_TE_anno <- filtered_TE_anno[-1]
	filtered_TE_anno <- do.call('rbind',filtered_TE_anno)
	filtered_TE_anno <- filtered_TE_anno %>% subset(raw_family!='DTC_ZM00081_consensus') %>% 
	dplyr::mutate(size=end-start)
	return(filtered_TE_anno)
}

filtered_TE_summary <- function(filtered_TE) {
	TE_summary <- filtered_TE %>% group_by(lineage) %>% 
	dplyr::summarise(TotalTE=n(),TotalMb=sum(size)/1000000)

	TE_summary <- filtered_TE %>% group_by(lineage,method) %>% 
	dplyr::summarise(TotalTE=n(),TotalMb=sum(size)/1000000) %>%
	pivot_wider(names_from=method,values_from=c(TotalTE,TotalMb),names_vary='slowest') %>%
	left_join(x=TE_summary,by='lineage')

	TE_summary <- filtered_TE %>% group_by(lineage,class) %>% 
	dplyr::summarise(TotalTE=n(),TotalMb=sum(size)/1000000) %>%
	pivot_wider(names_from=class,values_from=c(TotalTE,TotalMb),names_vary='slowest') %>%
	left_join(x=TE_summary,by='lineage')

	TE_summary <- filtered_TE %>% group_by(lineage,condense_superfamily) %>% 
	dplyr::summarise(TotalTE=n(),TotalMb=sum(size)/1000000) %>%
	pivot_wider(names_from=condense_superfamily,values_from=c(TotalTE,TotalMb),names_vary='slowest') %>%
	left_join(x=TE_summary,by='lineage')
	return(TE_summary)
} 

load_AW_data <- function(AW_dir) {
	sub_dirs <- list.dirs(AW_dir,recursive=F)
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
	all_AW_data <- all_AW_data %>% 
	mutate(ID_BlockSize = (ID_end-ID_start)+1,ASM_BlockSize=(ASM_end-ASM_start)+1 )

	all_AW_data <- all_AW_data %>% mutate(reduced_type = 
		case_when(Block_Type %in% c("Missing_Data","unalignable") ~ 'unalignable',
			Block_Type == "structural_insertion_inB73" ~ 'structural_insertion_inB73',
			Block_Type == 'alignable_region' ~ 'alignable_region',
			TRUE ~ 'structural_insertion_inNAM'
		))
	return(all_AW_data)
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

summarise_TE_Tag <- function(n_calls,tag_type) {
	summarised_calls <- n_calls %>%
	group_by(ASM_Comp,classification) %>% 
	dplyr::summarise(TE_ClassificationTotal=n()) %>% 
	dplyr::mutate(tag=tag_type)
}

process_gene_calls <- function(gene_dir,tag,type) {
	gene_calls <- lapply(list.files(gene_dir,full.names=T),fread)
	gene_file_order <- unlist(lapply(list.files(gene_dir,full.names=T),function(x) unlist(str_split(x,pattern='_'))[7]))

	id_exon_colnames <- c("id_name",'id_chr','id_start','id_end','exon_num','alignable_region','structural_insertion_inID','structural_insertion_inASM','unalignable','Missing_Data','AW_Blocks','classification')
	id_gene_colnames <- c("id_name",'id_chr','id_start','id_end','alignable_region','structural_insertion_inID','structural_insertion_inASM','unalignable','Missing_Data','AW_Blocks','classification')

	asm_exon_colnames <- c("asm_name",'asm_chr','asm_start','asm_end','exon_num','alignable_region','structural_insertion_inID','structural_insertion_inASM','unalignable','Missing_Data','AW_Blocks','classification')
	asm_gene_colnames <- c("asm_name",'asm_chr','asm_start','asm_end','alignable_region','structural_insertion_inID','structural_insertion_inASM','unalignable','Missing_Data','AW_Blocks','classification')
 
	if (tag=='id' & type =='exon') {
		gene_calls  <- lapply(gene_calls,setNames,id_exon_colnames)
	} else if (tag=='id' & type =='gene') {
		gene_calls  <- lapply(gene_calls,setNames,id_gene_colnames)
	} else if (tag=='asm' & type =='exon') {
		gene_calls  <- lapply(gene_calls,setNames,asm_exon_colnames)
	} else {
		gene_calls  <- lapply(gene_calls,setNames,asm_gene_colnames)
	}

	gene_calls <- mapply(cbind,gene_calls,'ASM_Comp'=gene_file_order,SIMPLIFY=F)
	all_gene_calls <- do.call('rbind',gene_calls)

	all_gene_calls$ASM_Comp <- factor(all_gene_calls$ASM_Comp,levels=NAM_level_orders)

	if (tag=='id') {
		all_gene_calls$id_chr <- factor(all_gene_calls$id_chr,levels=chr_level_orders)
		all_gene_calls$classification <- factor(all_gene_calls$classification,levels=classification_level_orders)
		all_gene_calls <- all_gene_calls %>% separate(id_name,into=c('gene_name','transcript'),sep='_')
	} else {
		all_gene_calls$asm_chr <- factor(all_gene_calls$asm_chr,levels=chr_level_orders)
			all_gene_calls$classification <- factor(all_gene_calls$classification,levels=classification_level_orders)
		all_gene_calls <- all_gene_calls %>% separate(asm_name,into=c('gene_name','transcript'),sep='_')
	}
	
	return(all_gene_calls)
}

obtain_gene_synteny <- function(dir_name) {
	syntenic_gene_files <- list.files(dir_name,full.names=T)
	syntenic_genes_data <- lapply(syntenic_gene_files,fread)
	syntenic_file_order <- unlist(lapply(syntenic_gene_files, function(x) unlist(str_split(unlist(str_split(x,pattern='/'))[8],pattern='-'))[2]))

	syntenic_genes_data <- mapply(cbind,syntenic_genes_data,'NAM_Line'=syntenic_file_order,SIMPLIFY=F)
	syntenic_colnames <- c("Sorghum_gene_name","NAM_gene_name",'NAM_Line')
	syntenic_genes_data <- lapply(syntenic_genes_data,setNames,syntenic_colnames)

	syntenic_genes_all <- do.call('rbind',syntenic_genes_data)
	return(syntenic_genes_all)
}

add_synteny_column <- function(gene.df,syntenic_genes.df) {
	gene.df <- gene.df %>% 
	dplyr::mutate(syntenic_gene = case_when(
		gene_name %in% syntenic_genes.df$NAM_gene_name ~ 'TRUE',
		TRUE ~ 'FALSE'
		))
	return(gene.df)
}

obtain_synteny_summary <- function(gene.df,synteny_tag,type_tag) {
	if (synteny_tag==TRUE & type_tag=='id') {
		summary_gene.df <- gene.df %>%
		subset(syntenic_gene==TRUE) %>% 
		group_by(ASM_Comp,classification) %>%
		dplyr::summarise(gene_ClassificationTotal = n()) %>% 
		dplyr::mutate(tag='B73')
	} else if (synteny_tag==TRUE & type_tag=='asm') {
		summary_gene.df <- gene.df %>%
		subset(syntenic_gene==TRUE) %>% 
		group_by(ASM_Comp,classification) %>%
		dplyr::summarise(gene_ClassificationTotal = n()) %>% 
		dplyr::mutate(tag='NAM')		
	} else if (synteny_tag==FALSE & type_tag=='id') {
		summary_gene.df <- gene.df %>%
		subset(syntenic_gene==FALSE) %>% 
		group_by(ASM_Comp,classification) %>%
		dplyr::summarise(gene_ClassificationTotal = n()) %>% 
		dplyr::mutate(tag='B73')
	} else {
		summary_gene.df <- gene.df %>%
		subset(syntenic_gene==FALSE) %>% 
		group_by(ASM_Comp,classification) %>%
		dplyr::summarise(gene_ClassificationTotal = n()) %>% 
		dplyr::mutate(tag='NAM')		
	}
	return(summary_gene.df)
}

shared_frequency <- function(call_dist.df) {
	call_dist.df <- call_dist.df %>% dplyr::mutate(degree_shared = 
	case_when(n==25 ~ 'Core',
		n %in% c(23,24) ~ 'Near Core',
		n %in% seq(1,22) ~ 'Variable',
		TRUE ~ 'Private'
		))
	call_dist.df$n <- as.factor(call_dist.df$n)
	call_dist.df$degree_shared <- factor(call_dist.df$degree_shared,levels=c('Private','Variable','Near Core','Core'))
	return(call_dist.df)
}

obtain_polymorphic_percent <- function(call.df,tag.string) {
	poly_percent <- call.df %>%
	group_by(ASM_Comp,classification,condense_superfamily,method) %>%
	dplyr::summarise(n=n()) %>% ungroup() %>% 
	group_by(ASM_Comp,condense_superfamily,method) %>% 
	dplyr::mutate(percent_total=n/sum(n)*100) %>% dplyr::mutate(tag=tag.string) %>% 
	subset(classification=='polymorphic') %>% 
	summarySE(measurevar="percent_total", groupvars=c("condense_superfamily","method")) %>%
	dplyr::mutate(tag=tag.string)
	return(poly_percent)
}

obtain_large_polymorphic_percent <- function(call.df,tag.string) {
	largefam_poly_percent <- call.df %>%
	group_by(ASM_Comp,raw_family) %>%
	dplyr::mutate(fam_size=n()) %>% subset(fam_size >= 20) %>% 
	ungroup() %>% group_by(ASM_Comp,raw_family,classification,fam_size)%>%
	dplyr::summarise(class_count=n()) %>% ungroup() %>% 
	group_by(ASM_Comp,raw_family) %>%
	dplyr::mutate(percent_total=class_count/sum(class_count)*100) %>% 
	subset(classification=='polymorphic') %>% mutate(tag=tag.string)
	return(largefam_poly_percent)
}

summarise_AW_CA <- function(CA.df,tag.string) {
	summary_AW_CA <- CA.df %>%
	pivot_longer(cols=c('n','Mb'),names_to='Type',values_to='Amount') %>%
	tidyr::separate(Category,into=c("tag",'number'),sep='_') %>%
	dplyr::mutate(Category=paste(tag,number,sep=' ')) %>%
	dplyr::mutate(Category=factor(Category,levels=c('Category 1','Category 2','Category 3','Category 4','Category 5')))%>%
	dplyr::mutate(Type=case_when(
		Type=='n'~'Number',
		TRUE ~ 'Cumulative Mb'))%>%
	dplyr::mutate(Type=factor(Type,levels=c('Number','Cumulative Mb')))%>% 
	select(Category,Type,Amount) %>% 
	dplyr::mutate(tag=tag.string)
}

fam_supfam_relate <- function(B73_TE,NAM_TE) {
	fam_supfam_relation <- rbind(B73_TE,NAM_TE) %>% 
	group_by(raw_family, condense_superfamily) %>% dplyr::summarise(n=n()) %>%
	ungroup() %>% group_by(raw_family) %>% arrange(desc(n)) %>% 
	filter(row_number()==1) %>% select(1,2)
	colnames(fam_supfam_relation)[1] <- 'TE_family'
	return(fam_supfam_relation)
}

TE_SV_relationship <- function(type_SNPD_TE_overlaps,type_SNPD_SV_categories) {
	polymorphic_TEs_in_SNPD <- type_SNPD_TE_overlaps %>%
	subset(TE_Classification=='polymorphic')
	colnames(polymorphic_TEs_in_SNPD)[9] <- 'TE_LC'

	polymorphic_TEs_in_SNPD.bed <- polymorphic_TEs_in_SNPD %>%
	select(6,7,8,9,10,11,12,13,14,15,16,5)
	colnames(polymorphic_TEs_in_SNPD.bed)[1:3] <- c('chr','start','end')

	SNPD_insertion_categories.bed <- type_SNPD_SV_categories %>%
	select(6,7,8,9,10,12,15,5)
	colnames(SNPD_insertion_categories.bed)[1:3] <- c('chr','start','end')

	type_polymorphicTE_InsertionCategory_SNPD <- list()

	for (i in 1:length(unique(polymorphic_TEs_in_SNPD.bed$TE_LC))) {
		LC <- unique(polymorphic_TEs_in_SNPD.bed$TE_LC)[i]
		LC_polymorphic_TEs <- polymorphic_TEs_in_SNPD.bed %>% subset(TE_LC==LC)
		LC_insertion_categories <- SNPD_insertion_categories.bed %>% subset(AW_Comp==LC)
		polymorphicTE_insertion_overlap <- bedtoolsr::bt.intersect(LC_polymorphic_TEs,LC_insertion_categories,wo=T)
		colnames(polymorphicTE_insertion_overlap) <- c('TE_chr','TE_start','TE_end','TE_LC','TE_name','TE_class','TE_superfamily','TE_family','TE_method','TE_classification','TE_size','TE_SNPD_ID','AW_chr','AW_start','AW_end','AW_BlockID','AW_LC','AW_Size','AW_Category','AW_SNPD_ID','bp_overlap')
		type_polymorphicTE_InsertionCategory_SNPD[[i]] <- polymorphicTE_insertion_overlap
	}
	polymorphicTE_InsertionCategory_SNPD <- do.call('rbind',type_polymorphicTE_InsertionCategory_SNPD)
	return(polymorphicTE_InsertionCategory_SNPD)
}

obtain_cumulMB <- function(AW_data) {
	# Pull B73 + NAM genome wide proportions of 
	# SVs/Total from Genome Wide AW Data
	B73_genome_cumulMb <- AW_data %>% 
	subset(reduced_type != 'structural_insertion_inNAM') %>% 
	group_by(reduced_type) %>% 
	dplyr::summarise(ID_cumulMb=sum(ID_BlockTotal_Mb))
	B73_total_cumulMb <- sum(B73_genome_cumulMb$ID_cumulMb)
	B73_insertion_cumulMb <- B73_genome_cumulMb$ID_cumulMb[2]

	NAM_genome_cumulMb <- AW_data %>% 
	subset(reduced_type != 'structural_insertion_inB73') %>% 
	group_by(reduced_type) %>% 
	dplyr::summarise(ASM_cumulMb=sum(ASM_BlockTotal_Mb))
	NAM_total_cumulMb <- sum(NAM_genome_cumulMb$ASM_cumulMb)
	NAM_insertion_cumulMb <- NAM_genome_cumulMb$ASM_cumulMb[2]

	total_cumulMb <- B73_total_cumulMb+NAM_total_cumulMb
	insertion_cumulMb <- B73_insertion_cumulMb+NAM_insertion_cumulMb

	cumulMb <- data.frame(group='Genome-Wide',perc_inserted=(insertion_cumulMb/total_cumulMb)*100,Mb_inserted=insertion_cumulMb)

	# Read in data for AW Blocks in SNP Depleted Regions
	# and do the same thing as genome wide
	B73_SNPD_AW_overlaps <- fread('/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SNP_Depleted_Summaries/B73_SNP_Depleted_AW_Overlaps.csv')
	NAM_SNPD_AW_overlaps <- fread('/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SNP_Depleted_Summaries/NAM_SNP_Depleted_AW_Overlaps.csv')

	B73_SNPD_cumulMb <- B73_SNPD_AW_overlaps %>% 
	group_by(AW_type) %>% 
	dplyr::summarise(ID_cumulMb=sum(bp_overlap)/1000000)
	B73_SNPD_total_cumulMb <- sum(B73_SNPD_cumulMb$ID_cumulMb)
	B73_SNPD_insertion_cumulMb <- B73_SNPD_cumulMb$ID_cumulMb[2]	

	NAM_SNPD_cumulMb <- NAM_SNPD_AW_overlaps %>% 
	group_by(AW_type) %>% 
	dplyr::summarise(ASM_cumulMb=sum(bp_overlap)/1000000)
	NAM_SNPD_total_cumulMb <- sum(NAM_SNPD_cumulMb$ASM_cumulMb)
	NAM_SNPD_insertion_cumulMb <- NAM_SNPD_cumulMb$ASM_cumulMb[2]

	SNPD_cumulMb <- B73_SNPD_total_cumulMb+NAM_SNPD_total_cumulMb
	SNPD_insertion_cumulMb <-B73_SNPD_insertion_cumulMb+NAM_SNPD_insertion_cumulMb

	#Bind Data Together For Plotting
	cumulMb <- rbind(cumulMb, data.frame(group='SNP Depleted',perc_inserted=(SNPD_insertion_cumulMb/SNPD_cumulMb)*100,Mb_inserted=SNPD_insertion_cumulMb))
	return(cumulMb)
}

obtain_TE_cumulMb <- function(B73_TE,NAM_TE) {
	B73_TE_cumulMb <- B73_TE %>% 
	dplyr::mutate(size=end-start) %>% 
	group_by(classification) %>% 
	dplyr::summarise(TE_cumulMB=sum(size)/1000000)	

	B73_totalTE_cumulMb <- sum(B73_TE_cumulMb$TE_cumulMB)
	B73_polymorphicTE_cumulMb <- B73_TE_cumulMb$TE_cumulMB[2]	

	NAM_TE_cumulMb <- NAM_TE %>% 
	dplyr::mutate(size=end-start) %>% 
	group_by(classification) %>% 
	dplyr::summarise(TE_cumulMB=sum(size)/1000000)	

	NAM_totalTE_cumulMb <- sum(NAM_TE_cumulMb$TE_cumulMB)
	NAM_polymorphicTE_cumulMb <- NAM_TE_cumulMb$TE_cumulMB[2]	

	totalTE_cumulMb <- B73_totalTE_cumulMb+NAM_totalTE_cumulMb
	polymorphicTE_cumulMb <- B73_polymorphicTE_cumulMb+NAM_polymorphicTE_cumulMb	

	TE_cumulMb <- data.frame(group='Genome-Wide',perc_polymorphic=(polymorphicTE_cumulMb/totalTE_cumulMb)*100,Mb_polymorphic=polymorphicTE_cumulMb)

	B73_SNPD_TE_overlaps <- fread('/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SNP_Depleted_Summaries/B73_SNP_Depleted_TE_Overlaps.csv')
	NAM_SNPD_TE_overlaps <- fread('/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SNP_Depleted_Summaries/NAM_SNP_Depleted_TE_Overlaps.csv')	

	B73_SNPD_TE_cumulMb <- B73_SNPD_TE_overlaps %>% 
	group_by(TE_Classification) %>% 
	dplyr::summarise(TE_cumulMB=sum(TE_size)/1000000)	

	B73_SNPD_totalTE_cumulMb <- sum(B73_SNPD_TE_cumulMb$TE_cumulMB)
	B73_SNPD_polymorphicTE_cumulMb <- B73_SNPD_TE_cumulMb$TE_cumulMB[2]	

	NAM_SNPD_TE_cumulMb <- NAM_SNPD_TE_overlaps %>% 
	group_by(TE_Classification) %>% 
	dplyr::summarise(TE_cumulMB=sum(TE_size)/1000000)	

	NAM_SNPD_totalTE_cumulMb <- sum(NAM_SNPD_TE_cumulMb$TE_cumulMB)
	NAM_SNPD_polymorphicTE_cumulMb <- NAM_SNPD_TE_cumulMb$TE_cumulMB[2]	

	SNPD_totalTE_cumulMb <- B73_SNPD_totalTE_cumulMb+NAM_SNPD_totalTE_cumulMb
	SNPD_polymorphicTE_cumulMb <- B73_SNPD_polymorphicTE_cumulMb+NAM_SNPD_polymorphicTE_cumulMb	

	TE_cumulMb <- rbind(TE_cumulMb, data.frame(group='SNP Depleted',perc_polymorphic=(SNPD_polymorphicTE_cumulMb/SNPD_totalTE_cumulMb)*100, Mb_polymorphic=SNPD_polymorphicTE_cumulMb))
	return(TE_cumulMb)
}


all_B73_full_calls %>% 
group_by(gene_name,classification,.drop=F) %>% 
dplyr::summarise(n=n()) %>% 
subset(classification=='shared')

B73_gene_synteny <- all_B73_full_calls %>% 
group_by(gene_name,id_chr,id_start,id_end,syntenic_gene) %>%
dplyr::summarise(n=n()) %>% dplyr::mutate(lineage='B73') %>%
select(id_chr,id_start,id_end,gene_name,lineage,syntenic_gene)
colnames(B73_gene_synteny) <- c('chr','start','end','gene_name','lineage','syntenic_gene')

NAM_gene_synteny <- all_NAM_full_calls %>% 
group_by(gene_name,asm_chr,asm_start,asm_end,syntenic_gene,ASM_Comp) %>%
dplyr::summarise(n=n()) %>%
select(asm_chr,asm_start,asm_end,gene_name,ASM_Comp,syntenic_gene)
colnames(NAM_gene_synteny) <- c('chr','start','end','gene_name','lineage','syntenic_gene')

gene_synteny <- rbind(B73_gene_synteny,NAM_gene_synteny)

for (i in NAM_level_orders) {
	lineage_gene_synteny <- gene_synteny %>% filter(lineage==i)
	lineage_TE <- filtered_TE_anno %>% filter(lineage==i)
}




######################################################

blast_SNPD_left_fledge_bed <- extract_needle_files(B73_SNPD_SV_categories %>% filter(AW_size >= 200),X=100,tag='left_blast')
blast_SNPD_right_fledge_bed <- extract_needle_files(B73_SNPD_SV_categories %>% filter(AW_size >= 200),X=200,tag='right_blast')
write_fledge_beds(blast_SNPD_left_fledge_bed,X=100,tag='blast',fledge_type='left_blast')
write_fledge_beds(blast_SNPD_right_fledge_bed,X=100,tag='blast',fledge_type='right_blast')

# From these bed files, we can extract the actual sequence of these fledges
# Extract sequence
# bedtools getfasta -fi /home/springer/shared/ns-genome/Zmays_B73v5/10.fasta -bed SNPD_I100_right_interior.bed -nameOnly -fo SNPD_I100_right_interior.fa.out
# Split into distinct fastas for each fledge
# awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }' ../extract_fastas/SNPD_I100_right_interior.fa.out

#bedtools getfasta -fi /home/springer/shared/ns-genome/Zmays_B73v5/10.fasta -bed blast_SNPD_right_blast.bed -nameOnly -fo blast_SNPD_right_blast.fa.out
#awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }' ../extract_fastas/blast_SNPD_left_blast.fa.out

# Generate a param file which will be passed to a bash script
# that runs all of the needle pairwise alignments in parallel
# that runs all of the dotplot pairwise figs in parallel
left_blast_fasta_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/blast_edges/blast_SNPD/blast_SNPD_left_fledges'
right_blast_fasta_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/blast_edges/blast_SNPD/blast_SNPD_right_fledges'
generate_blast_param(left_blast_fasta_dir,right_blast_fasta_dir,type='blastn')

blast_results_files <- list.files('/home/brandvai/mmunasin/TE_Intron/store_data/blast_edges/blast_SNPD/blast_outfiles',full.names=T)
blast_results <- lapply(blast_results_files,fread)

#Subset actual results
SNPD_blast_hits <- blast_results[sapply(blast_results, function(x) dim(x)[1]) > 0]
SNPD_blast_hits <- do.call('rbind',SNPD_blast_hits)
colnames(SNPD_blast_hits) <- c('query_file','subject_file','perc_identity','align_length','num_mistmatches','gap_opens','q_start','q_end','s_start','s_end','evalue','bit_score')
SNPD_blast_hits <- SNPD_blast_hits %>%
dplyr::mutate(orient=case_when(
	q_end > q_start & s_end > s_start ~ 'Alignment',
	q_start > q_end & s_start > s_end ~ 'Alignment',
	q_end > q_start & s_start > s_end ~ 'Inverted_Alignment',
	q_start > q_end & s_end > s_start ~ 'Inverted_Alignment'
	))

SNPD_blast_hits$AW_ID <- unlist(lapply(SNPD_blast_hits$query_file, function(x) paste(unlist(str_split(x,pattern=':'))[4],collapse='_')))
SNPD_blast_hits$AW_Comp <- unlist(lapply(SNPD_blast_hits$query_file, function(x) paste(unlist(str_split(x,pattern=':'))[5],collapse='_')))
unique_id <- paste(SNPD_blast_hits$AW_ID,SNPD_blast_hits$AW_Comp,sep='_')

SNPD_blast_hits <- B73_SNPD_SV_categories %>% select(AW_ID,AW_Comp,Category) %>%
left_join(x=SNPD_blast_hits,y=.,by=c('AW_ID','AW_Comp')) %>%
group_by(AW_ID,AW_Comp) %>%
filter(row_number()==1) %>%
select(-query_file,-subject_file)

TIR_put_insertions <- SNPD_blast_hits %>% 
filter(orient=='Inverted_Alignment') %>% 
select(AW_ID,AW_Comp,Category) %>% 
dplyr::mutate(unique_ID=paste(AW_ID,AW_Comp,sep='_'))


plot1 <- SNPD_blast_hits %>% 
select(AW_ID,AW_Comp,Category,perc_identity,align_length,orient) %>% 
dplyr::mutate(perc_identity=case_when(
	is.na(perc_identity)~0,
	TRUE ~ perc_identity)) %>%
dplyr::mutate(orient=case_when(
	is.na(orient)~'No_Alignment',
	TRUE ~ orient)) %>%
dplyr::mutate(orient=factor(orient,levels=c("No_Alignment","Inverted_Alignment","Alignment")))%>%
dplyr::mutate(Category=factor(Category,levels=c('Category_1','Category_2','Category_3','Category_4','Category_5'))) %>%
ggplot(aes(x=perc_identity,group=Category,fill=orient))+
geom_histogram(binwidth=5)+
facet_wrap(~Category)
ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig7X.png'
	,plot=plot1 ,width=5,height=3,bg='#ffffff')

plot2B <- SNPD_blast_hits %>% 
select(AW_ID,AW_Comp,Category,perc_identity,align_length,orient) %>% 
dplyr::mutate(perc_identity=case_when(
	is.na(perc_identity)~0,
	TRUE ~ perc_identity)) %>%
dplyr::mutate(orient=case_when(
	is.na(orient)~'No_Alignment',
	TRUE ~ orient)) %>%
dplyr::mutate(orient=factor(orient,levels=c("No_Alignment","Inverted_Alignment","Alignment")))%>%
dplyr::mutate(Category=factor(Category,levels=c('Category_1','Category_2','Category_3','Category_4','Category_5'))) %>%
ggplot(aes(x=align_length,y=perc_identity,group=Category,color=orient))+
geom_point()+
facet_wrap(~Category)
ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig7X2B.png'
	,plot=plot2B ,width=5,height=2.5,bg='#ffffff')


blast_results_files <- list.files('/home/brandvai/mmunasin/TE_Intron/store_data/blast_edges/blast_SNPD/L100_R500_blast_outfiles',full.names=T)
blast_results <- lapply(blast_results_files,fread)

#Subset actual results
SNPD_blast_hits <- blast_results[sapply(blast_results, function(x) dim(x)[1]) > 0]
SNPD_blast_hits <- do.call('rbind',SNPD_blast_hits)
colnames(SNPD_blast_hits) <- c('query_file','subject_file','perc_identity','align_length','num_mistmatches','gap_opens','q_start','q_end','s_start','s_end','evalue','bit_score')
SNPD_blast_hits <- SNPD_blast_hits %>%
dplyr::mutate(orient=case_when(
	q_end > q_start & s_end > s_start ~ 'Alignment',
	q_start > q_end & s_start > s_end ~ 'Alignment',
	q_end > q_start & s_start > s_end ~ 'Inverted_Alignment',
	q_start > q_end & s_end > s_start ~ 'Inverted_Alignment'
	))

SNPD_blast_hits$AW_ID <- unlist(lapply(SNPD_blast_hits$query_file, function(x) paste(unlist(str_split(x,pattern=':'))[4],collapse='_')))
SNPD_blast_hits$AW_Comp <- unlist(lapply(SNPD_blast_hits$query_file, function(x) paste(unlist(str_split(x,pattern=':'))[5],collapse='_')))
unique_id <- paste(SNPD_blast_hits$AW_ID,SNPD_blast_hits$AW_Comp,sep='_')

SNPD_blast_hits <- B73_SNPD_SV_categories %>% select(AW_ID,AW_Comp,Category) %>%
left_join(x=SNPD_blast_hits,y=.,by=c('AW_ID','AW_Comp')) %>%
group_by(AW_ID,AW_Comp) %>%
filter(row_number()==1) %>%
select(-query_file,-subject_file)

LTR_put_insertions <- SNPD_blast_hits %>% 
filter(align_length >= 50) %>% 
filter(orient=='Alignment') %>% 
select(AW_ID,AW_Comp,Category) %>% 
dplyr::mutate(unique_ID=paste(AW_ID,AW_Comp,sep='_'))

######################################################
#EMBOSS ANALYSIS
######################################################
# Extract Fledge Coordinates from SVs
SNPD_X100_left_fledge_bed <- extract_needle_files(B73_SNPD_SV_categories,X=100,tag='left_fledge')
SNPD_X100_right_fledge_bed <- extract_needle_files(B73_SNPD_SV_categories,X=100,tag='right_fledge')
write_fledge_beds(SNPD_X100_left_fledge_bed,X=100,tag='interior',fledge_type='left_fledge')
write_fledge_beds(SNPD_X100_right_fledge_bed,X=100,tag='interior',fledge_type='right_fledge')

SNPD_I100_left_fledge_bed <- extract_needle_files(B73_SNPD_SV_categories,X=100,tag='left_interior')
SNPD_I100_right_fledge_bed <- extract_needle_files(B73_SNPD_SV_categories,X=100,tag='right_interior')
write_fledge_beds(SNPD_I100_left_fledge_bed,X=100,tag='interior',fledge_type='left_interior')
write_fledge_beds(SNPD_I100_right_fledge_bed,X=100,tag='interior',fledge_type='right_interior')

left_fledge_fasta_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/blast_edges/SNPD_I100/SNPD_I100_left_fledges'
right_fledge_fasta_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/blast_edges/SNPD_I100/SNPD_I100_right_fledges'
generate_fledge_param(left_fledge_fasta_dir,right_fledge_fasta_dir,type='needle')
generate_fledge_param(left_fledge_fasta_dir,right_fledge_fasta_dir,type='dottup')

needle_identity <- fread('~/TE_Intron/store_data/blast_edges/SNPD_X100/SNPD_X100_all_identity.txt',sep=':')[,c(1,3)]
colnames(needle_identity) <- c('ID','Identity_Score')
needle_identity <- needle_identity %>% separate(ID,sep='[.]',into=c('AW_ID','ext')) %>%
select(AW_ID,Identity_Score) %>% separate(Identity_Score,sep = " (?=[^ ]+$)",into=c('Identity_Fraction','Percent')) %>%
dplyr::mutate(across(c('Percent'), substr, 2, nchar(Percent))) %>%
dplyr::mutate(across(c('Percent'), substr, 1, nchar(Percent)-2))
