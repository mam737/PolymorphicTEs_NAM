.libPaths('/home/brandvai/mmunasin/Rlibs4')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Rmisc)
library(forcats)
library(ggforce)
library(MetBrewer)

#Functions


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

obtain_polymorphicTE_inSNPD <- function(TE.df,insertion_cat.df) {
	polyTEs.df <- TE.df %>% subset(TE_Classification=='polymorphic')
	colnames(polyTEs.df)[9] <- 'TE_LC'
	polyTEs.bed <- polyTEs.df %>% select(6,7,8,9,10,11,12,13,14,15,16,5)
	colnames(polyTEs.bed)[1:3] <- c('chr','start','end')

	insert_cat.bed <- insertion_cat.df %>% select(6,7,8,9,10,12,15,5)
	colnames(insert_cat.bed)[1:3] <- c('chr','start','end')

	polyTE_insert_cat_SNPD <- list()

	for (i in 1:length(unique(polyTEs.bed$TE_LC))) {
		LC <- unique(polyTEs.bed$TE_LC)[i]
		LC_polymorphicTEs <- polyTEs.bed %>% subset(TE_LC==LC)
		
		LC_insertion_categories <- insert_cat.bed %>% subset(AW_Comp==LC)
	
		polymorphicTE_insertion_overlap <- bedtoolsr::bt.intersect(LC_polymorphicTEs,LC_insertion_categories,wo=T)
	
		colnames(polymorphicTE_insertion_overlap) <- c('TE_chr','TE_start','TE_end','TE_LC','TE_name','TE_class','TE_superfamily','TE_family','TE_method','TE_classification','TE_size','TE_SNPD_ID','AW_chr','AW_start','AW_end','AW_BlockID','AW_LC','AW_Size','AW_Category','AW_SNPD_ID','bp_overlap')
		polyTE_insert_cat_SNPD[[i]] <- polymorphicTE_insertion_overlap
	}
	polymorphicTE_InsertionCategory_SNPD <- do.call('rbind',polyTE_insert_cat_SNPD)
	return(polymorphicTE_InsertionCategory_SNPD)
}

################################################################################################
################################################################################################
############### color palette notes ############################################################
################################################################################################
################################################################################################

# We're going to use colors from the R NatParksPalette palette: DeathValley
# The following three colors refer to AW Block Assignment: alignable(#18315a), structural insertion(#cb6a2d), and unalignable (#68434e)
# The following three colors refer to TE Class Assignment: shared(#416191), polymorphic (#fbb25d), and ambiguous (#b57f84)

# This ordering of the NAM lines follows those in the NAM paper
NAM_level_orders <- c("B73","B97","Ky21","M162W","Ms71","Oh43","Oh7B","M37W","Mo18W","Tx303","HP301","P39","Il14H","CML52","CML69","CML103","CML228","CML247","CML277","CML322","CML333","Ki3","Ki11","NC350","NC358","Tzi8")
chr_level_orders <- c("chr1","chr2","chr3",'chr4','chr5','chr6','chr7','chr8','chr9','chr10')
classification_level_orders <- c("shared",'polymorphic','ambiguous')

################################################################################################
################################################################################################
############### Table S1 #################################################################
################################################################################################
################################################################################################

filtered_TE_anno_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/filtered_TE_anno'

filtered_TE_anno <- obtain_filteredTE_files(filtered_TE_anno_dir)

TE_summary <- filtered_TE_anno %>% group_by(lineage) %>% dplyr::summarise(TotalTE=n(),TotalMb=sum(size)/1000000)

TE_summary <- filtered_TE_anno %>% group_by(lineage,method) %>% 
dplyr::summarise(TotalTE=n(),TotalMb=sum(size)/1000000) %>%
pivot_wider(names_from=method,values_from=c(TotalTE,TotalMb),names_vary='slowest') %>%
left_join(x=TE_summary,by='lineage')

TE_summary <- filtered_TE_anno %>% group_by(lineage,class) %>% 
dplyr::summarise(TotalTE=n(),TotalMb=sum(size)/1000000) %>%
pivot_wider(names_from=class,values_from=c(TotalTE,TotalMb),names_vary='slowest') %>%
left_join(x=TE_summary,by='lineage')

TE_summary <- filtered_TE_anno %>% group_by(lineage,condense_superfamily) %>% 
dplyr::summarise(TotalTE=n(),TotalMb=sum(size)/1000000) %>%
pivot_wider(names_from=condense_superfamily,values_from=c(TotalTE,TotalMb),names_vary='slowest') %>%
left_join(x=TE_summary,by='lineage')

filtered_TE_anno %>% group_by(lineage,class,method)%>% 
dplyr::summarise(TotalTE=n(),TotalMb=sum(size)/1000000) %>%
pivot_wider(names_from=c(class,method),values_from=c(TotalTE,TotalMb),names_vary='slowest') %>%
ungroup() %>% select(-1) %>% colMeans()

fwrite(TE_summary,file='/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SupplementalTables/TableS1.csv')

FamSize_Dist <- filtered_TE_anno %>% 
group_by(lineage,raw_family) %>% 
dplyr::summarise(FamSize=n()) %>% 
ungroup() %>% group_by(raw_family) %>% 
dplyr::summarise(FamSize_Range=max(FamSize) - min(FamSize),
	FamSize_Median=median(FamSize),
	FamSize_RangeRatio=(max(FamSize)-min(FamSize))/median(FamSize),
	FamSize_MaxRatio=max(FamSize)/median(FamSize))

s1 <- ggplot(FamSize_Dist,aes(x=FamSize_Range))+
geom_histogram(binwidth=10,fill='lightblue')+
labs(title='Distribution of Range of Family Sizes',x='Range of Family Size Across NAM',y='Count')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/FamSize_Range_Histogram.png'
	,plot=s1 ,width=4.0,height=2.5)

s2 <- ggplot(FamSize_Dist,aes(x=FamSize_Median,y=FamSize_Range)) +
geom_point(color='lightblue',size=0.5)+
labs(title='Median Family Size vs Range of Family Size',x='Median Family Size',y='Range of Family Size')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/FamSize_Range_Point.png'
	,plot=s2 ,width=4.0,height=2.5)

s3 <- FamSize_Dist %>% subset(FamSize_Median > 100) %>%
ggplot(aes(x=FamSize_RangeRatio)) +
geom_histogram(binwidth=0.02,fill='lightblue')+
labs(title='Range of Fam Size (>100)/Median Fam Size Ratio Distribution',x='Fam Size Range/Median Fam Size',y='Count')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/FamSize_RangeMed_Ratio.png'
	,plot=s3 ,width=4.0,height=2.5)

s4 <- FamSize_Dist %>% subset(FamSize_Median > 100) %>%
ggplot(aes(x=FamSize_MaxRatio)) +
geom_histogram(binwidth=0.02,fill='lightblue')+
labs(title='Max Fam Size (>100)/Median Fam Size Ratio Distribution',x='Max Fam Size/Median Fam Size',y='Count')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/FamSize_MaxMed_Ratio.png'
	,plot=s4 ,width=4.0,height=2.5)

################################################################################################
################################################################################################
############### Figure 1B + S1 #################################################################
################################################################################################
################################################################################################

AW_head_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/summarised_AnchorWave_Regions'

all_AW_data <- obtain_AnchorWave_files(AW_head_dir)

summarised_all_AW_data <- all_AW_data %>% group_by(Lineage_Comp,reduced_type) %>% 
dplyr::summarise(ID_BlockTotal_Mb=sum(ID_BlockSize)/1000000, ASM_BlockTotal_Mb=sum(ASM_BlockSize)/1000000) %>% data.frame()

B73_AW_props <- summarised_all_AW_data %>% group_by(Lineage_Comp) %>% dplyr::mutate(IDBlock_ClassificationProp=ID_BlockTotal_Mb/sum(ID_BlockTotal_Mb)) %>% summarySE(measurevar='IDBlock_ClassificationProp',groupvars=c('reduced_type'))%>% subset(reduced_type != 'structural_insertion_inNAM')
B73_AW_total_means <- summarised_all_AW_data %>% group_by(Lineage_Comp) %>% summarySE(measurevar='ID_BlockTotal_Mb',groupvars=c('reduced_type'))%>% subset(reduced_type != 'structural_insertion_inNAM')

fig_1b_1 <- ggplot(B73_AW_total_means,aes(x=reduced_type,y=ID_BlockTotal_Mb,fill=reduced_type)) +
geom_bar(stat='identity',width=0.5) + 
geom_errorbar(aes(ymin=ID_BlockTotal_Mb-sd,ymax=ID_BlockTotal_Mb+sd),width=0.25,position=position_dodge(.9))+
scale_fill_manual(labels = c("Alignable", "Structural Insertion", "Unalignable"),values=c("#18315a",'#cb6a2d','#68434e')) +
theme_bw()+ ylim(c(0,1500))+
theme(legend.position = "none",axis.text.x=element_blank())+
scale_x_discrete(labels=c('Alignable', 'SV Sequence in B73', 'Unalignable'))+
labs(x='',y='Mb in Genome') + 
guides(fill = guide_legend(title = "AW Block Assignment")) + 
ggtitle("")

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig1b.pdf'
	,plot=fig_1b_1 ,width=1.5,height=2.0)

average_all_AW_data <- summarised_all_AW_data %>% group_by(reduced_type) %>% dplyr::summarise(ID_BlockTotal_Mb=mean(ID_BlockTotal_Mb),ASM_BlockTotal_Mb=mean(ASM_BlockTotal_Mb)) %>% dplyr::mutate(Lineage_Comp='Average') %>% select(4,1,2,3)

all_average_AW_data <- rbind(summarised_all_AW_data,average_all_AW_data)

all_average_AW_data$Lineage_Comp <- factor(all_average_AW_data$Lineage_Comp, levels=c(NAM_level_orders,'Average'))

figs1_1 <- all_average_AW_data %>% subset(reduced_type != 'structural_insertion_inNAM') %>%
ggplot(aes(x=Lineage_Comp,y=ID_BlockTotal_Mb,fill=reduced_type)) + 
geom_bar(position="fill", stat="identity")+ 
scale_fill_manual(labels = c("Alignable", "SV Sequence", "Unalignable"),values=c("#18315a",'#cb6a2d','#68434e')) +
theme(legend.position = "bottom",axis.text.x = element_text(angle = 90)) +
labs(x='NAM Line Comparison',y='Mb Prop in Genome') + 
guides(fill = guide_legend(title = "AW Block Assignment")) + 
ggtitle("B73 Total AnchorWave Classification By Type")

figs1_2 <- all_average_AW_data %>% subset(reduced_type != 'structural_insertion_inB73') %>%
ggplot(aes(x=Lineage_Comp,y=ASM_BlockTotal_Mb,fill=reduced_type)) + 
geom_bar(position="fill", stat="identity")+ 
scale_fill_manual(labels = c("Alignable", "SV Sequence", "Unalignable"),values=c("#18315a",'#cb6a2d','#68434e')) +
theme(legend.position = "bottom",axis.text.x = element_text(angle = 90)) +
labs(x='NAM Line Comparison',y='Mb Prop in Genome') + 
guides(fill = guide_legend(title = "AW Block Assignment")) + 
ggtitle("NAM Total AnchorWave Classification By Type")

figs1 <- ggarrange(figs1_1,figs1_2,ncol=1,nrow=2,common.legend=T,legend='bottom')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS1.pdf'
	,plot=figs1,width=8,height=6)

################################################################################################
################################################################################################
############### Figure 1C + S2 #################################################################
################################################################################################
################################################################################################

B73_TE_n_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/polymorphic_TE_calls/B73/TE_n'
NAM_TE_n_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/polymorphic_TE_calls/NAM/TE_n'

all_B73_TE_n_calls <- process_TE_calls(B73_TE_n_dir)
all_NAM_TE_n_calls <- process_TE_calls(NAM_TE_n_dir)

summarised_TE_B73 <- all_B73_TE_n_calls %>% 
group_by(ASM_Comp,classification) %>% 
dplyr::summarise(TE_ClassificationTotal=n()) %>% 
dplyr::mutate(tag='B73')

summarised_TE_NAM <- all_NAM_TE_n_calls %>% 
group_by(ASM_Comp,classification) %>% 
dplyr::summarise(TE_ClassificationTotal=n()) %>% 
dplyr::mutate(tag='NAM')

TE_totals <- summarySE(rbind(summarised_TE_B73,summarised_TE_NAM), measurevar="TE_ClassificationTotal", groupvars=c("tag","classification"))

B73_TE_props <- summarised_TE_B73 %>% 
group_by(ASM_Comp) %>% 
dplyr::mutate(TE_ClassificationProp=TE_ClassificationTotal/sum(TE_ClassificationTotal)) %>% 
summarySE(measurevar="TE_ClassificationProp", groupvars=c("classification")) %>% 
dplyr::mutate(tag='B73') %>%
dplyr::mutate(classification=factor(classification,levels=classification_level_orders))

fig1c_r1c1 <- ggplot(B73_TE_props,aes(x=tag,y=TE_ClassificationProp,fill=classification)) +
geom_bar(position='dodge',stat='identity',width=0.75) +
geom_errorbar(aes(ymin=TE_ClassificationProp-sd,ymax=TE_ClassificationProp+sd),width=0.25,position=position_dodge(.75))+
scale_fill_manual(labels = c("Shared","Polymorphic","Ambiguous"),values=c("#416191",'#fbb25d','#b57f84')) +
scale_x_discrete(labels=c('B73', 'NAM'))+ coord_cartesian(ylim = c(0, 1)) +
labs(x='',y='Mean Proportion') +
theme_bw()+ 
guides(fill = guide_legend(title = "Classification"))

B73_exon_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/polymorphic_gene_calls/B73/exon_calls'
B73_full_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/polymorphic_gene_calls/B73/full_calls'
NAM_exon_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/polymorphic_gene_calls/NAM/exon_calls'
NAM_full_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/polymorphic_gene_calls/NAM/full_calls'

all_B73_exon_calls <- process_gene_calls(B73_exon_dir,tag='id',type='exon')
all_B73_full_calls <- process_gene_calls(B73_full_dir,tag='id',type='gene')
all_NAM_exon_calls <- process_gene_calls(NAM_exon_dir,tag='asm',type='exon')
all_NAM_full_calls <- process_gene_calls(NAM_full_dir,tag='asm',type='gene')

syntenic_gene_list_dir <- '/home/springer/shared/ns-genome/Zmays_NAM_genomes/NAM_syntenic_genes'

syntenic_genes_all <- obtain_gene_synteny(syntenic_gene_list_dir)

all_B73_exon_calls <- add_synteny_column(all_B73_exon_calls,syntenic_genes_all)
all_B73_full_calls <- add_synteny_column(all_B73_full_calls,syntenic_genes_all)
all_NAM_exon_calls <- add_synteny_column(all_NAM_exon_calls,syntenic_genes_all)
all_NAM_full_calls <- add_synteny_column(all_NAM_full_calls,syntenic_genes_all)

summarised_nsg_B73_full <- obtain_synteny_summary(all_B73_full_calls,synteny_tag=FALSE,type_tag='id')
summarised_nsg_NAM_full <- obtain_synteny_summary(all_NAM_full_calls,synteny_tag=FALSE,type_tag='asm')

nsg_full_totals <- summarySE(rbind(summarised_nsg_B73_full,summarised_nsg_NAM_full), measurevar="gene_ClassificationTotal", groupvars=c("tag","classification"))

nsg_full_props <- summarised_nsg_B73_full %>% 
group_by(ASM_Comp) %>% 
dplyr::mutate(gene_ClassificationProp=gene_ClassificationTotal/sum(gene_ClassificationTotal)) %>% 
summarySE(measurevar="gene_ClassificationProp", groupvars=c("classification")) %>%
dplyr::mutate(tag='B73')

fig1c_r2c2 <- ggplot(nsg_full_props,aes(x=tag,y=gene_ClassificationProp,fill=classification)) +
geom_bar(position='dodge',stat='identity',width=0.75) +
geom_errorbar(aes(ymin=gene_ClassificationProp-sd,ymax=gene_ClassificationProp+sd),width=0.25,position=position_dodge(.75))+
scale_fill_manual(labels = c("Shared","Polymorphic","Ambiguous"),values=c("#416191",'#fbb25d','#b57f84')) +
scale_x_discrete(labels=c('B73'))+ coord_cartesian(ylim = c(0, 1)) +
labs(x='',y='Mean Proportion') + 
theme_bw()+
guides(fill = guide_legend(title = "Classification"))

summarised_sg_B73_full <- obtain_synteny_summary(all_B73_full_calls,synteny_tag=TRUE,type_tag='id')
summarised_sg_NAM_full <- obtain_synteny_summary(all_NAM_full_calls,synteny_tag=TRUE,type_tag='asm')

sg_full_totals <- summarySE(rbind(summarised_sg_B73_full,summarised_sg_NAM_full), measurevar="gene_ClassificationTotal", groupvars=c("tag","classification"))

sg_full_props <- summarised_sg_B73_full %>% group_by(ASM_Comp) %>% dplyr::mutate(gene_ClassificationProp=gene_ClassificationTotal/sum(gene_ClassificationTotal)) %>% summarySE(measurevar="gene_ClassificationProp", groupvars=c("classification")) %>% mutate(tag='B73')

fig1c_r2c3 <- ggplot(sg_full_props,aes(x=tag,y=gene_ClassificationProp,fill=classification)) +
geom_bar(position='dodge',stat='identity',width=0.75) +
geom_errorbar(aes(ymin=gene_ClassificationProp-sd,ymax=gene_ClassificationProp+sd),width=0.25,position=position_dodge(.75))+
scale_fill_manual(labels = c("Shared", "Polymorphic", "Ambiguous"),values=c("#416191",'#fbb25d','#b57f84')) +
scale_x_discrete(labels=c('B73'))+ coord_cartesian(ylim = c(0, 1)) +
labs(x='',y='Mean Proportion') + 
theme_bw()+
guides(fill = guide_legend(title = "Classification"))

fig1cde <- ggarrange(fig1c_r1c1,fig1c_r2c3,fig1c_r2c2,ncol=3,nrow=1,common.legend=T,legend='none')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig1cde.pdf'
	,plot=fig1cde,width=4,height=2.0,bg='#ffffff')

average_summarised_TE_B73 <- summarised_TE_B73 %>% group_by(classification) %>% dplyr::summarise(TE_ClassificationTotal=mean(TE_ClassificationTotal)) %>% dplyr::mutate(ASM_Comp='Average',tag='B73') %>% select(3,1,2,4)
all_average_summarised_TE_B73 <- rbind(summarised_TE_B73,average_summarised_TE_B73)
all_average_summarised_TE_B73$ASM_Comp <- factor(all_average_summarised_TE_B73$ASM_Comp, levels=c(NAM_level_orders,'Average'))
all_average_summarised_TE_B73$classification <- factor(all_average_summarised_TE_B73$classification,levels=classification_level_orders)

figs2_r1c1 <- all_average_summarised_TE_B73 %>%
ggplot(aes(x=ASM_Comp,y=TE_ClassificationTotal,fill=classification)) + 
geom_bar(position="fill", stat="identity")+ 
scale_fill_manual(labels = c("Shared","Polymorphic","Ambiguous"),values=c("#416191",'#fbb25d','#b57f84')) +
theme(legend.position = "bottom",axis.text.x = element_text(angle = 90)) +
labs(x='NAM Line Comparison',y='TE Proportions') + 
guides(fill = guide_legend(title = "Classification State")) + 
ggtitle("B73 TE Classification Totals")

average_summarised_TE_NAM <- summarised_TE_NAM %>% group_by(classification) %>% dplyr::summarise(TE_ClassificationTotal=mean(TE_ClassificationTotal)) %>% dplyr::mutate(ASM_Comp='Average',tag='NAM') %>% select(3,1,2,4)
all_average_summarised_TE_NAM <- rbind(summarised_TE_NAM,average_summarised_TE_NAM)
all_average_summarised_TE_NAM$ASM_Comp <- factor(all_average_summarised_TE_NAM$ASM_Comp, levels=c(NAM_level_orders,'Average'))
all_average_summarised_TE_NAM$classification <- factor(all_average_summarised_TE_NAM$classification,levels=classification_level_orders)

figs2_r1c2 <- all_average_summarised_TE_NAM %>%
ggplot(aes(x=ASM_Comp,y=TE_ClassificationTotal,fill=classification)) + 
geom_bar(position="fill", stat="identity")+ 
scale_fill_manual(labels = c("Shared","Polymorphic","Ambiguous"),values=c("#416191",'#fbb25d','#b57f84')) +
theme(legend.position = "bottom",axis.text.x = element_text(angle = 90)) +
labs(x='NAM Line Comparison',y='TE Proportions') + 
guides(fill = guide_legend(title = "Classification State")) + 
ggtitle("NAM TE Classification Totals")

summarised_B73_exon_calls <- all_B73_exon_calls %>% 
group_by(ASM_Comp,classification) %>%
dplyr::summarise(gene_ClassificationTotal=n()) %>% 
dplyr::mutate(tag='B73')

average_summarised_B73_exon_calls <- summarised_B73_exon_calls %>% group_by(classification) %>% dplyr::summarise(gene_ClassificationTotal=mean(gene_ClassificationTotal)) %>% dplyr::mutate(ASM_Comp='Average',tag='B73') %>% select(3,1,2,4)
all_average_summarised_B73_exon_calls <- rbind(summarised_B73_exon_calls,average_summarised_B73_exon_calls)
all_average_summarised_B73_exon_calls$ASM_Comp <- factor(all_average_summarised_B73_exon_calls$ASM_Comp, levels=c(NAM_level_orders,'Average'))

figs2_r2c1 <- all_average_summarised_B73_exon_calls %>%
ggplot(aes(x=ASM_Comp,y=gene_ClassificationTotal,fill=classification)) + 
geom_bar(position="fill", stat="identity")+ 
scale_fill_manual(labels = c("Shared","Polymorphic","Ambiguous"),values=c("#416191",'#fbb25d','#b57f84')) +
theme(legend.position = "bottom",axis.text.x = element_text(angle = 90)) +
labs(x='NAM Line Comparison',y='Gene Proportions') + 
guides(fill = guide_legend(title = "Classification State")) + 
ggtitle("B73 Exon-Only Gene Classification Totals")

summarised_NAM_exon_calls <- all_NAM_exon_calls %>% 
group_by(ASM_Comp,classification) %>%
dplyr::summarise(gene_ClassificationTotal=n()) %>% 
dplyr::mutate(tag='NAM')

average_summarised_NAM_exon_calls <- summarised_NAM_exon_calls %>% group_by(classification) %>% dplyr::summarise(gene_ClassificationTotal=mean(gene_ClassificationTotal)) %>% dplyr::mutate(ASM_Comp='Average',tag='NAM') %>% select(3,1,2,4)
all_average_summarised_NAM_exon_calls <- rbind(summarised_NAM_exon_calls,average_summarised_NAM_exon_calls)
all_average_summarised_NAM_exon_calls$ASM_Comp <- factor(all_average_summarised_NAM_exon_calls$ASM_Comp, levels=c(NAM_level_orders,'Average'))

figs2_r2c2 <- all_average_summarised_NAM_exon_calls %>%
ggplot(aes(x=ASM_Comp,y=gene_ClassificationTotal,fill=classification)) + 
geom_bar(position="fill", stat="identity")+ 
scale_fill_manual(labels = c("Shared","Polymorphic","Ambiguous"),values=c("#416191",'#fbb25d','#b57f84')) +
theme(legend.position = "bottom",axis.text.x = element_text(angle = 90)) +
labs(x='NAM Line Comparison',y='Gene Proportions') + 
guides(fill = guide_legend(title = "Classification State")) + 
ggtitle("NAM Exon-Only Gene Classification Totals")

summarised_B73_full_calls <- all_B73_full_calls %>% 
group_by(ASM_Comp,classification) %>%
dplyr::summarise(gene_ClassificationTotal=n()) %>% 
dplyr::mutate(tag='B73')

average_summarised_B73_full_calls <- summarised_B73_full_calls %>% 
group_by(classification) %>% 
dplyr::summarise(gene_ClassificationTotal=mean(gene_ClassificationTotal)) %>% 
dplyr::mutate(ASM_Comp='Average',tag='B73') %>% select(3,1,2,4)

all_average_summarised_B73_full_calls <- rbind(summarised_B73_full_calls,average_summarised_B73_full_calls)
all_average_summarised_B73_full_calls$ASM_Comp <- factor(all_average_summarised_B73_full_calls$ASM_Comp, levels=c(NAM_level_orders,'Average'))

figs2_r3c1 <- all_average_summarised_B73_full_calls %>%
ggplot(aes(x=ASM_Comp,y=gene_ClassificationTotal,fill=classification)) + 
geom_bar(position="fill", stat="identity")+ 
scale_fill_manual(labels = c("Shared","Polymorphic","Ambiguous"),values=c("#416191",'#fbb25d','#b57f84')) +
theme(legend.position = "bottom",axis.text.x = element_text(angle = 90)) +
labs(x='NAM Line Comparison',y='Gene Proportions') + 
guides(fill = guide_legend(title = "Classification State")) + 
ggtitle("B73 Full-Length Gene Classification Totals")

summarised_NAM_full_calls <- all_NAM_full_calls %>% 
group_by(ASM_Comp,classification) %>%
dplyr::summarise(gene_ClassificationTotal=n()) %>% 
dplyr::mutate(tag='NAM')

average_summarised_NAM_full_calls <- summarised_NAM_full_calls %>% 
group_by(classification) %>% 
dplyr::summarise(gene_ClassificationTotal=mean(gene_ClassificationTotal)) %>% 
dplyr::mutate(ASM_Comp='Average',tag='NAM') %>% select(3,1,2,4)

all_average_summarised_NAM_full_calls <- rbind(summarised_NAM_full_calls,average_summarised_NAM_full_calls)
all_average_summarised_NAM_full_calls$ASM_Comp <- factor(all_average_summarised_NAM_full_calls$ASM_Comp, levels=c(NAM_level_orders,'Average'))

figs2_r3c2 <- all_average_summarised_NAM_full_calls %>%
ggplot(aes(x=ASM_Comp,y=gene_ClassificationTotal,fill=classification)) + 
geom_bar(position="fill", stat="identity")+ 
scale_fill_manual(labels = c("Shared","Polymorphic","Ambiguous"),values=c("#416191",'#fbb25d','#b57f84')) +
theme(legend.position = "bottom",axis.text.x = element_text(angle = 90)) +
labs(x='NAM Line Comparison',y='Gene Proportions') + 
guides(fill = guide_legend(title = "Classification State")) + 
ggtitle("NAM Full-Length Gene Classification Totals")

figs2 <- ggarrange(figs2_r1c1,NULL,figs2_r1c2,figs2_r2c1,NULL,figs2_r2c2,figs2_r3c1,NULL,figs2_r3c2,nrow=3,ncol=3,widths=c(1,0.15,1),common.legend=T,legend='bottom')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS2.pdf'
	,plot=figs2,width=10,height=6)

################################################################################################
################################################################################################
############### Figure 2 #################################################################
################################################################################################
################################################################################################

B73_TE_call_dist <- all_B73_TE_n_calls %>% 
group_by(TE_name,classification,.drop=F) %>% 
dplyr::summarise(n=n()) %>% subset(classification=='shared')

B73_TE_call_dist <- shared_frequency(B73_TE_call_dist)

# Old 2A Color Palette - c('#212D51','#B7ABBC','#FEB424','#D8511D')

fig2A <- B73_TE_call_dist %>% group_by(n,degree_shared) %>% dplyr::summarise(count_n=n()) %>%
ggplot(aes(x=n,y=count_n,fill=degree_shared)) +
geom_bar(stat='identity') + 
scale_fill_manual(labels=c('Private','Variable','Near Core','Core'),values=c("#d1b252","#a97f2f","#7e5522","#472c0b"))+
labs(x='Times Classified as Shared Between B73 and NAM',y='Number of TEs') + 
theme_bw()+
theme(legend.position = "bottom",legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size=7),axis.text.x = element_text(angle = 90,hjust=1.10,vjust=0.95))+
guides(fill=guide_legend(title=''))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig2A.tif'
	,plot=fig2A,width=5,height=2,device='tiff',dpi=700,bg='#ffffff')

all_B73_full_calls$syntenic_gene <- factor(all_B73_full_calls$syntenic_gene,levels=c('TRUE','FALSE'))

B73_gene_call_dist <- all_B73_full_calls %>% 
group_by(gene_name,classification,.drop=F) %>% 
dplyr::summarise(n=n()) %>% 
subset(classification=='shared')

B73_gene_call_dist <- shared_frequency(B73_gene_call_dist)

structural_B73_TEs <- all_B73_TE_n_calls %>% 
subset(method=='structural') %>% pull(TE_name) %>% unique()

method_B73_TE_summary <- B73_TE_call_dist %>% 
mutate(method=case_when(
	TE_name %in% structural_B73_TEs ~ 'Structural', TRUE ~ 'Homology')) %>%
group_by(method,degree_shared) %>%
dplyr::summarise(total_count=n()) %>% select(2,3,1)
colnames(method_B73_TE_summary)[3] <- 'Grouping'

syntenic_B73_gene_summary <- B73_gene_call_dist %>% 
mutate(syntenic_gene=case_when(
	gene_name %in% syntenic_genes_all$NAM_gene_name ~ 'TRUE',TRUE ~ 'FALSE')) %>%
group_by(syntenic_gene,degree_shared) %>% dplyr::summarise(total_count=n()) %>% mutate(synteny_call=case_when(syntenic_gene=='TRUE' ~ 'Syntenic',TRUE ~ 'Non-Syntenic')) %>% ungroup() %>% select(-syntenic_gene)
colnames(syntenic_B73_gene_summary)[3] <- 'Grouping'

B73_summary <- rbind(method_B73_TE_summary,syntenic_B73_gene_summary)
B73_summary$Grouping <- factor(B73_summary$Grouping,levels=c('Non-Syntenic','Syntenic','Homology','Structural'))

fig2C <- ggplot(B73_summary,aes(x=Grouping,y=total_count,fill=fct_rev(degree_shared)))+
geom_bar(position='fill',stat='identity',width=0.4)+
scale_fill_manual(labels=rev(c('Private','Variable','Near Core','Core')),values=rev(c("#d1b252","#a97f2f","#7e5522","#472c0b")))+
labs(x='',y='Total %') + coord_flip()+
theme_bw()+
theme(legend.position = "none")+
guides(fill=guide_legend(title=''))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig2C.png'
	,plot=fig2C,width=5,height=1.5,bg='#ffffff')
ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig2C.tif'
	,plot=fig2C,width=4,height=2,device='tiff',dpi=700,bg='#ffffff')

structuralTE_identity <- fread('/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/structurally_annotated_LTRs_identity.tsv')
B73_structuralTE_identity <- structuralTE_identity %>% subset(lineage=='B73')
B73_identity_score <- B73_structuralTE_identity %>% select(name,ltr_identity)
colnames(B73_identity_score)[1] <- 'TE_name'

structuralLTR_B73_TE_call_dist <- B73_TE_call_dist %>% 
subset(TE_name %in% B73_structuralTE_identity$name)
structuralLTR_B73_TE_call_dist <- left_join(structuralLTR_B73_TE_call_dist,B73_identity_score,by='TE_name')
structuralLTR_B73_TE_call_dist <- structuralLTR_B73_TE_call_dist %>%
dplyr::mutate(ltr_age=case_when(
	ltr_identity==1 ~ 'Very Young',
	ltr_identity >= 0.95 & ltr_identity < 1 ~ 'Young',
	ltr_identity >= 0.9 & ltr_identity < 0.95 ~ 'Moderate',
	ltr_identity >= 0.8 & ltr_identity < 0.9 ~ 'Old')) %>%
dplyr::mutate(ltr_age=factor(ltr_age,levels=c('Very Young','Young','Moderate','Old')))

fig2D_alt2 <- structuralLTR_B73_TE_call_dist %>% 
group_by(degree_shared,ltr_age) %>% dplyr::summarise(n=n()) %>%
ggplot(aes(x=degree_shared,y=n,fill=ltr_age))+
geom_bar(position='fill',stat='identity',width=0.4)+
scale_fill_manual(values=rev(met.brewer("Hokusai2",n=4)))+
labs(x='',y='Total %') + coord_flip()+
theme_bw()+
theme(legend.position = "bottom")+
guides(fill=guide_legend(title='',reverse=TRUE))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig2D_alt2.png'
	,plot=fig2D_alt2,width=4,height=2,bg='#ffffff')
ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig2D_alt2.tif'
	,plot=fig2D_alt2,width=4,height=2,device='tiff',dpi=700,bg='#ffffff')

################################################################################################
################################################################################################
############### Figure 3 #################################################################
################################################################################################
################################################################################################

B73_polymorphic_percent <- obtain_polymorphic_percent(all_B73_TE_n_calls,tag.string='B73')
NAM_polymorphic_percent <- obtain_polymorphic_percent(all_NAM_TE_n_calls,tag.string='NAM')

polymorphic_percent <- rbind(B73_polymorphic_percent,NAM_polymorphic_percent)
polymorphic_percent$method <- factor(polymorphic_percent$method,levels=c('structural','homology'))
levels(polymorphic_percent$method)[levels(polymorphic_percent$method)=="structural"] <- "Structural"
levels(polymorphic_percent$method)[levels(polymorphic_percent$method)=="homology"] <- "Homology"

fig3A <- polymorphic_percent %>% subset(tag=='B73') %>%
ggplot(aes(x=condense_superfamily,y=percent_total,fill=method))+
geom_col(position=position_dodge2(width = 0.75, preserve = "single")) +
geom_errorbar(aes(ymin=percent_total-sd,ymax=percent_total+sd),width=0.25,position=position_dodge(.75))+
scale_fill_manual(values=c("#FBB25D","#E38E45"))+
labs(x='Superfamily',y='% Polymorphic')+
theme_bw()+
guides(fill = guide_legend(title = "TE Set"))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig3A.png'
	,plot=fig3A,width=6.0,height=1.5,bg='#ffffff')

B73_largefam_polymorphic_percent <- obtain_large_polymorphic_percent(all_B73_TE_n_calls,tag.string='B73')
NAM_largefam_polymorphic_percent <- obtain_large_polymorphic_percent(all_NAM_TE_n_calls,tag.string='NAM')

largefam_polymorphic_percent <- rbind(B73_largefam_polymorphic_percent,NAM_largefam_polymorphic_percent)

fig3B_alt1 <- B73_largefam_polymorphic_percent %>% group_by(raw_family,fam_size) %>%
dplyr::summarise(range_perc_total=max(percent_total)-min(percent_total)) %>% 
ggplot(aes(x=range_perc_total)) + geom_histogram(binwidth=5,fill='#e8b697' ,color='#CB6A2D')+
labs(x='% Polymorphic Range',y='Number of B73 TE Families')+
theme_bw()

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig3B_alt1.png'
	,plot=fig3B_alt1,width=4.0,height=3.5,bg='#ffffff')

fig3B_alt2 <- B73_largefam_polymorphic_percent %>% group_by(raw_family,fam_size) %>%
dplyr::summarise(range_perc_total=max(percent_total)-min(percent_total)) %>% 
ggplot(aes(x=fam_size,y=range_perc_total)) + geom_point(size=0.5,color='#CB6A2D')+
labs(x='TE Family Size',y='Range of Polymorphic %')+
theme_bw()

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig3B_alt2.png'
	,plot=fig3B_alt2,width=4.0,height=3.5,bg='#ffffff')

################################################################################################
################################################################################################
############### Figure S3 #################################################################
################################################################################################
################################################################################################

B73_superfam_percents <- all_B73_TE_n_calls %>%
group_by(ASM_Comp,classification,condense_superfamily) %>%
dplyr::summarise(n=n()) %>% ungroup() %>%
group_by(ASM_Comp,condense_superfamily) %>%
dplyr::mutate(percent_total=n/sum(n)*100) %>%
summarySE(measurevar="percent_total", groupvars=c("condense_superfamily","classification"))
colnames(B73_superfam_percents)[1] <- 'grouping'

B73_method_percents <- all_B73_TE_n_calls %>%
group_by(ASM_Comp,classification,method) %>%
dplyr::summarise(n=n()) %>% ungroup() %>%
group_by(ASM_Comp,method) %>%
dplyr::mutate(percent_total=n/sum(n)*100) %>%
summarySE(measurevar="percent_total", groupvars=c("method","classification")) %>%
dplyr::mutate(method=case_when(method=='homology' ~ 'Homology',TRUE ~ 'Structural'))
colnames(B73_method_percents)[1] <- 'grouping'

grouping_order <- c('Structural','Homology','DTA','DTC','DTH','DTM','DTT','RIL','RLC','RLG','RLX')
B73_grouping_percents <- rbind(B73_superfam_percents,B73_method_percents) %>%
dplyr::mutate(grouping=factor(grouping,levels=grouping_order)) %>%
dplyr::mutate(classification=case_when(classification=='shared' ~ 'Shared',
	classification=='ambiguous' ~ 'Ambiguous',
	TRUE ~ 'Polymorphic')) %>%
dplyr::mutate(classification=factor(classification,levels=c('Shared','Ambiguous','Polymorphic')))

figure_s3r1 <- ggplot(B73_grouping_percents,aes(x=grouping,y=percent_total,fill=classification))+
geom_bar(position='dodge',stat='identity',width=0.85)+
geom_errorbar(aes(ymin=percent_total-sd,ymax=percent_total+sd),width=0.25,position=position_dodge(0.75))+
facet_wrap(~classification,ncol=1,nrow=3,scales='free_y')+
labs(x='Group',y='% of Group',title='B73 TE Percent Of Classification By Group')+
scale_fill_manual(name='Classification',labels = c("Shared", "Ambiguous", "Polymorphic"),values=c("#416191",'#b57f84','#fbb25d'))+
theme_bw()+
theme(legend.position='bottom')

NAM_superfam_percents <- all_NAM_TE_n_calls %>%
group_by(ASM_Comp,classification,condense_superfamily) %>%
dplyr::summarise(n=n()) %>% ungroup() %>%
group_by(ASM_Comp,condense_superfamily) %>%
dplyr::mutate(percent_total=n/sum(n)*100) %>%
summarySE(measurevar="percent_total", groupvars=c("condense_superfamily","classification"))
colnames(NAM_superfam_percents)[1] <- 'grouping'

NAM_method_percents <- all_NAM_TE_n_calls %>%
group_by(ASM_Comp,classification,method) %>%
dplyr::summarise(n=n()) %>% ungroup() %>%
group_by(ASM_Comp,method) %>%
dplyr::mutate(percent_total=n/sum(n)*100) %>%
summarySE(measurevar="percent_total", groupvars=c("method","classification")) %>%
dplyr::mutate(method=case_when(method=='homology' ~ 'Homology',TRUE ~ 'Structural'))
colnames(NAM_method_percents)[1] <- 'grouping'

grouping_order <- c('Structural','Homology','DTA','DTC','DTH','DTM','DTT','RIL','RLC','RLG','RLX')
NAM_grouping_percents <- rbind(NAM_superfam_percents,NAM_method_percents) %>%
dplyr::mutate(grouping=factor(grouping,levels=grouping_order)) %>%
dplyr::mutate(classification=case_when(classification=='shared' ~ 'Shared',
	classification=='ambiguous' ~ 'Ambiguous',
	TRUE ~ 'Polymorphic')) %>%
dplyr::mutate(classification=factor(classification,levels=c('Shared','Ambiguous','Polymorphic')))

figure_s3r2 <- ggplot(NAM_grouping_percents,aes(x=grouping,y=percent_total,fill=classification))+
geom_bar(position='dodge',stat='identity',width=0.85)+
geom_errorbar(aes(ymin=percent_total-sd,ymax=percent_total+sd),width=0.25,position=position_dodge(0.75))+
facet_wrap(~classification,ncol=1,nrow=3,scales='free_y')+
labs(x='Group',y='% of Group',title='NAM TE Percent Of Classification By Group')+
scale_fill_manual(name='Classification',labels = c("Shared", "Ambiguous", "Polymorphic"),values=c("#416191",'#b57f84','#fbb25d'))+
theme_bw()+
theme(legend.position='bottom')

figure_s3 <- ggarrange(figure_s3r1,figure_s3r2,nrow=2,ncol=1,common.legend=TRUE,legend='bottom')
ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS3.png'
	,plot=figure_s3,width=7.0,height=6.5,bg='#ffffff')

################################################################################################
################################################################################################
############### Figure S4 #################################################################
################################################################################################
################################################################################################

B73_TE_superfam_calls <- all_B73_TE_n_calls %>% select(TE_name,condense_superfamily) %>% distinct()

B73_TE_call_dist_superfam <- left_join(B73_TE_call_dist,B73_TE_superfam_calls,by='TE_name')

B73_TE_superfam_summary <- B73_TE_call_dist_superfam %>%
group_by(condense_superfamily,degree_shared) %>%
dplyr::summarise(total_count=n()) %>%
dplyr::mutate(condense_superfamily=factor(condense_superfamily,levels=c("DTA","DTC",'DTH','DTM','DTT','RIL','RLC','RLG',"RLX")))

figure_S4 <- ggplot(B73_TE_superfam_summary,aes(x=fct_rev(condense_superfamily),y=total_count,fill=fct_rev(degree_shared))) +
geom_bar(position='fill',stat='identity',width=0.4)+
scale_fill_manual(labels=rev(c('Private','Variable','Near Core','Core')),values=rev(c('#212D51','#B7ABBC','#FEB424','#D8511D')))+
labs(x='',y='Total %',title='B73 TE Groupings By Superfamily Across NAM Comparisons') + coord_flip()+
theme_bw()+
theme(legend.position = "bottom",plot.title = element_text(size=10))+
guides(fill=guide_legend(title=''))


ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS4.png'
	,plot=figure_S4,width=5,height=2.5,bg='#ffffff')

################################################################################################
################################################################################################
############### Figure 4 #################################################################
################################################################################################
################################################################################################

AWBlocks_CategoryAssignment <- fread("/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/all_AWBlocks_CategoryAssignment.tsv")

B73_AWBlocks_CA <- AWBlocks_CategoryAssignment %>% subset(Block_Type=='structural_insertion_inB73')
NAM_AWBlocks_CA <- AWBlocks_CategoryAssignment %>% subset(!(Block_Type %in% c('structural_insertion_inB73','alignable_region','unalignable','Missing_Data')))

summary_B73_AWBlocks_CA <- B73_AWBlocks_CA %>% group_by(Category) %>%
dplyr::summarise(n=n(),Mb=sum(ID_BlockSize)/1000000)

summary_NAM_AWBlocks_CA <- NAM_AWBlocks_CA %>% group_by(Category) %>%
dplyr::summarise(n=n(),Mb=sum(ASM_BlockSize)/1000000)

B73_summary_AW_CA <- summarise_AW_CA(summary_B73_AWBlocks_CA,'B73') 
NAM_summary_AW_CA <- summarise_AW_CA(summary_NAM_AWBlocks_CA,'NAM')

summary_AW_CA <- rbind(B73_summary_AW_CA,NAM_summary_AW_CA)

fig_4C <- ggplot(B73_summary_AW_CA,aes(x=Type,y=Amount,fill=Category)) +
geom_bar(position='fill',stat='identity',width=0.4)+
labs(x='Structural Insertions',y='Proportion',fill='Category')+
scale_fill_manual(values=met.brewer("Hokusai3", 5))+
theme_bw()+
theme(legend.position = "none")

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig4C.png'
	,plot=fig_4C ,width=2,height=4,bg='#ffffff')

################################################################################################
################################################################################################
############### Figure 5 #################################################################
################################################################################################
################################################################################################

B73_genome_cumulMb <- summarised_all_AW_data %>% 
subset(reduced_type != 'structural_insertion_inNAM') %>% 
group_by(reduced_type) %>% 
dplyr::summarise(ID_cumulMb=sum(ID_BlockTotal_Mb))
B73_total_cumulMb <- sum(B73_genome_cumulMb$ID_cumulMb)
B73_insertion_cumulMb <- B73_genome_cumulMb$ID_cumulMb[2]

NAM_genome_cumulMb <- summarised_all_AW_data %>% 
subset(reduced_type != 'structural_insertion_inB73') %>% 
group_by(reduced_type) %>% 
dplyr::summarise(ASM_cumulMb=sum(ASM_BlockTotal_Mb))
NAM_total_cumulMb <- sum(NAM_genome_cumulMb$ASM_cumulMb)
NAM_insertion_cumulMb <- NAM_genome_cumulMb$ASM_cumulMb[2]

total_cumulMb <- B73_total_cumulMb+NAM_total_cumulMb
insertion_cumulMb <- B73_insertion_cumulMb+NAM_insertion_cumulMb

cumulMb <- data.frame(group='Genome-Wide',perc_inserted=(insertion_cumulMb/total_cumulMb)*100,Mb_inserted=insertion_cumulMb)
 
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

cumulMb <- rbind(cumulMb, data.frame(group='SNP Depleted',perc_inserted=(SNPD_insertion_cumulMb/SNPD_cumulMb)*100,Mb_inserted=SNPD_insertion_cumulMb))

Fig5A <- ggplot(cumulMb,aes(x=group,y=perc_inserted)) +
geom_bar(aes(fill = group), position = "dodge", stat="identity",width=0.4)+
scale_fill_manual(labels = c("Genome-Wide","SNP Depleted"),values=c("#2b4655","#738e8e")) +
labs(x='',y='% Cumulative Mb of Inserted Sequence')+ylim(c(0,40))+
theme_bw()+
theme(legend.position = "none")

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig5A.png'
	,plot=Fig5A ,width=1.5,height=1.5,bg='#ffffff')

B73_TE_cumulMb <- all_B73_TE_n_calls %>% 
dplyr::mutate(size=end-start) %>% 
group_by(classification) %>% 
dplyr::summarise(TE_cumulMB=sum(size)/1000000)

B73_totalTE_cumulMb <- sum(B73_TE_cumulMb$TE_cumulMB)
B73_polymorphicTE_cumulMb <- B73_TE_cumulMb$TE_cumulMB[2]

NAM_TE_cumulMb <- all_NAM_TE_n_calls %>% 
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

Fig5B <- ggplot(TE_cumulMb,aes(x=group,y=perc_polymorphic)) +
geom_bar(aes(fill = group), position = "dodge", stat="identity",width=0.4)+
scale_fill_manual(labels = c("Genome-Wide","SNP Depleted"),values=c("#2b4655","#738e8e")) +
labs(x='',y='% Cumulative Mb of Polymorphic TE Sequence')+ylim(c(0,40))+
theme_bw()+
theme(legend.position = "none")

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig5B.png'
	,plot=Fig5B ,width=1.5,height=1.5,bg='#ffffff')

fig5AB <- ggarrange(Fig5A,Fig5B,nrow=1,ncol=2)

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig5AB.png'
	,plot=fig5AB ,width=4,height=4,bg='#ffffff')

combined_summary_AW_CA <- summary_AW_CA %>% group_by(Category,Type) %>%
dplyr::summarise(Amount=sum(Amount)) %>%
dplyr::mutate(tag='Genome Wide')

insertions_B73_SNPD_AW_overlaps <- B73_SNPD_AW_overlaps %>% subset(AW_type=='structural_insertion_inB73')

reduced_B73_AWBlocks_CA <- B73_AWBlocks_CA %>% select(AW_BlockID,Lineage_Comp,Category)
colnames(reduced_B73_AWBlocks_CA) <- c('AW_ID','AW_Comp','Category')

B73_SNPD_insertion_categories <- left_join(insertions_B73_SNPD_AW_overlaps, reduced_B73_AWBlocks_CA,by=c('AW_ID','AW_Comp'))
#fwrite(B73_SNPD_insertion_categories,'/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SNP_Depleted_Summaries/B73_SNP_Depleted_Insertion_Categories.csv')

SNPD_summary_B73_AWBlocks_CA <- B73_SNPD_insertion_categories %>% group_by(Category) %>%
dplyr::summarise(n=n(),Mb=sum(AW_size)/1000000)

summarise_SNPD_B73_insertion_CA <- summarise_AW_CA(SNPD_summary_B73_AWBlocks_CA,'B73') 

insertions_NAM_SNPD_AW_overlaps <- NAM_SNPD_AW_overlaps %>% subset(AW_type=='structural_insertion_inNAM')

reduced_NAM_AWBlocks_CA <- NAM_AWBlocks_CA %>% select(AW_BlockID,Lineage_Comp,Category)
colnames(reduced_NAM_AWBlocks_CA) <- c('AW_ID','AW_Comp','Category')

NAM_SNPD_insertion_categories <- left_join(insertions_NAM_SNPD_AW_overlaps, reduced_NAM_AWBlocks_CA,by=c('AW_ID','AW_Comp'))
#fwrite(NAM_SNPD_insertion_categories,'/home/brandvai/mmunasin/TE_Intron/store_data/summary_data_files/SNP_Depleted_Summaries/NAM_SNP_Depleted_Insertion_Categories.csv')

SNPD_summary_NAM_AWBlocks_CA <- NAM_SNPD_insertion_categories %>% group_by(Category) %>%
dplyr::summarise(n=n(),Mb=sum(AW_size)/1000000)
summarise_SNPD_NAM_insertion_CA <- summarise_AW_CA(SNPD_summary_NAM_AWBlocks_CA,'NAM') 

SNPD_summary_AW_CA <- rbind(summarise_SNPD_B73_insertion_CA,summarise_SNPD_NAM_insertion_CA)

combined_SNPD_summary_AW_CA <- SNPD_summary_AW_CA %>% group_by(Category,Type) %>%
dplyr::summarise(Amount=sum(Amount)) %>% dplyr::mutate(tag='SNP Depleted')

combined_AW_CA <- rbind(combined_summary_AW_CA,combined_SNPD_summary_AW_CA)

fig_5C <- ggplot(combined_AW_CA,aes(x=tag,y=Amount,fill=Category)) +
geom_bar(position='fill',stat='identity',width=0.3)+
facet_wrap(~Type,nrow=1,ncol=2)+
labs(x='Structural Insertions',y='Proportion',fill='Category')+
scale_fill_manual(values=met.brewer("Hokusai3", 5))+
theme_bw()+
theme(legend.position = "none")

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig5C.png'
	,plot=fig_5C ,width=2.5,height=2,bg='#ffffff')

fig_5D_updated <- combined_AW_CA %>% 
subset(Category!='Category 1') %>% 
dplyr::mutate(SV_type=case_when(
	Category=='Category 3' ~ 'Insertion',
	TRUE ~ 'Deletion')) %>%
dplyr::mutate(SV_type=factor(SV_type,levels=c("Insertion","Deletion"))) %>%
subset(Type=='Cumulative Mb') %>%
subset(tag=='SNP Depleted') %>%
ggplot(aes(x=SV_type,y=Amount,fill=Category)) +
geom_bar(stat='identity',position='stack',width=0.6) +
labs(x='SV Type',y='Cumulative Mb',fill='')+ 
scale_fill_manual(values=met.brewer("Hokusai3", 5)[2:5])+
theme_bw()+
theme(legend.position='none')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/fig5D_updated.png'
	,plot=fig_5D_updated ,width=2,height=2,bg='#ffffff')

family_superfamily_relation <- rbind(all_B73_TE_n_calls,all_NAM_TE_n_calls) %>% 
group_by(raw_family, condense_superfamily) %>% dplyr::summarise(n=n()) %>%
ungroup() %>% group_by(raw_family) %>% arrange(desc(n)) %>% 
filter(row_number()==1) %>% select(1,2)
colnames(family_superfamily_relation)[1] <- 'TE_family'

superfamily_palette <- met.brewer('Tam',n=9)
names(superfamily_palette) <- c('DTA','DTC','DTH','DTM','DTT','RIL','RLC','RLG','RLX')

B73_polymorphicTE_InsertionCategory_SNPD <- obtain_polymorphicTE_inSNPD(B73_SNPD_TE_overlaps,B73_SNPD_insertion_categories)
NAM_polymorphicTE_InsertionCategory_SNPD <- obtain_polymorphicTE_inSNPD(NAM_SNPD_TE_overlaps,NAM_SNPD_insertion_categories)

B73_Cat3_polymorphicTE_SNPD <- B73_polymorphicTE_InsertionCategory_SNPD %>%
subset(AW_Category=='Category_3') %>% dplyr::mutate(tag='B73_TE')
NAM_Cat3_polymorphicTE_SNPD <- NAM_polymorphicTE_InsertionCategory_SNPD %>%
subset(AW_Category=='Category_3')%>% dplyr::mutate(tag='NAM_TE')
Cat3_polymorphicTE_SNPD <- rbind(B73_Cat3_polymorphicTE_SNPD,NAM_Cat3_polymorphicTE_SNPD)

FigS5_R1_C1 <- Cat3_polymorphicTE_SNPD %>% group_by(tag,TE_superfamily) %>% 
dplyr::summarise(n=n()) %>%
ggplot(aes(x=tag,y=n,fill=TE_superfamily))+
geom_bar(position='stack',stat='identity')+
labs(x='Lineage',y='Count')+
scale_fill_manual(values=superfamily_palette)+
theme_bw()+
theme(legend.position='none')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS5_R1C1.pdf'
	,plot=FigS5_R1_C1,width=2,height=3,bg='#ffffff')

FigS5_R1_C2 <- Cat3_polymorphicTE_SNPD %>% subset(tag=='NAM_TE') %>%
group_by(AW_LC,TE_superfamily) %>% 
dplyr::summarise(n=n()) %>%
ggplot(aes(x=AW_LC,y=n,fill=TE_superfamily))+
geom_bar(position='stack',stat='identity')+
labs(x='NAM Line',y='Count')+
scale_fill_manual(values=superfamily_palette)+
theme_bw()+
theme(legend.position='none',axis.text.x = element_text(angle = 90,size=5, vjust = 0.5, hjust=1))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS5_R1C2.pdf'
	,plot=FigS5_R1_C2,width=4,height=3,bg='#ffffff')

FigS5_R2_C1 <- Cat3_polymorphicTE_SNPD %>%
dplyr::mutate(TE_chr=factor(TE_chr,levels=chr_level_orders)) %>%
group_by(tag,TE_chr) %>% dplyr::summarise(n=n()) %>%
ggplot(aes(x=tag,y=n,fill=TE_chr))+
geom_bar(position='stack',stat='identity')+
labs(x='Lineage',y='Count')+
scale_fill_manual(values=met.brewer('Pillement',n=10))+
theme_bw()+
theme(legend.position='none')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS5_R2C1.pdf'
	,plot=FigS5_R2_C1,width=2,height=3,bg='#ffffff')

FigS5_R2_C2 <- Cat3_polymorphicTE_SNPD %>% subset(tag=='NAM_TE') %>%
dplyr::mutate(TE_chr=factor(TE_chr,levels=chr_level_orders)) %>%
group_by(AW_LC,TE_chr) %>% 
dplyr::summarise(n=n()) %>%
ggplot(aes(x=AW_LC,y=n,fill=TE_chr))+
geom_bar(position='stack',stat='identity')+
labs(x='NAM Line',y='Count')+
scale_fill_manual(values=met.brewer('Pillement',n=10))+
theme_bw()+
theme(legend.position='none',axis.text.x = element_text(angle = 90,size=5, vjust = 0.5, hjust=1))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS5_R2C2.pdf'
	,plot=FigS5_R2_C2,width=4,height=3,bg='#ffffff')

top4_SNPD_families <- c("CRM2_7577nt","DTA_ZM00383_consensus","ji_AC204382_8228","ji_AC215728_13156")


figS6_R1C1 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=='CRM2_7577nt') %>% 
subset(TE_size>7500) %>%
ggplot(aes(x=as.numeric(TE_size))) +
geom_histogram(binwidth=5,fill='#5C1E5C') +
labs(x='TE Size',y='Count') +
theme_bw()
ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R1C1.pdf'
	,plot=figS6_R1C1,width=2,height=3,bg='#ffffff')

figS6_R1C2 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=='CRM2_7577nt') %>% 
subset(TE_size>7500) %>% 
group_by(tag,TE_chr) %>% dplyr::summarise(n=n()) %>%
dplyr::mutate(TE_chr=factor(TE_chr,levels=chr_level_orders)) %>%
ggplot(aes(x=tag,y=n,fill=TE_chr)) +
geom_bar(position='stack',stat='identity')+
labs(x='Lineage',y='Count')+
scale_fill_manual(values=met.brewer('Pillement',n=10))+
theme_bw()+
theme(legend.position='none')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R1C2.pdf'
	,plot=figS6_R1C2,width=2,height=3,bg='#ffffff')

figS6_R1C3 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=='CRM2_7577nt') %>% 
subset(TE_size>7500) %>% subset(tag=='NAM_TE') %>%
group_by(AW_LC,TE_chr) %>% dplyr::summarise(n=n()) %>%
dplyr::mutate(TE_chr=factor(TE_chr,levels=chr_level_orders)) %>%
ggplot(aes(x=AW_LC,y=n,fill=TE_chr)) +
geom_bar(position='stack',stat='identity')+
labs(x='NAM Line',y='Count')+
scale_fill_manual(values=met.brewer('Pillement',n=10))+
theme_bw()+
theme(legend.position='none',axis.text.x = element_text(angle = 90,size=5, vjust = 0.5, hjust=1))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R1C3.pdf'
	,plot=figS6_R1C3,width=3,height=3,bg='#ffffff')

figS6_R2C1 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=="DTA_ZM00383_consensus") %>% 
ggplot(aes(x=as.numeric(TE_size))) +
geom_histogram(binwidth=1,fill='#FFD353') +
labs(x='TE Size',y='Count') +
scale_x_continuous(breaks=seq(395,407,by=2))+
theme_bw()
ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R2C1.pdf'
	,plot=figS6_R2C1,width=2,height=3,bg='#ffffff')

figS6_R2C2 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=="DTA_ZM00383_consensus") %>% 
group_by(tag,TE_chr) %>% dplyr::summarise(n=n()) %>%
dplyr::mutate(TE_chr=factor(TE_chr,levels=chr_level_orders)) %>%
ggplot(aes(x=tag,y=n,fill=TE_chr)) +
geom_bar(position='stack',stat='identity')+
labs(x='Lineage',y='Count')+
scale_fill_manual(values=met.brewer('Pillement',n=10))+
theme_bw()+
theme(legend.position='none')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R2C2.pdf'
	,plot=figS6_R2C2,width=2,height=3,bg='#ffffff')

figS6_R2C3 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=="DTA_ZM00383_consensus") %>% 
subset(tag=='NAM_TE') %>%
group_by(AW_LC,TE_chr) %>% dplyr::summarise(n=n()) %>%
dplyr::mutate(TE_chr=factor(TE_chr,levels=chr_level_orders)) %>%
ggplot(aes(x=AW_LC,y=n,fill=TE_chr)) +
geom_bar(position='stack',stat='identity')+
labs(x='NAM Line',y='Count')+
scale_fill_manual(values=met.brewer('Pillement',n=10))+
theme_bw()+
theme(legend.position='none',axis.text.x = element_text(angle = 90,size=5, vjust = 0.5, hjust=1))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R2C3.pdf'
	,plot=figS6_R2C3,width=3,height=3,bg='#ffffff')

figS6_R3C1 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=="ji_AC204382_8228") %>% 
subset(TE_size > 8400) %>%
ggplot(aes(x=as.numeric(TE_size))) +
geom_histogram(binwidth=50,fill='#341648') +
labs(x='TE Size',y='Count') +
theme_bw()+
theme(legend.position='none',axis.text.x = element_text(angle = 90,size=5, vjust = 0.5, hjust=1))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R3C1.pdf'
	,plot=figS6_R3C1,width=2,height=3,bg='#ffffff')

figS6_R3C2 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=="ji_AC204382_8228") %>% 
subset(TE_size > 8400) %>% 
group_by(tag,TE_chr) %>% dplyr::summarise(n=n()) %>%
dplyr::mutate(TE_chr=factor(TE_chr,levels=chr_level_orders)) %>%
ggplot(aes(x=tag,y=n,fill=TE_chr)) +
geom_bar(position='stack',stat='identity')+
labs(x='Lineage',y='Count')+
scale_fill_manual(values=met.brewer('Pillement',n=10))+
theme_bw()+
theme(legend.position='none')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R3C2.pdf'
	,plot=figS6_R3C2,width=2,height=3,bg='#ffffff')

figS6_R3C3 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=="ji_AC204382_8228") %>% 
subset(TE_size > 8400) %>% 
subset(tag=='NAM_TE') %>%
group_by(AW_LC,TE_chr) %>% dplyr::summarise(n=n()) %>%
dplyr::mutate(TE_chr=factor(TE_chr,levels=chr_level_orders)) %>%
ggplot(aes(x=AW_LC,y=n,fill=TE_chr)) +
geom_bar(position='stack',stat='identity')+
labs(x='NAM Line',y='Count')+
scale_fill_manual(values=met.brewer('Pillement',n=10))+
theme_bw()+
theme(legend.position='none',axis.text.x = element_text(angle = 90,size=5, vjust = 0.5, hjust=1))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R3C3.pdf'
	,plot=figS6_R3C3,width=3,height=3,bg='#ffffff')

figS6_R4C1 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=="ji_AC215728_13156") %>% 
ggplot(aes(x=as.numeric(TE_size))) +
geom_histogram(binwidth=5,fill='#8F2957') +
labs(x='TE Size',y='Count') +
theme_bw()+
theme(legend.position='none',axis.text.x = element_text(angle = 90,size=5, vjust = 0.5, hjust=1))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R4C1.pdf'
	,plot=figS6_R4C1,width=2,height=3,bg='#ffffff')

figS6_R4C2 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=="ji_AC215728_13156") %>% 
group_by(tag,TE_chr) %>% dplyr::summarise(n=n()) %>%
dplyr::mutate(TE_chr=factor(TE_chr,levels=chr_level_orders)) %>%
ggplot(aes(x=tag,y=n,fill=TE_chr)) +
geom_bar(position='stack',stat='identity')+
labs(x='Lineage',y='Count')+
scale_fill_manual(values=met.brewer('Pillement',n=10))+
theme_bw()+
theme(legend.position='none')

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R4C2.pdf'
	,plot=figS6_R4C2,width=2,height=3,bg='#ffffff')

figS6_R4C3 <- Cat3_polymorphicTE_SNPD %>% 
subset(TE_family=="ji_AC215728_13156") %>% 
subset(tag=='NAM_TE') %>%
group_by(AW_LC,TE_chr) %>% dplyr::summarise(n=n()) %>%
dplyr::mutate(TE_chr=factor(TE_chr,levels=chr_level_orders)) %>%
ggplot(aes(x=AW_LC,y=n,fill=TE_chr)) +
geom_bar(position='stack',stat='identity')+
labs(x='NAM Line',y='Count')+
scale_fill_manual(values=met.brewer('Pillement',n=10))+
theme_bw()+
theme(legend.position='none',axis.text.x = element_text(angle = 90,size=5, vjust = 0.5, hjust=1))

ggsave('/home/brandvai/mmunasin/TE_Intron/store_data/Rplots/08.31.22_figures/figS6_R4C3.pdf'
	,plot=figS6_R4C3,width=3,height=3,bg='#ffffff')
