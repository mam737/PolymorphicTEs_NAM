.libPaths('/path/to/mmunasin/Rlibs4')

library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(Rmisc)
library(forcats)
library(ggforce)
library(MetBrewer)

source("~/path/to/figures/generate_figure_functions.R")

################################################################################################
################################################################################################
############### color palette notes ############################################################
################################################################################################
################################################################################################

# We're going to use colors from the R NatParksPalette palette: DeathValley
# The following three colors refer to AW Block Assignment: alignable(#18315a), structural insertion(#cb6a2d), and unalignable (#68434e)
# The following three colors refer to TE Class Assignment: shared(#416191), polymorphic (#fbb25d), and ambiguous (#b57f84)
# The following four colors refer to shared frequency classifications: Private ("#d1b252"), Variable ("#a97f2f"), Near-Core("#7e5522"), Core("#472c0b")
# This ordering of the NAM lines follows those in the NAM paper

NAM_level_orders <- c("B73","B97","Ky21","M162W","Ms71","Oh43","Oh7B","M37W","Mo18W","Tx303","HP301","P39","Il14H","CML52","CML69","CML103","CML228","CML247","CML277","CML322","CML333","Ki3","Ki11","NC350","NC358","Tzi8")
chr_level_orders <- c("chr1","chr2","chr3",'chr4','chr5','chr6','chr7','chr8','chr9','chr10')
classification_level_orders <- c("shared",'polymorphic','ambiguous')

################################################################################################
################################################################################################
##################### Table S1 #################################################################
################################################################################################
################################################################################################

filtered_TE_anno_dir <- '/path/to/filtered_TE_anno'
filtered_TE_anno <- load_filtered_TE_data(filtered_TE_anno_dir)

TE_summary <- filtered_TE_summary(filtered_TE_anno)

fwrite(TE_summary,file='/path/to/summary_data_files/SupplementalTables/TableS1.csv')

B73_TEs <- filtered_TE_anno %>% filter(lineage=='B73') %>%
select(chr,start,end,name,class,condense_superfamily,raw_family,method)

################################################################################################
################################################################################################
############### Figure 1B + S1 #################################################################
################################################################################################
################################################################################################

head_AW_dir <- '/path/to/summarised_AnchorWave_Regions'
all_AW_data <- load_AW_data(head_AW_dir)

summarised_all_AW_data <- all_AW_data %>% 
group_by(Lineage_Comp,reduced_type) %>% 
dplyr::summarise(ID_BlockTotal_Mb=sum(ID_BlockSize)/1000000, ASM_BlockTotal_Mb=sum(ASM_BlockSize)/1000000) %>% 
data.frame()

B73_AW_total_means <- summarised_all_AW_data %>% 
group_by(Lineage_Comp) %>% 
summarySE(measurevar='ID_BlockTotal_Mb',groupvars=c('reduced_type'))%>% 
subset(reduced_type != 'structural_insertion_inNAM')

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

ggsave('/path/to/figures/fig1b.pdf'
	,plot=fig_1b_1 ,width=1.5,height=2.0)

average_all_AW_data <- summarised_all_AW_data %>% 
group_by(reduced_type) %>% 
dplyr::summarise(ID_BlockTotal_Mb=mean(ID_BlockTotal_Mb),ASM_BlockTotal_Mb=mean(ASM_BlockTotal_Mb)) %>% 
dplyr::mutate(Lineage_Comp='Average') %>% select(4,1,2,3)

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

ggsave('/path/to/figures/figS1.pdf'
	,plot=figs1,width=8,height=6)

################################################################################################
################################################################################################
############### Figure 1C #################################################################
################################################################################################
################################################################################################

B73_TE_n_dir <- '/path/to/store_data/polymorphic_TE_calls/B73/'
NAM_TE_n_dir <- '/path/to/store_data/polymorphic_TE_calls/NAM/'

all_B73_TE_n_calls <- process_TE_calls(B73_TE_n_dir)
all_NAM_TE_n_calls <- process_TE_calls(NAM_TE_n_dir)

summarised_TE_B73 <- summarise_TE_Tag(all_B73_TE_n_calls,'B73')
summarised_TE_NAM <- summarise_TE_Tag(all_NAM_TE_n_calls,'NAM')

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

B73_exon_dir <- '/path/to/polymorphic_gene_calls/B73/exon_calls'
B73_full_dir <- '/path/to/polymorphic_gene_calls/B73/full_calls'
NAM_exon_dir <- '/path/to/polymorphic_gene_calls/NAM/exon_calls'
NAM_full_dir <- '/path/to/polymorphic_gene_calls/NAM/full_calls'

all_B73_exon_calls <- process_gene_calls(B73_exon_dir,tag='id',type='exon')
all_B73_full_calls <- process_gene_calls(B73_full_dir,tag='id',type='gene')
all_NAM_exon_calls <- process_gene_calls(NAM_exon_dir,tag='asm',type='exon')
all_NAM_full_calls <- process_gene_calls(NAM_full_dir,tag='asm',type='gene')

syntenic_gene_list_dir <- '/path/to/NAM_syntenic_genes'
syntenic_genes_all <- obtain_gene_synteny(syntenic_gene_list_dir)

all_B73_exon_calls <- add_synteny_column(all_B73_exon_calls,syntenic_genes_all)
all_B73_full_calls <- add_synteny_column(all_B73_full_calls,syntenic_genes_all)
all_NAM_exon_calls <- add_synteny_column(all_NAM_exon_calls,syntenic_genes_all)
all_NAM_full_calls <- add_synteny_column(all_NAM_full_calls,syntenic_genes_all)

summarised_nsg_B73_full <- obtain_synteny_summary(all_B73_full_calls,synteny_tag=FALSE,type_tag='id')
summarised_nsg_NAM_full <- obtain_synteny_summary(all_NAM_full_calls,synteny_tag=FALSE,type_tag='asm')
summarised_sg_B73_full <- obtain_synteny_summary(all_B73_full_calls,synteny_tag=TRUE,type_tag='id')
summarised_sg_NAM_full <- obtain_synteny_summary(all_NAM_full_calls,synteny_tag=TRUE,type_tag='asm')

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

ggsave('/path/to/figures/fig1cde.pdf'
	,plot=fig1cde,width=4,height=2.0,bg='#ffffff')

################################################################################################
################################################################################################
#################### Figure S2 #################################################################
################################################################################################
################################################################################################

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

ggsave('/path/to/figures/figS2.pdf'
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

fig2A <- B73_TE_call_dist %>% group_by(n,degree_shared) %>% dplyr::summarise(count_n=n()) %>%
ggplot(aes(x=n,y=count_n,fill=degree_shared)) +
geom_bar(stat='identity') + 
scale_fill_manual(labels=c('Private','Variable','Near Core','Core'),values=c("#d1b252","#a97f2f","#7e5522","#472c0b"))+
labs(x='Times Classified as Shared Between B73 and NAM',y='Number of TEs') + 
theme_bw()+
theme(legend.position = "bottom",legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size=7),axis.text.x = element_text(angle = 90,hjust=1.10,vjust=0.95))+
guides(fill=guide_legend(title=''))

ggsave('/path/to/figures/fig2A.tif'
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

ggsave('/path/to/figures/fig2C.tif'
	,plot=fig2C,width=4,height=2,device='tiff',dpi=700,bg='#ffffff')

structuralTE_identity <- fread('/path/to/summary_data_files/structurally_annotated_LTRs_identity.tsv')
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

ggsave('/path/to/figures/fig2D_alt2.png'
	,plot=fig2D_alt2,width=4,height=2,bg='#ffffff')
ggsave('/path/to/figures/fig2D_alt2.tif'
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

ggsave('/path/to/figures/fig3A.png'
	,plot=fig3A,width=6.0,height=1.5,bg='#ffffff')

B73_largefam_polymorphic_percent <- obtain_large_polymorphic_percent(all_B73_TE_n_calls,tag.string='B73')
NAM_largefam_polymorphic_percent <- obtain_large_polymorphic_percent(all_NAM_TE_n_calls,tag.string='NAM')

largefam_polymorphic_percent <- rbind(B73_largefam_polymorphic_percent,NAM_largefam_polymorphic_percent)

fig3B_alt1 <- B73_largefam_polymorphic_percent %>% group_by(raw_family,fam_size) %>%
dplyr::summarise(range_perc_total=max(percent_total)-min(percent_total)) %>% 
ggplot(aes(x=range_perc_total)) + geom_histogram(binwidth=5,fill='#e8b697' ,color='#CB6A2D')+
labs(x='% Polymorphic Range',y='Number of B73 TE Families')+
theme_bw()

ggsave('/path/to/figures/fig3B_alt1.png'
	,plot=fig3B_alt1,width=4.0,height=3.5,bg='#ffffff')

fig3B_alt2 <- B73_largefam_polymorphic_percent %>% group_by(raw_family,fam_size) %>%
dplyr::summarise(range_perc_total=max(percent_total)-min(percent_total)) %>% 
ggplot(aes(x=fam_size,y=range_perc_total)) + geom_point(size=0.5,color='#CB6A2D')+
labs(x='TE Family Size',y='Range of Polymorphic %')+
theme_bw()

ggsave('/path/to/figures/fig3B_alt2.png'
	,plot=fig3B_alt2,width=4.0,height=3.5,bg='#ffffff')

largefam_range <- B73_largefam_polymorphic_percent %>% group_by(raw_family,fam_size) %>%
dplyr::summarise(range_perc_total=max(percent_total)-min(percent_total))

B73_structuralLTR_age <- B73_structuralTE_identity %>% 
select(name,ltr_identity,condense_superfamily,raw_family) %>%
group_by(raw_family)

fig3X_ltrage <- left_join(B73_structuralLTR_age,largefam_range,by=c("raw_family")) %>%
filter(!(is.na(fam_size))) %>% group_by(raw_family,fam_size,range_perc_total) %>%
dplyr::summarise(avg_struct_age=mean(ltr_identity))%>% 
dplyr::mutate(avg_ltr_age=case_when(
	avg_struct_age==1 ~ 'Very Young',
	avg_struct_age >= 0.95 & avg_struct_age < 1 ~ 'Young',
	avg_struct_age >= 0.9 & avg_struct_age < 0.95 ~ 'Moderate',
	avg_struct_age >= 0.8 & avg_struct_age < 0.9 ~ 'Old')) %>%
dplyr::mutate(avg_ltr_age=factor(avg_ltr_age,levels=c('Very Young','Young','Moderate','Old'))) %>%
ggplot(aes(x=fam_size,y=range_perc_total,color=avg_ltr_age)) + 
geom_point(size=0.5)+
labs(x='TE Family Size',y='Range of Polymorphic %')+
scale_color_manual(values=rev(met.brewer("Hokusai2",n=4)))+
theme_bw()

ggsave('/path/to/figures/fig3X_ltrage.png'
	,plot=fig3X_ltrage,width=4.0,height=3.5,bg='#ffffff')

fig3X2_ltrage <- left_join(B73_structuralLTR_age,largefam_range,by=c("raw_family")) %>%
filter(!(is.na(fam_size))) %>% group_by(raw_family,fam_size,range_perc_total) %>%
dplyr::summarise(avg_struct_age=mean(ltr_identity))%>% 
dplyr::mutate(avg_ltr_age=case_when(
	avg_struct_age==1 ~ 'Very Young',
	avg_struct_age >= 0.95 & avg_struct_age < 1 ~ 'Young',
	avg_struct_age >= 0.9 & avg_struct_age < 0.95 ~ 'Moderate',
	avg_struct_age >= 0.8 & avg_struct_age < 0.9 ~ 'Old')) %>%
dplyr::mutate(avg_ltr_age=factor(avg_ltr_age,levels=c('Very Young','Young','Moderate','Old'))) %>%
ggplot(aes(x=fam_size,y=range_perc_total,color=avg_ltr_age)) + 
geom_point(size=0.5)+facet_wrap(~avg_ltr_age)+
labs(x='TE Family Size',y='Range of Polymorphic %')+
scale_color_manual(values=rev(met.brewer("Hokusai2",n=4)))+
theme_bw()

ggsave('/path/to/figures/fig3X2_ltrage.png'
	,plot=fig3X2_ltrage,width=5.0,height=5,bg='#ffffff')

left_join(B73_structuralLTR_age,largefam_range,by=c("raw_family")) %>%
filter(!(is.na(fam_size))) %>% group_by(raw_family,fam_size,range_perc_total) %>%
dplyr::summarise(avg_struct_age=mean(ltr_identity))%>% 
dplyr::mutate(avg_ltr_age=case_when(
	avg_struct_age==1 ~ 'Very Young',
	avg_struct_age >= 0.95 & avg_struct_age < 1 ~ 'Young',
	avg_struct_age >= 0.9 & avg_struct_age < 0.95 ~ 'Moderate',
	avg_struct_age >= 0.8 & avg_struct_age < 0.9 ~ 'Old')) %>%
dplyr::mutate(avg_ltr_age=factor(avg_ltr_age,levels=c('Very Young','Young','Moderate','Old'))) %>%
group_by(avg_ltr_age) %>% dplyr::summarise(min=min(range_perc_total),max=max(range_perc_total),n=n())

################################################################################################
################################################################################################
############### Figure 4 #################################################################
################################################################################################
################################################################################################

AWBlocks_CategoryAssignment <- fread("/path/to/summary_data_files/all_AWBlocks_CategoryAssignment.tsv")

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

ggsave('/path/to/figures/fig4C.png'
	,plot=fig_4C ,width=2,height=4,bg='#ffffff')

################################################################################################
################################################################################################
############### Figure 5 #################################################################
################################################################################################
################################################################################################

# Read in data for AW Blocks in SNP Depleted Regions
B73_SNPD_AW_overlaps <- fread('/path/to/summary_data_files/SNP_Depleted_Summaries/B73_SNP_Depleted_AW_Overlaps.csv')
NAM_SNPD_AW_overlaps <- fread('/path/to/summary_data_files/SNP_Depleted_Summaries/NAM_SNP_Depleted_AW_Overlaps.csv')

cumulMb <- obtain_cumulMB(summarised_all_AW_data,B73_SNPD_AW_overlaps,NAM_SNPD_AW_overlaps)

Fig5A <- ggplot(cumulMb,aes(x=group,y=perc_inserted)) +
geom_bar(aes(fill = group), position = "dodge", stat="identity",width=0.4)+
scale_fill_manual(labels = c("Genome-Wide","SNP Depleted"),values=c("#2b4655","#738e8e")) +
labs(x='',y='% Cumulative Mb of Inserted Sequence')+ylim(c(0,40))+
theme_bw()+
theme(legend.position = "none")

ggsave('/path/to/figures/fig5A.png'
	,plot=Fig5A ,width=1.5,height=1.5,bg='#ffffff')

TE_cumulMb <- obtain_TE_cumulMb(all_B73_TE_n_calls,all_NAM_TE_n_calls)

Fig5B <- ggplot(TE_cumulMb,aes(x=group,y=perc_polymorphic)) +
geom_bar(aes(fill = group), position = "dodge", stat="identity",width=0.4)+
scale_fill_manual(labels = c("Genome-Wide","SNP Depleted"),values=c("#2b4655","#738e8e")) +
labs(x='',y='% Cumulative Mb of Polymorphic TE Sequence')+ylim(c(0,40))+
theme_bw()+
theme(legend.position = "none")

ggsave('/path/to/figures/fig5B.png'
	,plot=Fig5B ,width=1.5,height=1.5,bg='#ffffff')

B73_SNPD_SV_categories <- SNPD_SV_categories(B73_SNPD_AW_overlaps,B73_AWBlocks_CA,tag='B73')
NAM_SNPD_SV_categories <- SNPD_SV_categories(NAM_SNPD_AW_overlaps,NAM_AWBlocks_CA,tag='NAM')

combined_AW_CA <- AW_by_CA(summary_AW_CA,B73_SNPD_SV_categories,NAM_SNPD_SV_categories)

fig_5C <- ggplot(combined_AW_CA,aes(x=tag,y=Amount,fill=Category)) +
geom_bar(position='fill',stat='identity',width=0.3)+
facet_wrap(~Type,nrow=1,ncol=2)+
labs(x='Structural Insertions',y='Proportion',fill='Category')+
scale_fill_manual(values=met.brewer("Hokusai3", 5))+
theme_bw()+
theme(legend.position = "none")

ggsave('/path/to/figures/fig5C.png'
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

ggsave('/path/to/figures/fig5D_updated.png'
	,plot=fig_5D_updated ,width=2,height=2,bg='#ffffff')

B73_SNPD_size_dist <- B73_SNPD_SV_categories %>% 
filter(Category !='Category_1') %>% 
dplyr::mutate(SV_type=case_when(
	Category=='Category_3' ~ 'Insertion',
	TRUE ~ 'Deletion')) %>%
select(AW_chr,AW_start,AW_end,AW_ID,AW_Comp,AW_size,Category,SV_type)

NAM_SNPD_size_dist <- NAM_SNPD_SV_categories %>% 
filter(Category !='Category_1') %>% 
dplyr::mutate(SV_type=case_when(
	Category=='Category_3' ~ 'Insertion',
	TRUE ~ 'Deletion')) %>%
select(AW_chr,AW_start,AW_end,AW_ID,AW_Comp,AW_size,Category,SV_type)

size_dist <- rbind(B73_SNPD_size_dist,NAM_SNPD_size_dist)

################################################################################################
################################################################################################
############### Figure S3 + S4 #################################################################
################################################################################################
################################################################################################

SNP_dir <- '/path/to/SNP_density_summaries'
SNP_data <- lapply(list.files(SNP_dir,full.names=T),fread)

file_order <- unlist(lapply(list.files(SNP_dir,full.names=T),function(x) unlist(str_split(x,pattern='_'))[6]))
lineage_names <- substring(file_order,1,nchar(file_order)-3)
SNP_data <- mapply(cbind,SNP_data,'ASM_Comp'=lineage_names,SIMPLIFY=F)
SNP_data <- do.call('rbind',SNP_data)
B73_SNP_data <- SNP_data %>% subset(ASM_Comp=='B97')

SNPD_regions <- fread('/path/to/summary_data_files/SupplementalTables/TableS2.csv')
B73_SNPD_regions <- SNPD_regions %>% subset(NAM_Comparison=='B97')
B73_SNPD_highlights <- B73_SNPD_regions %>% select(B73_Start_Pos,B73_End_Pos,Chr) %>% dplyr::summarise(start=B73_Start_Pos/1000000,end=B73_End_Pos/1000000,chr=Chr) %>% dplyr::mutate(group=1:n())

FigS3 <- ggplot(B73_SNP_data,aes(x=BinStart/1000000,y=norm_SNP_count)) +
geom_line()+ ylim(c(0,0.025)) + geom_rect(data=B73_SNPD_highlights,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=0,ymax=0.025,group=group),color='transparent',fill='orange',alpha=0.3)+
facet_wrap(~chr,nrow=5,ncol=2,scales='free_x')+
labs(x='Sliding Window Mb Start Position',y='Normalized SNP Count',title='Normalized SNP Rate for B73 vs B97')

ggsave('/path/to/figures/figS3_final.pdf'
	,plot=FigS3 ,width=7,height=10,bg='#ffffff')

SNPD_size_dist <- data.frame(SNPD_size=c(SNPD_regions$B73_Size_Mb,SNPD_regions$NAM_Size_Mb))
FigS4 <- ggplot(SNPD_size_dist,aes(x=SNPD_size)) +
geom_histogram(binwidth=0.5,fill='lightblue',color='darkblue')+labs(x='SNP Depleted Region Size (Mb)',y='Count')

ggsave('/path/to/figures/figS4_final.pdf'
	,plot=FigS4 ,width=3,height=2,bg='#ffffff')
################################################################################################
################################################################################################
############### Figure 6 #################################################################
################################################################################################
################################################################################################

family_superfamily_relation <- fam_supfam_relate(all_B73_TE_n_calls,all_NAM_TE_n_calls)

superfamily_palette <- met.brewer('Tam',n=9)
names(superfamily_palette) <- c('DTA','DTC','DTH','DTM','DTT','RIL','RLC','RLG','RLX')

B73_polymorphicTE_InsertionCategory_SNPD <- TE_SV_relationship(B73_SNPD_TE_overlaps,B73_SNPD_insertion_categories)
NAM_polymorphicTE_InsertionCategory_SNPD <- TE_SV_relationship(NAM_SNPD_TE_overlaps,NAM_SNPD_insertion_categories)

B73_Cat3_polymorphicTE_SNPD <- B73_polymorphicTE_InsertionCategory_SNPD %>%
subset(AW_Category=='Category_3') %>% dplyr::mutate(tag='B73_TE')
NAM_Cat3_polymorphicTE_SNPD <- NAM_polymorphicTE_InsertionCategory_SNPD %>%
subset(AW_Category=='Category_3')%>% dplyr::mutate(tag='NAM_TE')
Cat3_polymorphicTE_SNPD <- rbind(B73_Cat3_polymorphicTE_SNPD,NAM_Cat3_polymorphicTE_SNPD)

Cat3_polymorphicTE_SNPD <- rbind(B73_Cat3_polymorphicTE_SNPD,NAM_Cat3_polymorphicTE_SNPD)
TableS3 <- Cat3_polymorphicTE_SNPD %>%
select(tag,TE_LC,TE_chr,TE_start,TE_end,TE_name,TE_class,TE_superfamily,TE_family,TE_method,TE_size,AW_chr,AW_start,AW_end,AW_BlockID,AW_Size,AW_SNPD_ID)
colnames(TableS3)[1:2] <- c("TE_Geno","AW_Comparison_Geno")
colnames(TableS3)[17] <- c('SNPD_BlockID')
fwrite(TableS3,file='/path/to/summary_data_files/SupplementalTables/TableS3.csv')

Cat3_Family_SNPD <- Cat3_polymorphicTE_SNPD %>% 
group_by(TE_family) %>% dplyr::summarise(n=n()) %>% 
left_join(family_superfamily_relation,by='TE_family')

Large_Cat3_Family_SNPD <- Cat3_Family_SNPD %>% subset(n > 5)
Small_Cat3_Family_Summary_SNPD <- Cat3_Family_SNPD %>% subset(n <= 5) %>%
group_by(condense_superfamily) %>% dplyr::summarise(n=sum(n)) %>% 
dplyr::mutate(TE_family='Small Family Summary',.before=condense_superfamily)

Family_Order <- c("CRM2_7577nt","DTA_ZM00383_consensus","ji_AC204382_8228","ji_AC215728_13156","DTC_ZM00089_consensus","Small Family Summary")
Fig6A <- rbind(Large_Cat3_Family_SNPD,Small_Cat3_Family_Summary_SNPD) %>%
dplyr::mutate(TE_family=factor(TE_family,levels=Family_Order)) %>%
ggplot(aes(x=TE_family,y=n,fill=condense_superfamily))+
geom_bar(position='stack',stat='identity',width=0.5)+
scale_fill_manual(name='condense_superfamily',values=superfamily_palette)+
labs(x='TE Family',y='')+
theme_classic()+
theme(legend.position="none",axis.text.x=element_blank())

ggsave('/path/to/fig6A.png'
	,plot=Fig6A,width=6,height=1.5,bg='#ffffff')


#########
B73_MultiTE_SNPD <- B73_SNPD_insertion_categories %>% filter(Category=='Category_4') %>%
select(AW_chr,AW_start,AW_end,AW_ID,AW_Comp,SNPD_ID,SNPD_chr,SNPD_start,SNPD_end)

NAM_MultiTE_SNPD <- NAM_SNPD_insertion_categories %>% filter(Category=='Category_4')

B73_MultiTE_SNPD_overlaps <- bedtoolsr::bt.intersect(B73_MultiTE_SNPD,B73_TEs,wo=T)
colnames(B73_MultiTE_SNPD_overlaps) <- c('AW_chr','AW_start','AW_end','AW_ID','AW_Comp','SNPD_ID','SNPD_chr','SNPD_start','SNPD_end','TE_chr','TE_start','TE_end','TE_name','TE_class','TE_superfamily','TE_family','TE_Method','bp_overlap')

B73_MultiTE_SNPD_overlaps <- B73_MultiTE_SNPD_overlaps %>%
dplyr::mutate(left_overlap = case_when(
	TE_start < AW_start & TE_end > AW_start ~ TRUE,
	TRUE ~ FALSE)) %>%
dplyr::mutate(right_overlap=case_when(
	TE_start < AW_end & TE_end > AW_end ~ TRUE,
	TRUE ~ FALSE)) %>% group_by(AW_ID,AW_Comp) %>%
dplyr::mutate(edge_case=case_when(
	TRUE %in% left_overlap & TRUE %in% right_overlap ~ 'both_edges',
	TRUE %in% left_overlap & !(TRUE %in% right_overlap) ~ 'left_edges',
	!(TRUE %in% left_overlap) & TRUE %in% right_overlap ~ 'right_edges',
	TRUE ~ 'no_edge'
	))

B73_MultiTE_SNPD_edge_matches <- B73_MultiTE_SNPD_overlaps %>% filter(edge_case=='both_edges') %>%
filter(left_overlap==TRUE | right_overlap==TRUE) %>% group_by(AW_ID,AW_Comp) %>%
dplyr::mutate(n=n()) %>% filter(n==2) %>% 
dplyr::mutate(edge_fammatch=case_when(
	length(unique(TE_family))==1 ~ TRUE,
	TRUE ~ FALSE)) %>% 
dplyr::mutate(edge_supermatch=case_when(
	length(unique(TE_superfamily))==1 ~ TRUE,
	TRUE ~ FALSE))


left_bounds <- B73_MultiTE_SNPD %>% 
dplyr::mutate(l_lbound = AW_start - 50, l_rbound = AW_start + 50, r_lbound = AW_end - 50, r_rbound = AW_end + 50) %>%
select(AW_chr,l_lbound,l_rbound,AW_ID,AW_start,AW_end,AW_Comp, SNPD_ID)

B73_MultiTE_SNPD_lbound_overlaps <- bedtoolsr::bt.intersect(left_bounds,B73_TEs,wo=T)
colnames(B73_MultiTE_SNPD_lbound_overlaps) <- c('AW_chr','l_bound','r_bound','AW_ID','AW_start','AW_end','AW_Comp','SNPD_ID','TE_chr','TE_start','TE_end','TE_name','TE_class','TE_superfamily','TE_family','TE_Method', 'TE_identity','bp_overlap')

lbound_test <- B73_MultiTE_SNPD_lbound_overlaps %>% 
group_by(AW_ID,AW_Comp) %>% dplyr::mutate(n=n()) %>%
dplyr::mutate(fl_edge_overlap=case_when(
	n==1 ~ TRUE,
	length(unique(TE_family))==1 ~ TRUE,
	TRUE ~ FALSE)) #%>%
#ungroup() %>% group_by(AW_ID,AW_Comp,fl_edge_overlap) %>%
#select(AW_ID,AW_Comp,fl_edge_overlap)%>%
#distinct(AW_ID,AW_Comp,fl_edge_overlap,.keep_all=T)

right_bounds <- B73_MultiTE_SNPD %>% 
dplyr::mutate(l_lbound = AW_start - 50, l_rbound = AW_start + 50, r_lbound = AW_end - 50, r_rbound = AW_end + 50) %>%
select(AW_chr,r_lbound,r_rbound,AW_ID,AW_start,AW_end,AW_Comp, SNPD_ID)

B73_MultiTE_SNPD_rbound_overlaps <- bedtoolsr::bt.intersect(right_bounds,B73_TEs,wo=T)
colnames(B73_MultiTE_SNPD_rbound_overlaps) <- c('AW_chr','l_bound','r_bound','AW_ID','AW_start','AW_end','AW_Comp','SNPD_ID','TE_chr','TE_start','TE_end','TE_name','TE_class','TE_superfamily','TE_family','TE_Method', 'TE_identity','bp_overlap')

rbound_test <- B73_MultiTE_SNPD_rbound_overlaps %>% 
group_by(AW_ID,AW_Comp) %>% dplyr::mutate(n=n()) %>%
dplyr::mutate(fl_edge_overlap=case_when(
	n==1 ~ TRUE,
	length(unique(TE_family))==1 ~ TRUE,
	TRUE ~ FALSE)) #%>%
#ungroup() %>% group_by(AW_ID,AW_Comp,fl_edge_overlap) %>%
#select(AW_ID,AW_Comp,fl_edge_overlap)%>%
#distinct(AW_ID,AW_Comp,fl_edge_overlap,.keep_all=T)

lbound_join <- lbound_test %>% filter(fl_edge_overlap==TRUE) %>% 
select(AW_ID,AW_start,AW_end,AW_Comp,TE_family,fl_edge_overlap) %>%
distinct()
colnames(lbound_join)[5:6] <- c('l_TE_family','l_fl_edge_overlap')

rbound_join <- rbound_test %>% filter(fl_edge_overlap==TRUE) %>% 
select(AW_ID,AW_start,AW_end,AW_Comp,TE_family,fl_edge_overlap) %>%
distinct()
colnames(rbound_join)[5:6] <- c('r_TE_family','r_fl_edge_overlap')

joined <- left_join(rbound_join,lbound_join,by=c('AW_ID','AW_start','AW_end','AW_Comp'))

B73_MultiTE_SNPD_edge_class <- B73_MultiTE_SNPD_overlaps %>% 
group_by(AW_ID,AW_Comp) %>% dplyr::summarise(edge_case=unique(edge_case))
