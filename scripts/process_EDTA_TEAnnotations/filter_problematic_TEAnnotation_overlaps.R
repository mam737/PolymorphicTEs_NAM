#### filter_problematic_TEAnnotation_overlaps.R
#### Manisha Munasinghe - Last Updated - 01/11/23
#### Certain overlaps between TE Annotations
#### Reflect situations that are biologically unfeasible
#### and should consequently be removed from TE Annotation File

.libPaths('/path/to/Rlibs/')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(ggplot2)
# Must Run Before Using bedtoolsr
# export PATH=/path/to/bedtools2/bin/:$PATH 
# Using bedtools v2.30.0
library(bedtoolsr) 
library(ggpubr)
library(grid)
set.seed(1994)

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
#[1] grid      stats     graphics  grDevices utils     datasets  methods
#[8] base#

#other attached packages:
#[1] ggpubr_0.4.0       bedtoolsr_2.30.0-4 ggplot2_3.3.5      tidyr_1.2.0
#[5] data.table_1.14.2  stringr_1.4.0      dplyr_1.0.9#

#loaded via a namespace (and not attached):
# [1] pillar_1.7.0     compiler_4.0.4   tools_4.0.4      lifecycle_1.0.1
# [5] tibble_3.1.6     gtable_0.3.0     pkgconfig_2.0.3  rlang_1.0.2
# [9] DBI_1.1.2        cli_3.3.0        withr_2.5.0      generics_0.1.2
#[13] vctrs_0.4.1      tidyselect_1.1.2 glue_1.6.2       R6_2.5.1
#[17] rstatix_0.7.0    fansi_1.0.3      carData_3.0-5    purrr_0.3.4
#[21] car_3.0-12       magrittr_2.0.3   scales_1.2.0     backports_1.4.1
#[25] ellipsis_0.3.2   assertthat_0.2.1 abind_1.4-5      colorspace_2.0-3
#[29] ggsignif_0.6.3   utf8_1.2.2       stringi_1.7.6    munsell_0.5.0
#[33] broom_0.8.0      crayon_1.5.1

# Start with Initially Processed TE Annotation with nonTE annotations removed
# These files have been split by chromosome + are stored within one directory
args <- commandArgs(trailingOnly = TRUE)
TE_dir_name <- args[1]

lineage <- unlist(str_split(TE_dir_name,pattern='/'))[8]

# Combine all TE files across chromosomes 
TE_files <- list.files(TE_dir_name,full.names=T)
TE_data <- lapply(TE_files,fread)
TE_data <- do.call('rbind',TE_data)
TE_data <- TE_data %>% separate(chr,into=c("lineage",'chr'),sep='_')

# Store TE data in appropriate bed format
TE_bed <- TE_data[,c(2,3,4,5,9,11,13,14)]
all_TE_names <- unique(TE_bed$name)

class_I <- c("LTR/Ty3","LTR/Copia","LTR/unknown")
class_II <- c("DTM/MITE","DTA/DNA","DTM/DNA","DTC/DNA","DTA/MITE","DTH/MITE","DTT/MITE","DTC/MITE","DTT/DNA","DTH/DNA")

starting_summary_method <- TE_bed %>% mutate(size=end-start) %>% group_by(method) %>% 
summarise(n=n(),Mb=sum(size)/1000000) %>% rbind(c('Total',sum(.$n),sum(.$Mb))) %>% 
pivot_wider(names_from=method,values_from=c(n,Mb),names_glue='{method}_{.value}',names_vary='slowest') %>%
select(Total_n,Total_Mb,structural_n,structural_Mb,homology_n,homology_Mb) %>% data.frame()

starting_summary_class <- TE_bed %>% mutate(size=end-start) %>% 
mutate(class = case_when(raw_superfamily %in% class_I ~ 'classI', TRUE ~'classII')) %>% 
group_by(class) %>% summarise(n=n(),Mb=sum(size)/1000000) %>% 
pivot_wider(names_from=class,values_from=c(n,Mb),names_glue='{class}_{.value}',names_vary='slowest') %>%
select(classI_n,classI_Mb,classII_n,classII_Mb) %>% data.frame()

starting_summary <- cbind(starting_summary_method,starting_summary_class) %>% mutate(Step='Step0',.before=Total_n)

###### STEP 1 - Remove all TEs that are "identical" - same size, supefamily, family classification ######

overlap_prop <- bedtoolsr::bt.intersect(TE_bed,TE_bed,wo=T)
colnames(overlap_prop) <- c('a_chr','a_start','a_end','a_name','a_type','a_superfamily','a_family','a_method','b_chr','b_start','b_end','b_name','b_type','b_superfamily','b_family','b_method','bp_overlap')
overlap_prop <- overlap_prop %>% 
mutate(a_size=a_end-a_start,.before='a_name') %>% 
mutate(b_size=b_end-b_start,.before='b_name') %>% 
subset(a_name!=b_name) %>% 
mutate(a_prop_overlap=bp_overlap/a_size,b_prop_overlap=bp_overlap/b_size)

identical_match <- overlap_prop %>% 
subset(a_size == b_size) %>% 
subset(a_chr==b_chr) %>% 
subset(a_start==b_start) %>% subset(a_end==b_end)

## First, let's match by any structural + homology pair, we will always retain the structural pair
remove_TE_list <- c()
for (i in unique(identical_match$a_name)) {
	sub_overlap <- identical_match %>% subset(a_name==i)

	if (sub_overlap$a_method=='structural' & sub_overlap$b_method=='homology') {
		remove_TE_list <- c(remove_TE_list,sub_overlap$b_name)
	} else if (sub_overlap$a_method=='homology' & sub_overlap$b_method=='structural') {
		remove_TE_list <- c(remove_TE_list,sub_overlap$a_name)

	} else {
		# We do this so we don't accidently randomly sample two different ones which would keep the duplicate in
		if (!(sub_overlap$a_name %in% remove_TE_list) & !(sub_overlap$b_name %in% remove_TE_list) ) {	
		remove_TE_list <- c(remove_TE_list, sample(c(sub_overlap$a_name,sub_overlap$b_name),1) ) 
		}
	}
}

unique_TE_bed <- TE_bed %>% subset(!(name %in% remove_TE_list))
all_TE_names <- all_TE_names[!(all_TE_names %in% remove_TE_list)]

step1_summary_method <- unique_TE_bed %>% mutate(size=end-start) %>% group_by(method) %>% 
summarise(n=n(),Mb=sum(size)/1000000) %>% rbind(c('Total',sum(.$n),sum(.$Mb))) %>% 
pivot_wider(names_from=method,values_from=c(n,Mb),names_glue='{method}_{.value}',names_vary='slowest') %>%
select(Total_n,Total_Mb,structural_n,structural_Mb,homology_n,homology_Mb) %>% data.frame()

step1_summary_class <- unique_TE_bed %>% mutate(size=end-start) %>% 
mutate(class = case_when(raw_superfamily %in% class_I ~ 'classI', TRUE ~'classII')) %>% 
group_by(class) %>% summarise(n=n(),Mb=sum(size)/1000000) %>% 
pivot_wider(names_from=class,values_from=c(n,Mb),names_glue='{class}_{.value}',names_vary='slowest') %>%
select(classI_n,classI_Mb,classII_n,classII_Mb) %>% data.frame()

step1_summary <- cbind(step1_summary_method,step1_summary_class) %>% mutate(Step='Step1',.before=Total_n)

###### STEP 2 - Extract all distinct TEs that do not overlap another TE! ######

overlap_prop <- bedtoolsr::bt.intersect(unique_TE_bed,unique_TE_bed,wo=T)
colnames(overlap_prop) <- c('a_chr','a_start','a_end','a_name','a_type','a_superfamily','a_family','a_method','b_chr','b_start','b_end','b_name','b_type','b_superfamily','b_family','b_method','bp_overlap')
overlap_prop <- overlap_prop %>% mutate(a_size=a_end-a_start,.before='a_name') %>% mutate(b_size=b_end-b_start,.before='b_name') %>% subset(a_name!=b_name) %>% mutate(a_prop_overlap=bp_overlap/a_size,b_prop_overlap=bp_overlap/b_size)

overlap_TE_names <- unique(overlap_prop$a_name)

distinct_TE_names <- all_TE_names[!(all_TE_names %in% overlap_TE_names)]

overlap_TE_bed <- TE_bed %>% subset(name %in% overlap_TE_names)

distinct_TE_bed <- TE_bed %>% subset(name %in% distinct_TE_names)

###### STEP 3 - Extract Structural TEs + Remove Edge Cases + Homology Overlaps ######
## 'Edge-Case' overlaps as an annotation whose start position is within 5bp of the
##   start of a structural TE annotation or similarly
##   whose end position is within 5bp of the end of a structural TE annotation

structural_overlap_TE_bed <- overlap_TE_bed %>% subset(method=='structural')

structural_overlap_prop <- bedtoolsr::bt.intersect(structural_overlap_TE_bed,unique_TE_bed,wo=T)
colnames(structural_overlap_prop) <- c('a_chr','a_start','a_end','a_name','a_type','a_superfamily','a_family','a_method','b_chr','b_start','b_end','b_name','b_type','b_superfamily','b_family','b_method','bp_overlap')
structural_overlap_prop <- structural_overlap_prop %>% mutate(a_size=a_end-a_start,.before='a_name') %>% mutate(b_size=b_end-b_start,.before='b_name') %>% subset(a_name!=b_name) %>% mutate(a_prop_overlap=bp_overlap/a_size,b_prop_overlap=bp_overlap/b_size)

start_edge_case_overlap <- structural_overlap_prop %>% subset(between(b_start,a_start,a_start+5))

start_remove_edge_case <- c()

for (i in 1:nrow(start_edge_case_overlap)) {
	sub_overlap <- start_edge_case_overlap[i,]

	# By definition, a must be structural while either b can be either structural/homology
	if (sub_overlap$a_method == 'structural' & sub_overlap$b_method == 'homology') {
		start_remove_edge_case <- c(start_remove_edge_case,sub_overlap$b_name)
	} else if (sub_overlap$a_size > sub_overlap$b_size) { # From here out, a is structural and b is structural
		start_remove_edge_case  <- c(start_remove_edge_case,sub_overlap$b_name) # remove smaller
	} else if (sub_overlap$a_size < sub_overlap$b_size) {
		start_remove_edge_case <- c(start_remove_edge_case,sub_overlap$a_name)
	} else { # if equally sized structural, randomly select one with care to not double sample as you could then retain overlaps
		if (!(sub_overlap$a_name %in% start_remove_edge_case) & !(sub_overlap$b_name %in% start_remove_edge_case)) {
			start_remove_edge_case <- c(start_remove_edge_case, sample(c(sub_overlap$a_name,sub_overlap$b_name),1) ) 
		}
	}
}

end_edge_case_overlap <- structural_overlap_prop %>% subset(between(b_end,a_end-5,a_end))

end_remove_edge_case <- c()

for (i in 1:nrow(end_edge_case_overlap)) {
	sub_overlap <- end_edge_case_overlap[i,]

	if (sub_overlap$a_method == 'structural' & sub_overlap$b_method == 'homology') {
		end_remove_edge_case <- c(end_remove_edge_case,sub_overlap$b_name)
	} else if (sub_overlap$a_size > sub_overlap$b_size) { # All subsequent cases are whe  a is structural and b is structural
		end_remove_edge_case  <- c(end_remove_edge_case,sub_overlap$b_name)
	} else if (sub_overlap$a_size < sub_overlap$b_size) {
		end_remove_edge_case <- c(end_remove_edge_case,sub_overlap$a_name)
	} else {
		if (!(sub_overlap$a_name %in% end_remove_edge_case) & !(sub_overlap$b_name %in% end_remove_edge_case)) {
			end_remove_edge_case <- c(end_remove_edge_case, sample(c(sub_overlap$a_name,sub_overlap$b_name),1) ) 
		}
	}
}

remove_edge_case <- unique(c(start_remove_edge_case,end_remove_edge_case))

current_retain_TE_bed <- unique_TE_bed %>% subset(!(name %in% remove_edge_case))

current_retain_overlap_TE_bed <- overlap_TE_bed %>% subset(!(name %in% remove_edge_case))

step3_summary_method <- current_retain_TE_bed %>% mutate(size=end-start) %>% group_by(method) %>% 
summarise(n=n(),Mb=sum(size)/1000000) %>% rbind(c('Total',sum(.$n),sum(.$Mb))) %>% 
pivot_wider(names_from=method,values_from=c(n,Mb),names_glue='{method}_{.value}',names_vary='slowest') %>%
select(Total_n,Total_Mb,structural_n,structural_Mb,homology_n,homology_Mb) %>% data.frame()

step3_summary_class <- current_retain_TE_bed %>% mutate(size=end-start) %>% 
mutate(class = case_when(raw_superfamily %in% class_I ~ 'classI', TRUE ~'classII')) %>% 
group_by(class) %>% summarise(n=n(),Mb=sum(size)/1000000) %>% 
pivot_wider(names_from=class,values_from=c(n,Mb),names_glue='{class}_{.value}',names_vary='slowest') %>%
select(classI_n,classI_Mb,classII_n,classII_Mb) %>% data.frame()

step3_summary <- cbind(step3_summary_method,step3_summary_class) %>% mutate(Step='Step3',.before=Total_n)

###### STEP 4 - Remove homology TEs fully contained within another homolgoy TEs ######

retained_homo_TE_bed <- current_retain_overlap_TE_bed  %>% subset(method=='homology')

homo_fully_contained <- bedtoolsr::bt.intersect(retained_homo_TE_bed,retained_homo_TE_bed,wo=T,F=1.0)
colnames(homo_fully_contained) <- c('a_chr','a_start','a_end','a_name','a_type','a_superfamily','a_family','a_method','b_chr','b_start','b_end','b_name','b_type','b_superfamily','b_family','b_method','bp_overlap')
homo_fully_contained <- homo_fully_contained %>% mutate(a_size=a_end-a_start,.before='a_name') %>% mutate(b_size=b_end-b_start,.before='b_name') %>% subset(a_name!=b_name) %>% mutate(a_prop_overlap=bp_overlap/a_size,b_prop_overlap=bp_overlap/b_size)

remove_fully_contain_homo <- unique(homo_fully_contained$b_name)

current_retain_TE_bed <- current_retain_TE_bed %>% subset(!(name %in% remove_fully_contain_homo))

current_retain_overlap_TE_bed <- current_retain_overlap_TE_bed %>% subset(!(name %in% remove_fully_contain_homo))

step4_summary_method <- current_retain_TE_bed %>% mutate(size=end-start) %>% group_by(method) %>% 
summarise(n=n(),Mb=sum(size)/1000000) %>% rbind(c('Total',sum(.$n),sum(.$Mb))) %>% 
pivot_wider(names_from=method,values_from=c(n,Mb),names_glue='{method}_{.value}',names_vary='slowest') %>%
select(Total_n,Total_Mb,structural_n,structural_Mb,homology_n,homology_Mb) %>% data.frame()

step4_summary_class <- current_retain_TE_bed %>% mutate(size=end-start) %>% 
mutate(class = case_when(raw_superfamily %in% class_I ~ 'classI', TRUE ~'classII')) %>% 
group_by(class) %>% summarise(n=n(),Mb=sum(size)/1000000) %>% 
pivot_wider(names_from=class,values_from=c(n,Mb),names_glue='{class}_{.value}',names_vary='slowest') %>%
select(classI_n,classI_Mb,classII_n,classII_Mb) %>% data.frame()

step4_summary <- cbind(step4_summary_method,step4_summary_class) %>% mutate(Step='Step4',.before=Total_n)

###### STEP 5 - Remove homology TEs with greater than 10% overlap with another homology TE ######

retained_homo_TE_bed <- current_retain_overlap_TE_bed  %>% subset(method=='homology')

homo_partial_overlap <- bedtoolsr::bt.intersect(retained_homo_TE_bed,retained_homo_TE_bed,wo=T)
colnames(homo_partial_overlap) <- c('a_chr','a_start','a_end','a_name','a_type','a_superfamily','a_family','a_method','b_chr','b_start','b_end','b_name','b_type','b_superfamily','b_family','b_method','bp_overlap')
homo_partial_overlap  <- homo_partial_overlap  %>% mutate(a_size=a_end-a_start,.before='a_name') %>% mutate(b_size=b_end-b_start,.before='b_name') %>% subset(a_name!=b_name) %>% mutate(a_prop_overlap=bp_overlap/a_size,b_prop_overlap=bp_overlap/b_size)

greater_10_overlap_TEs <- homo_partial_overlap %>% subset(b_prop_overlap > 0.1) %>% pull(b_name) %>% unique()

current_retain_TE_bed <- current_retain_TE_bed %>% subset(!(name %in% greater_10_overlap_TEs))

current_retain_overlap_TE_bed <- current_retain_overlap_TE_bed %>% subset(!(name %in% greater_10_overlap_TEs))

step5_summary_method <- current_retain_TE_bed %>% mutate(size=end-start) %>% group_by(method) %>% 
summarise(n=n(),Mb=sum(size)/1000000) %>% rbind(c('Total',sum(.$n),sum(.$Mb))) %>% 
pivot_wider(names_from=method,values_from=c(n,Mb),names_glue='{method}_{.value}',names_vary='slowest') %>%
select(Total_n,Total_Mb,structural_n,structural_Mb,homology_n,homology_Mb) %>% data.frame()

step5_summary_class <- current_retain_TE_bed %>% mutate(size=end-start) %>% 
mutate(class = case_when(raw_superfamily %in% class_I ~ 'classI', TRUE ~'classII')) %>% 
group_by(class) %>% summarise(n=n(),Mb=sum(size)/1000000) %>% 
pivot_wider(names_from=class,values_from=c(n,Mb),names_glue='{class}_{.value}',names_vary='slowest') %>%
select(classI_n,classI_Mb,classII_n,classII_Mb) %>% data.frame()

step5_summary <- cbind(step5_summary_method,step5_summary_class) %>% mutate(Step='Step5',.before=Total_n)

###### STEP 6 - Remove homology TEs with partial but not complete overlap with a structural TE ######
### if > thatn 5% of the homology TE annotation overlapped a structural annotation
### filter it out
retained_homo_TE_bed <- current_retain_overlap_TE_bed  %>% subset(method=='homology')
retained_struct_TE_bed <- current_retain_overlap_TE_bed  %>% subset(method=='structural')

struct_homo_partial_overlap <- bedtoolsr::bt.intersect(retained_struct_TE_bed,retained_homo_TE_bed,wo=T) # a is structural, b is homology
colnames(struct_homo_partial_overlap) <- c('a_chr','a_start','a_end','a_name','a_type','a_superfamily','a_family','a_method','b_chr','b_start','b_end','b_name','b_type','b_superfamily','b_family','b_method','bp_overlap')
struct_homo_partial_overlap  <- struct_homo_partial_overlap  %>% mutate(a_size=a_end-a_start,.before='a_name') %>% mutate(b_size=b_end-b_start,.before='b_name') %>% subset(a_name!=b_name) %>% mutate(a_prop_overlap=bp_overlap/a_size,b_prop_overlap=bp_overlap/b_size)

remove_struct_homo_partial <- struct_homo_partial_overlap %>% subset(b_prop_overlap!=1.0) %>% subset(b_prop_overlap > 0.05) %>% pull(b_name) %>% unique()

current_retain_TE_bed <- current_retain_TE_bed %>% subset(!(name %in% remove_struct_homo_partial))

current_retain_overlap_TE_bed <- current_retain_overlap_TE_bed %>% subset(!(name %in% remove_struct_homo_partial))

step6_summary_method <- current_retain_TE_bed %>% mutate(size=end-start) %>% group_by(method) %>% 
summarise(n=n(),Mb=sum(size)/1000000) %>% rbind(c('Total',sum(.$n),sum(.$Mb))) %>% 
pivot_wider(names_from=method,values_from=c(n,Mb),names_glue='{method}_{.value}',names_vary='slowest') %>%
select(Total_n,Total_Mb,structural_n,structural_Mb,homology_n,homology_Mb) %>% data.frame()

step6_summary_class <- current_retain_TE_bed %>% mutate(size=end-start) %>% 
mutate(class = case_when(raw_superfamily %in% class_I ~ 'classI', TRUE ~'classII')) %>% 
group_by(class) %>% summarise(n=n(),Mb=sum(size)/1000000) %>% 
pivot_wider(names_from=class,values_from=c(n,Mb),names_glue='{class}_{.value}',names_vary='slowest') %>%
select(classI_n,classI_Mb,classII_n,classII_Mb) %>% data.frame()

step6_summary <- cbind(step6_summary_method,step6_summary_class) %>% mutate(Step='Step6',.before=Total_n)

###### STEP 7 - Remove structural TEs with partial but not complete overlap with a structural TE ######
## Remove structural overlaps that overlap another structural annotation
## but are not completely contained (i.e., we retain nested annotations)
retained_struct_TE_bed <- current_retain_overlap_TE_bed  %>% subset(method=='structural')

struct_struct_partial_overlap <- bedtoolsr::bt.intersect(retained_struct_TE_bed,retained_struct_TE_bed,wo=T)
colnames(struct_struct_partial_overlap) <- c('a_chr','a_start','a_end','a_name','a_type','a_superfamily','a_family','a_method','b_chr','b_start','b_end','b_name','b_type','b_superfamily','b_family','b_method','bp_overlap')
struct_struct_partial_overlap  <- struct_struct_partial_overlap  %>% mutate(a_size=a_end-a_start,.before='a_name') %>% mutate(b_size=b_end-b_start,.before='b_name') %>% subset(a_name!=b_name) %>% mutate(a_prop_overlap=bp_overlap/a_size,b_prop_overlap=bp_overlap/b_size)

struct_struct_partial_overlap <- struct_struct_partial_overlap %>% subset(b_prop_overlap !=1.0)

remove_struct_struct_partial <- c()

for (i in 1:nrow(struct_struct_partial_overlap)) {
	sub_overlap <- struct_struct_partial_overlap[i,]


	if (sub_overlap$a_superfamily %in% class_I & sub_overlap$b_superfamily %in% class_II) {
		remove_struct_struct_partial <- c(remove_struct_struct_partial,sub_overlap$b_name)
	} else if (sub_overlap$a_superfamily %in% class_II & sub_overlap$b_superfamily %in% class_I) {
		remove_struct_struct_partial <- c(remove_struct_struct_partial,sub_overlap$a_name)
	} else if (sub_overlap$a_size > sub_overlap$b_size) { # <- At this point, we know they must be the same class
		remove_struct_struct_partial <- c(remove_struct_struct_partial,sub_overlap$b_name)
	} else if (sub_overlap$a_size < sub_overlap$b_size) {
		remove_struct_struct_partial <- c(remove_struct_struct_partial,sub_overlap$a_name)
	} else {
		if (!(sub_overlap$a_name %in% remove_struct_struct_partial ) & !(sub_overlap$b_name %in% remove_struct_struct_partial )) {
			remove_struct_struct_partial  <- c(remove_struct_struct_partial , sample(c(sub_overlap$a_name,sub_overlap$b_name),1) ) 
		}		
	}
}

final_retain_TE_bed <- current_retain_TE_bed %>% subset(!(name %in% remove_struct_struct_partial))

final_summary_method <- final_retain_TE_bed %>% mutate(size=end-start) %>% group_by(method) %>% 
summarise(n=n(),Mb=sum(size)/1000000) %>% rbind(c('Total',sum(.$n),sum(.$Mb))) %>% 
pivot_wider(names_from=method,values_from=c(n,Mb),names_glue='{method}_{.value}',names_vary='slowest') %>%
select(Total_n,Total_Mb,structural_n,structural_Mb,homology_n,homology_Mb) %>% data.frame()

final_summary_class <- final_retain_TE_bed %>% mutate(size=end-start) %>% 
mutate(class = case_when(raw_superfamily %in% class_I ~ 'classI', TRUE ~'classII')) %>% 
group_by(class) %>% summarise(n=n(),Mb=sum(size)/1000000) %>% 
pivot_wider(names_from=class,values_from=c(n,Mb),names_glue='{class}_{.value}',names_vary='slowest') %>%
select(classI_n,classI_Mb,classII_n,classII_Mb) %>% data.frame()

final_summary <- cbind(final_summary_method,final_summary_class) %>% mutate(Step='Step7',.before=Total_n)

complete_summary.df <- rbind(starting_summary, step1_summary,step3_summary,step4_summary,step5_summary,step6_summary,final_summary) %>% mutate(Lineage=lineage,.before=Step)

outfile_name <- paste('/path/to/filtered_TE_anno/summary_files/',lineage,'_filtering_summary_outputs.txt',sep='')
write.table(complete_summary.df,file=outfile_name,quote=F,sep='\t',col.names=T,row.names=F)

final_TE_full <- TE_data %>% subset(name %in% final_retain_TE_bed$name) %>% 
mutate(class = case_when(raw_superfamily %in% class_I ~ 'classI', TRUE ~'classII')) %>%
mutate(condense_superfamily = case_when(
	raw_superfamily == 'LTR/Copia' ~ 'RLC',
	raw_superfamily == 'LTR/Ty3' ~ 'RLG',
	raw_superfamily == 'LTR/unknown' ~ 'RLX',
	raw_superfamily %in% c('LINE/L1','LINE/RTE','LINE/unknown') ~ 'RIL',
	raw_superfamily %in% c('DTA/MITE','DTA/DNA') ~ 'DTA',
	raw_superfamily %in% c('DTC/MITE','DTC/DNA') ~ 'DTC',
	raw_superfamily %in% c('DTH/MITE','DTH/DNA') ~ 'DTH',
	raw_superfamily %in% c('DTM/MITE','DTM/DNA') ~ 'DTM',
	raw_superfamily %in% c('DTT/MITE','DTT/DNA') ~ 'DTT' ))

dir.create(paste('/path/to/filtered_TE_anno/',lineage,sep=''))

for (i in unique(final_TE_full$chr)) {
	chr_specific <- final_TE_full %>% subset(chr==i)
	file_path <- paste('/path/to/filtered_TE_anno/',lineage,'/',lineage,'_',i,'_EDTA_filtered.bed',sep='')
	write.table(chr_specific,file=file_path,quote=F,sep='\t',col.names=T,row.names=F)
}