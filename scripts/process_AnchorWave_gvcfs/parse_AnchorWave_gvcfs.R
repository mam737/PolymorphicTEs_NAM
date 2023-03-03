#### parse_AnchorWave_gvcfs.R
#### Manisha Munasinghe - Last Updated - 01/11/23
#### Take Raw AnchorWave GVCF and Bin Entries
#### as either nonvariant, SNP, InDel (<50bp), or Structural Variant (>50bp)
.libPaths('/path/to/Rlibs/')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)

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
#[1] tidyr_1.2.0       data.table_1.14.2 stringr_1.4.0     dplyr_1.0.9#

#loaded via a namespace (and not attached):
# [1] fansi_1.0.3      assertthat_0.2.1 utf8_1.2.2       crayon_1.5.1
# [5] R6_2.5.1         DBI_1.1.2        lifecycle_1.0.1  magrittr_2.0.3
# [9] pillar_1.7.0     stringi_1.7.6    rlang_1.0.2      cli_3.3.0
#[13] vctrs_0.4.1      generics_0.1.2   ellipsis_0.3.2   tools_4.0.4
#[17] glue_1.6.2       purrr_0.3.4      compiler_4.0.4   pkgconfig_2.0.3
#[21] tidyselect_1.1.2 tibble_3.1.6

source('/path/to/function_scripts/AnchorWave_process_functions.R')

# Pass in AnchorWave GVCF Pairwise Alignment Between B73 and a NAM line
# Remove Headers from AnchorWave GVCF prior to this step
args <- commandArgs(trailingOnly = TRUE)
gvcf_filename <- args[1]

ID_lineage <- 'B73'
ASM_lineage <- unlist(str_split(unlist(str_split(gvcf_filename,pattern='/'))[8],pattern='_'))[1]

gvcf <- fread(gvcf_filename)
names(gvcf)[1] <- 'CHROM'

# If the string "END" is in the INFO column, that means that entry is
# a nonvariant entry between the two genotypes
nonvariant <- gvcf[grep("END",gvcf$INFO),]

nonvariant <- nonvariant %>% 
tidyr::separate(col=INFO,into=c('ASM_CHR','ASM_END','ASM_Start','ASM_Strand', 'END'),sep=';')

# Split the attribute string
nonvariant <- update_ASM_columns(nonvariant,tag='nv')

nv_regions <- data.frame(nonvariant$CHROM,nonvariant$POS,nonvariant$END,nonvariant$ASM_CHR,nonvariant$ASM_Start,nonvariant$ASM_END)

df_col_names <- make_column_names(ID_lineage,ASM_lineage)

colnames(nv_regions) <- df_col_names
nv_regions$Type <- 'nonvariant_region'

# All the remaining regions are variant regions classified as either SNP, InDel, or SV
variant <- gvcf[!grep("END",gvcf$INFO),]
variant$ref_bp_size <- unlist(lapply(variant$REF,nchar))
variant$REF_END <- (variant$POS + variant$ref_bp_size)-1

variant <- variant %>% 
tidyr::separate(col=INFO,into=c('ASM_CHR','ASM_END','ASM_Start','ASM_Strand'),sep=';')

variant <- update_ASM_columns(variant,tag='v')

v_regions <- data.frame(variant$CHROM,variant$POS,variant$REF_END, variant$ASM_CHR,variant$ASM_Start, variant$ASM_END)

colnames(v_regions) <- df_col_names

# Run for Oh43
# v_regions <- v_regions %>% filter(B73_StartPos!=234530060)
v_regions$Type <- apply(v_regions,1,determine_variant_type)

full_region <- rbind(nv_regions,v_regions) %>% arrange(.[[1]],.[[2]])

output_dir <- '/path/to/output_dir/Parsed_AnchorWave_Regions/'
output_subdir <- paste(ID_lineage,ASM_lineage,'regions',sep='_')
final_output_dir <- paste(output_dir,output_subdir,sep='')
dir.create(final_output_dir)

# Write Parsed AnchorWave 
for (i in unique(full_region[[1]])) {
	sub_full_region <- full_region %>% subset(.[[1]]==i)
	chr_header <- paste('chr',i,sep='')
	filename <- paste(chr_header,ID_lineage,ASM_lineage,'regions.tsv',sep='_')

	output_path <- paste(final_output_dir,filename,sep='/')
	write.table(sub_full_region,output_path,quote=F,col.names=T,row.names=F)
	
}
