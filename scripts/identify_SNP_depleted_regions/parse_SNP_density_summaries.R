### parse_SNP_density_summaries.R
### Manisha Munasinghe - Last Updated - 03/03/23
### From sliding window SNP density estimates
### identify SNP depleted regions
.libPaths('/path/to/Rlibs/')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)

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

source('/path/to/function_scripts/SNP_density_functions.R')

NAM_level_orders <- c("B73","B97","Ky21","M162W","Ms71","Oh43","Oh7B","M37W","Mo18W","Tx303","HP301","P39","Il14H","CML52","CML69","CML103","CML228","CML247","CML277","CML322","CML333","Ki3","Ki11","NC350","NC358","Tzi8")
chr_level_orders <- c("chr1","chr2","chr3",'chr4','chr5','chr6','chr7','chr8','chr9','chr10')

#Load + Format Rolling Window SNP Density Across Lineages
SNP_density_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/SNP_density_summaries'
full_SNP_df <- obtain_SNP_densities(SNP_density_dir)

#Identify regions >2Mbp where nonSVbp >= 900,000 + SNP_Count <=100
snp_depleted_regions <- find_SNPdepleted(full_SNP_df)

#We know the SNP depleted coordinates in B73
#(remember the start_intersect is offset by 100,000bp)
#We now need to extract the SNP depleted coordinates in the compared NAM genotype
snp_depleted_regions <- extract_comp_intersect_coords(snp_depleted_regions)

B73_snp_depleted_bed <- make_SNPdepleted_bed(snp_depleted_regions,tag='B73')
NAM_snp_depleted_bed <- make_SNPdepleted_bed(snp_depleted_regions,tag='NAM')

AW_head_dir <- '/path/to/summarised_AnchorWave_Regions'
all_AW_data <- obtain_AnchorWave_files(AW_head_dir)

B73_AW_bed <- make_AW_bed(all_AW_data,tag='B73')
NAM_AW_bed <- make_AW_bed(all_AW_data,tag='NAM')

#Extract AnchorWave Blocks that overlap SNP depleted regions
#and save that information for later use
B73_SNPD_AW_overlaps <- pull_SNPDepleted_AW_overlaps(B73_snp_depleted_bed,B73_AW_bed)
NAM_SNPD_AW_overlaps <- pull_SNPDepleted_AW_overlaps(NAM_snp_depleted_bed,NAM_AW_bed)

fwrite(B73_SNPD_AW_overlaps, '/path/to/summary_data_files/SNP_Depleted_Summaries/B73_SNP_Depleted_AW_Overlaps.csv')
fwrite(NAM_SNPD_AW_overlaps, '/path/to/summary_data_files/SNP_Depleted_Summaries/NAM_SNP_Depleted_AW_Overlaps.csv')

B73_TE_n_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/polymorphic_TE_calls/B73/TE_n'
NAM_TE_n_dir <- '/home/brandvai/mmunasin/TE_Intron/store_data/polymorphic_TE_calls/NAM/TE_n'

all_B73_TE_n_calls <- process_TE_calls(B73_TE_n_dir)
all_NAM_TE_n_calls <- process_TE_calls(NAM_TE_n_dir)

B73_TE_bed <- make_TE_Bed(all_B73_TE_n_calls)
NAM_TE_bed <-make_TE_Bed(all_NAM_TE_n_calls)

#Extract TE annotations that overlap SNP depleted regions
#and save that information for later use
B73_SNPD_TE_overlaps <- pull_SNPDepleted_AW_overlaps(B73_snp_depleted_bed,B73_TE_bed)
NAM_SNPD_TE_overlaps <- pull_SNPDepleted_AW_overlaps(NAM_snp_depleted_bed,NAM_TE_bed)

fwrite(B73_SNPD_TE_overlaps, '/path/to/summary_data_files/SNP_Depleted_Summaries/B73_SNP_Depleted_TE_Overlaps.csv')
fwrite(NAM_SNPD_TE_overlaps, '/path/to/summary_data_files/SNP_Depleted_Summaries/NAM_SNP_Depleted_TE_Overlaps.csv')
