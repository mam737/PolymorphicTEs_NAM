# For a given feature, extract which AW Blocks it overlaps
pull_AW_blocks <- function(id_name_val,bed_df,comp) {
	if (comp=='id') {
		AW_Blocks <- bed_df %>% subset(id_name==id_name_val) %>% 
		pull(AW_name) %>% unique()
	}
	if (comp=='asm') {
		AW_Blocks <- bed_df %>% subset(asm_name==id_name_val) %>% 
		pull(AW_name) %>% unique()
	}
	AW_Blocks <- paste(AW_Blocks,collapse=',')
	return(AW_Blocks)
}

# There are technically 5 possibly groupings
# alignable, structural variation in B73, structural variation in NAM, unalignable, or missing data
# ensure column for each of these groupings for consistency across all files
check_missing_cols <- function(final_df,ID_struct_tag,ASM_struct_tag) {	
	if (!('alignable_region' %in% colnames(final_df))) {
		final_df$alignable_region <- 0
	}
	if (!(ID_struct_tag %in% colnames(final_df))) {
		final_df[[ID_struct_tag]] <- 0
	}
	if (!(ASM_struct_tag %in% colnames(final_df))) {
		final_df[[ASM_struct_tag]] <- 0
	}
	if(!('unalignable' %in% colnames(final_df))) {
		final_df$unalignable <- 0
	}
	if(!('Missing_Data' %in% colnames(final_df))) {
		final_df$Missing_Data <- 0
	}
	return(final_df)
}

# Obtain full length gene annotation
generate_full_length <- function(gff.df) {
	start_gff.df <- gff.df %>% group_by(Canonical_Tx) %>% slice(1) %>% select(Chr,Start,Canonical_Tx)
	end_gff.df <- gff.df %>% group_by(Canonical_Tx) %>% slice(n()) %>% select(Chr,End,Canonical_Tx)

	full_gff.df <- cbind(start_gff.df,End=end_gff.df$End)
	full_gff.df <- full_gff.df[,c(1,2,4,3)]
	return(full_gff.df)
}

# Extract exonID
extract_exon_classification <- function(attribute_string) {
  classification <- substring(attribute_string, regexpr("*exon_id=", attribute_string),regexpr(";r", attribute_string)-1)
  classification <- unlist(str_split(classification,pattern='='))[2]
  return(classification)
}

# Pull Gene Start Coord
add_gene_start <- function(id_name_val,gff_exon) {
	gene_start <- gff_exon %>% subset(Canonical_Tx==id_name_val) %>% slice(1) %>% select(Start)
	return(gene_start)
}

# Pull Gene End Coord
add_gene_end <- function(id_name_val,gff_exon) {
	gene_end <- gff_exon %>% subset(Canonical_Tx==id_name_val) %>% slice(n()) %>% select(End)
	return(gene_end)
}

# Pull Exon Number for Exon in Gene
add_gene_exon_number <- function(id_name_val,gff_exon) {
	exon_number <- gff_exon %>% subset(Canonical_Tx==id_name_val) %>% nrow()
	return(exon_number)
}
