# Parse Out Extra Text Around Strings
remove_extra_txt <- function(col_string) {
	rel_text <- unlist(str_split(col_string,pattern='='))[2]
	return(rel_text)
}

# Classify Variants As Either SNPs, Small InDels <= 50 bp, or SV > 50 bp
determine_variant_type <- function(row.df) {
	ID_size <- (as.numeric(row.df[3]) - as.numeric(row.df[2]))+1
	ASM_size <- (as.numeric(row.df[6]) - as.numeric(row.df[5])) +1

	size_dif <- abs(ID_size - ASM_size)

	if (size_dif == 0 ) {
		return('SNP')
	}
	if (size_dif > 0 & size_dif <= 50) {
		return('InDel')
	} 
	if (size_dif > 50) {
		return('structural_variant')
	}
}

update_ASM_columns <- function(gvcf.df,tag) {
	gvcf.df$ASM_CHR <- unlist(lapply(gvcf.df$ASM_CHR,remove_extra_txt))
	gvcf.df$ASM_END <- unlist(lapply(gvcf.df$ASM_END,remove_extra_txt))
	gvcf.df$ASM_Start <- unlist(lapply(gvcf.df$ASM_Start,remove_extra_txt))
	gvcf.df$ASM_Strand <- unlist(lapply(gvcf.df$ASM_Strand,remove_extra_txt))
	if (tag=='nv') {
		gvcf.df$END <- unlist(lapply(gvcf.df$END,remove_extra_txt))	
	}
	return(gvcf.df)
}

make_column_names <- function(ID,ASM) {
	col1 <- paste(ID,'Chr',sep='_')
	col2 <- paste(ID,'StartPos',sep='_')
	col3 <- paste(ID,'EndPos',sep='_')
	col4 <- paste(ASM,'Chr',sep='_')
	col5 <- paste(ASM,'StartPos',sep='_')
	col6 <- paste(ASM,'EndPos',sep='_')	

	col_vec <- c(col1,col2,col3,col4,col5,col6)
	return(col_vec)
}

# We restrict Structural Variants to SVs with size = 0 in one genotype
# and size > 50 bp in the other
reformat_SV_type <- function(row) {
	ID_lineage <- unlist(str_split(names(row)[1],pattern='_'))[1]
	ASM_lineage <- unlist(str_split(names(row)[4],pattern='_'))[1]

	supplied_row <- as.character(row)
	ID_size <- as.numeric(supplied_row[3]) - as.numeric(supplied_row[2])
	ASM_size <- as.numeric(supplied_row[6]) - as.numeric(supplied_row[5])

	if (ID_size > 0 & ASM_size > 0) {
		updated_type <- 'unalignable'
	} else if (ID_size > 0 & ASM_size==0) {
		updated_type <- paste('structural_insertion_in',ID_lineage,sep='')
	} else if (ID_size==0 & ASM_size >0) {
		updated_type <- paste('structural_insertion_in',ASM_lineage,sep='')
	}
	return(updated_type)
}

# There are gaps in the AnchorWave alignment that represent Missing Data
# We identify them so we are aware
identify_missing_blocks <- function(df) {
	missing_blocks <- list()
	for (i in 1:(nrow(df)-1) ) {
		if (df[[i,3]] + 1 != df[[i+1,2]]) {

			id_chr <- df[[i,1]]
			id_start_pos <- df[[i,3]] + 1
			id_end_pos <- df[[i+1,2]] - 1
			asm_chr <- df[[i,4]]
			asm_start_pos <- df[[i,6]] + 1
			asm_end_pos <- df[[i+1,5]]-1
			type <- 'Missing_Data'
			missing_row <- c(id_chr,id_start_pos,id_end_pos,asm_chr,asm_start_pos,asm_end_pos,type)
			missing_blocks <- list.append(missing_blocks,missing_row)
		}
	}

	full_missing.df <- do.call('rbind',missing_blocks) %>% data.frame()
	if (nrow(full_missing.df) != 0) {
		colnames(full_missing.df) <- colnames(df)
		full_missing.df[,1] <- as.integer(full_missing.df[,1])
		full_missing.df[,2] <- as.integer(full_missing.df[,2])
		full_missing.df[,3] <- as.integer(full_missing.df[,3])
		full_missing.df[,4] <- as.integer(full_missing.df[,4])
		full_missing.df[,5] <- as.integer(full_missing.df[,5])
		full_missing.df[,6] <- as.integer(full_missing.df[,6])		
	} 
	return(full_missing.df)
}