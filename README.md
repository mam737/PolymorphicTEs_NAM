# Combined analysis of transposable elements and structural variation in maize genomes reveals genome contraction outpaces expansion

This repository contains scripts used to identify polymorphic TEs between B73 and
the remaining 25 maize inbred founder lines for the maize Nested Association Mapping (NAM) population. It also contains all several relevant datasets generated and used in this analysis that are available here for public use. 

Subfolders within the scripts folder correlate to specific steps of the analysis.

1. proccess_EDTA_TEAnnotations - Scripts used to filter and perform QC on the raw 
panEDTA TE Annotations generated in [Ou et al.](https://www.biorxiv.org/content/10.1101/2022.10.09.511471v1) (bioRxiv 2022). These TE Annotations are publicly available
on [MaizeGDB](https://maizegdb.org/NAM_project)
	- EDTA_gff_to_bed.sh: Shell script to change GFF files to BED files
	- process_initial_EDTA_bed.R: Remove annotations for nonTE, helitron, and specific features of structurally annotated LTRs
	- filter_problematic_TEAnnotation_overlaps.R: Filter overlaps between TE annotations that reflect situations that are biologically unfeasible and should consequently be filtered from TE Annotation File

2. process_AnchorWave_gvcfs - Scripts used to parse AnchorWave gvcf pairwise alignments into alignable, structural variant, and unalignable sequence
	- parse_AnchorWave_gvcfs.R: Parse AnchorWave GVCF pairwise alignment into nonvariant, SNP, InDel (<50bp), and structural variant (>50bp) sequence
	- generated_summarisedAW.R: Using parsed AnchorWave outputs bin regions into either alignable (nonvariant, SNP, and Indel), structural variant (>50 bp in one genome, 0 bp in other genome), and unalignable sequence

3. classify_feature_annotations - Scripts used to classify features by intersecting annotations with summarised AnchorWave alignments
	- classify_polymorphic_TEAnnotations.R: Classify TE Annotations
	- classify_polymorphic_geneAnnotations.R: Classify Gene Annotations (both exon-only and full-length)