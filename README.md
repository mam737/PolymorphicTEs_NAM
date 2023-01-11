# Combined analysis of transposable elements and structural variation in maize genomes reveals genome contraction outpaces expansion

This repository contains scripts used to identify polymorphic TEs between B73 and
the remaining 25 maize inbred founder lines for the maize Nested Association Mapping (NAM)
population.

Subfolders within the scripts folder correlate to specific steps of the analysis.

1. proccess_EDTA_TEAnnotations - Scripts used to filter and perform QC on the raw 
panEDTA TE Annotations generated in [Ou et al.](https://www.biorxiv.org/content/10.1101/2022.10.09.511471v1) (bioRxiv 2022). These TE Annotations are publicly available
on [MaizeGDB](https://maizegdb.org/NAM_project)
	- EDTA_gff_to_bed.sh: Shell script to change GFF files to BED files
	- process_initial_EDTA_bed.R: Remove annotations for nonTE, helitron, and specific features of structurally annotated LTRs
	- filter_problematic_TEAnnotation_overlaps.R: Filter overlaps between TE annotations that reflect situations that are biologically unfeasible and should consequently be filtered from TE Annotation File
