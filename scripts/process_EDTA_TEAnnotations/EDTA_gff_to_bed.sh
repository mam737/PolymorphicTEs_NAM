#!/bin/bash

#### EDTA_gff_to_bed.sh
#### Manisha Munasinghe - Last Updated - 01/11/23
#### Take Raw panEDTA TE Annotations
#### and convert from GFF3 Format To Bed Format

raw_dir=/path/to/NAM_EDTA/raw_panEDTA_gff_annotations

for entry in "$raw_dir"/*.gff3
do
  echo "$entry"
#  # Remove gff headers
  #sed -i '/###*/d' *.gff3
  outfile="${entry%%.gff3}"
  outfile+=".bed"
  
  cut -f9 "$entry" | cut -d';' -f1 | cut -d'=' -f2 > temp.bed
  paste $entry temp.bed > temp2.bed
  awk '{print $1,$4,$5,$10,$6,$7,$2,$3,$9}' OFS='\t' temp2.bed > $outfile

done

rm temp.bed
rm temp2.bed

