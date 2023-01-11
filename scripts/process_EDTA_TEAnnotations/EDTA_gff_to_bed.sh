#!/bin/bash

raw_dir=/home/brandvai/mmunasin/NAM_EDTA/raw_gff_files

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

