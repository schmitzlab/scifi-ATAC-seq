#!/bin/bash

bedfiles=

# functions
runQC(){

	# input
	bed=scifi-ATAC-seq.tn5.bed.gz
	id=$( basename $bed )
	name=$( echo $id | cut -d'.' -f1 )
	ref=~/GenomeRef/Zea_mays_v5/v5/Zm-B73-REFERENCE-NAM-5.0.fa.fai.sorted
	ann=~/GenomeRef/Zea_mays_v5/v5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.sorted.gff3
	
	# verbose
	echo " - running Socrates scATAC-seq QC analysis for $name ..."
	# run 
	Rscript ./QC_scATAC_data_peaks.R $bed $name $ann $ref
	Rscript ./filter_low_qual_cells.R $name.raw.soc.rds $name
    Rscript ./scATAC_meta_QC-xz_fixtn5-1k.R ${name}.updated_metadata.txt $name
}
export -f runQC

parallel -j $threads runQC {} ::: $bedfiles
