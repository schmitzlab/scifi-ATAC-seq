#!/usr/bin/bash

#load modules
runSouporcell(){
array=($1)
sample=${array[0]}
barcode=${array[1]}
fasta=~/reference_genomes/Zmays/v5/Zm-B73-REFERENCE-NAM-5.0.fa
vcf=./NAM_7geno.vcf

# inputs
bam=./BAM/scifi-merged.mq10.BC.rmdup.mm.CB.bam
meta=./scifi-merged.updated_metadata_v5.txt
min=50
threads=8
id=$(basename $meta .txt)_${sample}_mindp${min}_${barcode}

# functions
	#id=$1
	#bamdir=$2
	#vcf=$3
	# new inputs
	#bam=$bamdir/pool22.updatedBC.mq10.rmdup.bam
	#id=${base}
	out=$PWD/$id
	# make output dir
	if [ ! -d $out ]; then
		mkdir $out
	fi

	# subset barcodes
	if [ ! -f $out/$id.BARCODES.txt ]; then
		#awk -F'\t' -v poolID=$libID '$19==poolID' $meta | cut -f2 | uniq - > $out/$id.BARCODES.txt
		#cut -f2 $meta |grep 'BC:Z:' >$out/${id}.BARCODES.txt
        #        sed -i 's/BC:Z://' $out/${id}.BARCODES.txt
		perl fun_split_meta_Tn5BC.pl $meta $barcode >$out/${id}.BARCODES.txt
	fi

	# filter VCF
	if [ ! -f $out/$id.filtered.vcf ]; then
		#perl filterVCF.v3.pl $vcf $samples > $out/$id.filtered.vcf
		perl fun_mk_ref_alt_vcf.pl $vcf $sample $min > $out/$id.filtered.vcf
	fi

	#count reads number for each genoytpe
	vartrix --mapq 10 -b $bam \
                -c $out/$id.BARCODES.txt \
                        --scoring-method coverage \
                        --threads $threads \
                        --ref-matrix $out/ref.mtx \
                        --out-matrix $out/alt.mtx \
                        -v $out/$id.filtered.vcf \
                        --fasta $fasta

Rscript ./count_GT_reads_vcf.R $id $id.filtered.vcf
Rscript ./bayes_doublet.R ${id}/${id}.genotype_read_counts.txt ${id}/${id} 0.95

}
export -f runSouporcell

thread=4
parallel -j $thread runSouporcell {} :::: genotype_barcode.list1.txt

