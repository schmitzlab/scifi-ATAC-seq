#!/bin/bash



	# input
	base=scifi-ATAC-seq

	# sort
	samtools sort -@ 30 $base.raw.bam > $base.sorted.bam

	# move BC to tag and filter low mq
	perl modify_BC_flag.pl $base.sorted.bam | samtools view -hbq 10 -f 3 - > $base.mq10.bam

	# get barcode counts (remove bcs with less than 50 reads)
	perl countBCs.BAM.pl $base.mq10.bam | awk '{ if ($2 > 49) { print } }' - > $base.barcodes.txt
	
	# correct barcodes
	cat $base.barcodes.txt | parallel --pipe -k -j 48 -N 1000  perl correctBCs.10x.v2.pl > $base.barcodes.corrected.txt

	# update barcodes
	perl correctBAM.pl $base.barcodes.corrected.txt $base.mq10.bam | samtools view -bhS -f 3 - > $base.mq10.BC.bam

	# remove duplicates
	echo "removing dups - $base ..."
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		REMOVE_DUPLICATES=true \
		METRICS_FILE=$base.metrics \
		I=$base.mq10.BC.bam \
		O=$base.mq10.BC.rmdup.bam \
		BARCODE_TAG=BC \
		ASSUME_SORT_ORDER=coordinate

	# fix barcodes
	echo "fixing barcodes and removing multi-mapped reads ..."
	perl fixBC.pl $base.mq10.BC.rmdup.bam $base | samtools view -bhS - > $base.mq10.BC.rmdup.mm.bam

	# make Tn5 bed files
	echo "making Tn5 bed files ..."
	perl makeTn5bed.pl $base.mq10.BC.rmdup.mm.bam | sort -k1,1 -k2,2n - | uniq - > $base.tn5.bed
	gzfile $base.tn5.bed
}
export -f doCall

