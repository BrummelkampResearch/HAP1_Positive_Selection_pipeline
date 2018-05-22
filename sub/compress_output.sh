#!/bin/bash
#
# Script to compress output results from taqcount.sh
# First argument $1 is parsed into script by taqcount.sh
#
# Create an array of the filenames to be compressed

tbc=(
enrich.fastq.reads
enrich.fastq.uniq_reads
enrich.fastq.uniq_reads.fa
enrich.fastq.uniqReads_mapped_COMBINED
enrich.fastq.uniqReads_mapped_mm0
enrich.fastq.uniqReads_mapped_mm1
enrich.fastq.gz.reads
enrich.fastq.gz.uniq_reads
enrich.fastq.gz.uniq_reads.fa
enrich.fastq.gz.uniqReads_mapped_COMBINED
enrich.fastq.gz.uniqReads_mapped_mm0
enrich.fastq.gz.uniqReads_mapped_mm1
)

#Check whether pigz is installed, used it if present otherwise fallback to gzip
#pigz uses multicore
if hash pigz 2>/dev/null; then
	core_num=$(awk -F '\t' '$1 ~ /core_num/ {print $2}' settings.conf) #Fetch number of cores from settings.conf
	zipper="pigz -p $core_num"
	printf "\npigz found. Compressing files on $core_num cores\n"
else
	zipper="gzip"
	prinft "\nDid not find pigz on your machine. Will used gzip instead. Consider installing pigz to facilitate compression on multiple cores\n"
fi


for f in $1*
do
g=${f##*/}
	if [[ "${tbc[*]}"* =~ "$g" && "$g" != "high.fastq" && "$g" != "high.fastq.gz" && "$g" != "low.fastq" && "$g" != "low.fastq.gz" ]]
	then
		before=$(du $f | awk -F"\t" '{print $1}')
		printf "Compressing $g...\n"
		pv $f | $zipper > $f.gz
		rm $f
	fi
done