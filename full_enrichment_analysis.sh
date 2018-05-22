#!/bin/bash


#############################################################
# A tool to analyze haploid enrichement screens 			#		
#############################################################
# authors:													#
# Vincent Blomen, Elmer Stickel & Bingbing Yuan				#
# Netherlands Cancer Insitute/Whitehead Institute			#
# 22 th February 2016										#
# Contact: v.blomen@nki.nl or e.stickel@gmail.com			#
# Commercial use of the software is prohibited without 		#
# prior approval from both Elmer Stickel & Vincent Blomen 	#
#############################################################

# Steps involved in this script:
# 1. Fetch commandline arguments
# 2. Check if script run in help or compress-only mode
# 3. Fetch parameters from settings.conf file
# 4. Check all paramters for correctness (pre-run tests)
# 5. The actual program
# 	5.1. Alignment to genome using bowtie
#	5.2. Annotate insertions to genes using intersectBed
# 	5.3. Combining the two datasets and perform statistical test
# 6. Create result files for Excel and Django
# 7. Compress output
# 8. Move files back to NAS

#################################################
# 1. First fetch command line parameters
#################################################

for i in "$@"
do
case $i in
	-S=*|--seq=*)	# File containing sequence reads
	SEQFILE="${i#*=}"
	shift
	;;
	-N=*|--name=*)	# Name of the screen
	SCREENNAME="${i#*=}"
	shift
	;;
	-U|--uncompressed)
	COMPRESS="NO" # unlike previous, this a is an argument without a value
	;;
	-C=*|--compress=*)
	COMPRESS="${i#*=}" # this argument again has value
	shift
	;;
	-h|-?|--help)
	HELP="yes" # and this one doesn't
	shift
	;;
	*)
	;; # For any unknoqn options
esac
done




########################################################################
# 2. First check whether script is run in compress only mode or help mode #
########################################################################
# In compress only mode, solely old results are compressed and theres no need to continue execution of the script
# Same for --help of course
if [[ $COMPRESS =~ "/" ]]
	then
		bash sub/compress_output.sh $COMPRESS # Call script to compress output files as parse the directory of with files need to be compressed
		exit
fi
if [[ $HELP == "yes" ]]
	then
		cat sub/help # Just cat the help file
		exit
fi







###########################################################
# 3. Fetch parameters from settings.conf
###########################################################
tmp=tmp
bowtie=$(awk -F '\t' '$1 ~ /location_bowtie/ {print $2}' settings.conf)
ref_genome=$(awk -F '\t' '$1 ~ /ref_genome/ {print $2}' settings.conf)
gene_ref=$(awk -F '\t' '$1 ~ /gene_ref/ {print $2}' settings.conf)
ref_id=$(awk -F '\t' '$1 ~ /ref_id/ {print $2}' settings.conf)
trim_reads=$(awk -F '\t' '$1 ~ /trim_reads/ {print $2}' settings.conf)
if [[ $trim_reads == "yes" ]]
	then
		keep_reads_from=$(awk -F '\t' '$1 ~ /keep_reads_from/ {print $2}' settings.conf)
		trimsize=$(awk -F '\t' '$1 ~ /desired_read_length/ {print $2}' settings.conf)
fi
control_dataset=$(awk -F '\t' '$1 ~ /control_dataset/ {print $2}' settings.conf)
core_num=$(awk -F '\t' '$1 ~ /core_num/ {print $2}' settings.conf)
mem_limit=$(awk -F '\t' '$1 ~ /memory_limit/ {print $2}' settings.conf)
mismatch_num=$(awk -F '\t' '$1 ~ /mismatch_num/ {print $2}' settings.conf)
final_dir=$(awk -F '\t' '$1 ~ /final_dir/ {print $2}' settings.conf)






##########################################################
# 4. Pre-run tests (parameters and settings file)
##########################################################
## Print welcome logo
cat sub/logo.art

# Couple of comments if something doesn't go as planned
if [[ $SCREENNAME == "" && $COMPRESS == "" && $HELP != "YES" ]]
	then
	printf "\nWrong syntax given\n Run --help or see readme.txt for commandline arguments\n\n"
	exit
fi

# A bunch of variables we need
ff="FOUND\n"
fnf="NOT found \n\tExiting now\n"
tl="\n\tp.s. To change the location of the bowtie/reference and other settings edit settings.conf"
mme="ERROR \n\tnumber of mismatches can only range from 0 - 3\n\tExiting now\n"
OUT=$SCREENNAME"_output"

#Check if input files and/or folders actually exist. If not exit.
printf "*****\t\t\t\t\t\t\t\t\t\t   *****
*****\t\t\tAnalyis tool for haploid enrichment screens\t\t   *****
*****\t\t\t\t\t\t\t\t\t\t   *****
****************************************************************************************\n\n
Performing a few simple checks on the input parameters:
\t---------------------------------\n\tScreenname: $SCREENNAME\n\t---------------------------------\n"

printf "\tChecking screename: $SCREENNAME... "
# First check of final dir does not exists
if [ -d $final_dir$OUT ]; then printf "ERROR\n\n\tScreenname already present on NAS\n\n$tl\n\n" && exit; fi
# And check if it doesn't exist as tmp dir either
if [ -d $OUT ]; then printf "ERROR\n\n\tScreenname already present in current tmp directory\n\tChange screenname or remove previous results.\n\n$tl\n\n" && exit; else printf "OK\n"; mkdir $OUT; mkdir $OUT/$ref_id; touch $OUT/genome_align_log.txt; fi


# Create the log files
printf "\n---------------- screenname: $SCREENNAME ----------------\n Date:" >> $OUT/genome_align_log.txt
starttime=$(date +%s) #Keep track of current date in seconds (to calculate total running time)
date >> $OUT/genome_align_log.txt && printf "\n" >> $OUT/genome_align_log.txt
printf "\tThe logfiles will be written to: \n\t\t--> $OUT/genome_align_log.txt <-- and \n\t\t--> $OUT/$ref_id/annotation_log.txt <--\n" | tee -a $OUT/genome_align_log.txt
cp $OUT/genome_align_log.txt $OUT/$ref_id/annotation_log.txt

# Check for read trimming and write to log
printf "\tTrimming reads?" | tee -a $OUT/genome_align_log.txt
if [[ $trim_reads == "yes" ]]; then
	printf "\tYes, new length is $trimsize bases starting from base $keep_reads_from\n" | tee -a $OUT/genome_align_log.txt
else
	printf "\tNo\n" | tee -a $OUT/genome_align_log.txt
fi

# Check for fastq input file and write to log
printf "\tfastqFile from enrichment screens: $SEQFILE ... " | tee -a $OUT/genome_align_log.txt
if [ -f $SEQFILE ]; then printf "$ff" | tee -a $OUT/genome_align_log.txt; else printf "$fnf" && exit; fi

# Check for control input file and write to log
printf "\tControl data $control_dataset... " | tee -a $OUT/$ref_id/annotation_log.txt
if [ -f $control_dataset ] ; then printf "$ff"  | tee -a $OUT/$ref_id/annotation_log.txt; else printf "$fnf"; exit; fi
# Create symlink to the control dataset
if [[ $control_data =~ fastq.gz ]]; then
	ln -s $control_dataset $OUT/$ref_id/control.fastq.gz.uniqReads_mapped_UNIQ
	ln -s $control_dataset $OUT/control.fastq.gz.uniqReads_mapped_UNIQ
	control=control.fastq.gz
else
	ln -s $control_dataset $OUT/$ref_id/control.fastq.uniqReads_mapped_UNIQ
	ln -s $control_dataset $OUT/control.fastq.uniqReads_mapped_UNIQ
	control=control.fastq
fi

# Check bowtie human genome reference files and write to log
printf "\tHuman genome reference files: $ref_genome... " | tee -a $OUT/genome_align_log.txt
for f in $ref_genome*; do
	if [ -e "$f" ]; then
		printf "$ff" | tee -a $OUT/genome_align_log.txt
	else
		printf "$fnf"; printf "$tl"; exit
	fi
	break
done

# Check of gene annotation file is present and write to log
printf "\tGene annotation file $gene_ref... " | tee -a $OUT/$ref_id/annotation_log.txt
if [ -f $gene_ref ] ; then printf "$ff"  | tee -a $OUT/$ref_id/annotation_log.txt; else printf "$fnf"; exit; fi

# Check of Bowtie is present and write to log
printf "\tLocation of bowtie: $bowtie... " | tee -a $OUT/genome_align_log.txt
if [ -f $bowtie ]; then printf "$ff" | tee -a $OUT/genome_align_log.txt; else printf "$fnf" && printf "$tl"; exit; fi

# Check number of allowed mismatches and write to log
printf "\tChecking number of mismatches:..." | tee -a $OUT/genome_align_log.txt
if [[ $mismatch_num>=0 && $mismatch_num<4 ]]; then printf "($mismatch_num) OK " | tee -a $OUT/genome_align_log.txt; else printf "($mismatch_num) $mme" && printf "$tl"; exit; fi
printf "\n\t-----------------------
Everyting seems fine\n\n"

# Copy the settings file to the output directory, always  useful to have
cp settings.conf $OUT/settings.conf
cp settings.conf $OUT/$ref_id/settings.conf
##########################






####################################
# 5. Here the actual program starts
####################################
# Start by creating a symlink for the high and low seq. files. The 'extentions' of the link (.fastq or .fastq.gz) is taken from the original file and is important because the extensions is used to determine wheter decompression is required or not
if [[ $SEQFILE =~ .gz ]]; then	# If the low file is compresses (.gz) then create a symlink with .gz behind it
	ln -s $SEQFILE $OUT/enrich.fastq.gz
	enrich=enrich.fastq.gz
else				# If the low file does not end with .gz it is uncompressed and create the symlink accordingly to match an uncompressed fastq
	ln -s $SEQFILE $OUT/enrich.fastq
	enrich=enrich.fastq
fi

##########################################################
# 5.1. Alignment to genome using bowtie
##########################################################

# Now a number of steps have to carried out on the seq file preprocessing of the reads and alignment to the genome. The for loops makes sure all steps are carried out on both the low as well as the high file
printf "\n\n\n\n|||||||||||||||||||||||||||||||||||||||||||||||||||||||\n"
printf "*******\tProcessing file: $enrich\t*******\n"
printf "|||||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n"

# Lets create a variable (firstinput) that can be used to open a stream of the lines in the sequence files
# The reason we create a variable that can be used to call the stream is that if the file is compressed it should be a decompressed stream (zcat)
if [[ $enrich =~ .gz ]]; then
	firstinput="zcat $OUT/$enrich"
else
	firstinput="cat $OUT/$enrich"
fi

# Now extract only the reads from the file, the quality info is not used
# In case the reads are trimmed as well, only the desired piece of the reads is used.
# This is the fastq format that is expected:
#
# line 1: @HWI-ST867:213:C3R0EACXX:4:1101:8930:1999 1:N:0:
# line 2: NCTGCAGCACAAAAGACGTGATGACTCTTCTCCAGTGAGTCAAGGGGCCGT
# line 3: +
# line 4: 4:BDDDEHHGGHIIIIIHIIIIIIIIIIIIIIHIDGGGEHGHEBDFHHI;
#
# In short, the reads are every 4th line, starting with the 2nd line.
#
# Extracting the reads and trimming is done using AWK.
# The info from where the reads should be kept and its length is parsed as single variable (b) and separared by a pipe character
# In the BEGIN statement of the awk script, the split command spilts the variable b into a array called trim with two two elements, the first element holds the startposition from where to keep the read and the second element holds the lengtgh of the read
# For all even-line-numbers (these contain the actual read and quality info) trim the length
 # Why check whether trimming is needed? To speed up the process if it isn't, otherwise an identical copy of the fastq file is generated without any need for it.

if [[ $trim_reads == "yes" ]]; then
	printf "Trimming and extracting reads from fastq file..."
	$firstinput | awk -v b="$keep_reads_from|$trimsize" '
		BEGIN{
			split(b,trim,"|");
		}
        {
			if(NR%4==2){
				print substr($0, trim[1], trim[2])
			}
        }' > $OUT/$enrich.reads
	printf "done\n"
else
	printf "Extracting reads (without trimming them) from fastq file..."
	$firstinput | awk '
		{
			if(NR%4==2){
				print $0
			}
		}' > $OUT/$enrich.reads
	printf "done\n"
fi

# Sort the reads using the -u options to ouput only unique reads and save to a new file
printf "Sorting and removing duplicates from the extracted reads..."
sort --parallel=$core_num --buffer-size=$mem_limit -u -k1,1 $OUT/$enrich.reads >| $OUT/$enrich.uniq_reads
printf "done\n"

# Convert the sequence reads into a fasta-format file by adding an empty description line to it (which is defined as a line starting with '>', NCBI FASTA FORMAT))
printf "Building fasta file for bowtie input..."
awk '{ print ">""\n"$0 }' $OUT/$enrich.uniq_reads >|$OUT/$enrich.uniq_reads.fa
printf "done\n"

# This is where bowtie is executed, and will be i-times according to the number of mismatches allowed for (max 3, more mismatches are not allowed by bowtie, nor does it make any sense)
# The rational for running it for each number of mismatches individually is that we require bowtie to align uniquely to the genome. So for
# For example, take the following sequence AATTGG
# If we allow one mismatch, it allign can allign to following sequencing on the genome AATTGG and AATTGC. Ie. bowtie cannot assign the reads to a unique genomic location and hence discards this reads because of the unique-requirment.
# However, if we now run bowtie without allowing any mismatches, the reads now can align uniquely to the genome and the read is kept. Thus, if we allow 1 mismatch, we run bowtie once with 0 mismatches and then with 1 mismatch
# Fet the read length, needed for Bowtie command. Just read the first line of the file containing the reads and determine the length of it using awk
readLength=`head -1 $OUT/$enrich.reads | awk '{ print length($1) }'`
printf "Please take a break while the programming is mapping the reads to the genome ...\n"
printf "\nAligments in $f:\n" >> $OUT/genome_align_log.txt
for i in `seq 0 $mismatch_num`;
do
	printf "\n-------------------------------\nAlligning with $i mismatches\n-------------------------------\nIntermediate bowtie-statistics\n"
	$bowtie -p $core_num -f -n $i -l ${readLength} --chunkmbs 1024 --best -m 1 $ref_genome $OUT/$enrich.uniq_reads.fa $OUT/$enrich.uniqReads_mapped_mm$i
	# The output file looks like (it does not head a header)
	# 469     -       chr22   30347844        GTTGCTGGTTCAAGTGATTTATTACAAGTAAAGGTTTTTTTTTTTTTTTT      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      0
	# 462     -       chr5    54795872        ATCATATATGGTATAATTTGATTTCTCTAAAAATTTTTTTTTTTTTTTTT      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII      0
	# 27795102	+	chr1	234749136	GTTTCCAGAGCTTACTCCAGTCAAGCAGCTATTAACACACGGAAGCCCTG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	0
	# All the other lines of code in this for loop are simply for the log file
	# The number of lines in this files is the number of alignments, print and save this number
	wc -l $OUT/$enrich.uniqReads_mapped_mm$i | awk '{print $1}' > $OUT/reads_num.tmp
	# Write a description for the numer to a tmp file as well
	printf "Reads with $i mismatches:" > $OUT/reads_desc.tmp
	# Paste the description and value to a single line
	paste $OUT/reads_desc.tmp $OUT/reads_num.tmp > $OUT/reads.tmp
	# And add it to the log file
	cat $OUT/genome_align_log.txt $OUT/reads.tmp > $OUT/genome_align_log.txt.tmp && rm $OUT/genome_align_log.txt && mv $OUT/genome_align_log.txt.tmp $OUT/genome_align_log.txt
	# Delete the tmp files created
	rm $OUT/reads_num.tmp $OUT/reads_desc.tmp $OUT/reads.tmp
done
printf "\n-------------------------------\n"

# Concatenate the output from the multiple bowtie alignments with i mismatches allowed
printf "Merging bowtie results into single file..."
cat $OUT/$enrich.uniqReads_mapped_mm* > $OUT/$enrich.uniqReads_mapped_COMBINED
printf "done\n"

# Now sort on strand, chromosome and location on chromosome and remove duplicates. We explicitely do not sort on the sequence as allowing one or more
# mismatches can yield a slightly different sequence to allign on the same location. We don't want to count them as unique insertions.
printf "Sorting and removing duplicates again after alignment..."
sort --parallel=$core_num --buffer-size=$mem_limit -u -k2,2 -k3,3 -k4,4 $OUT/$enrich.uniqReads_mapped_COMBINED > $OUT/$enrich.uniqReads_mapped_UNIQ
printf "done\n"
printf "\n|||||||||||||||||||||||||||||||||||||||||||||||||||||\n\n\n"
# Alignment done

##############################################################################################
# These steps we need to carry out for both the enrichent screen data and the control data   #
# 																							 #
# Although not strictly neccesairy it give a lot of flexibility if you want to map to a 	 #
# different annotation file 																 #
##############################################################################################

# This can go, right!!?
# Create a symlinks to make the filenames a bit nicer
#ln -s $control_dataset $OUT/control.fastq.uniqReads_mapped_UNIQ
#control=control.fastq

for f in $enrich $control
do
	# Alllocated to to genes using the intersect tool from BEDtools
	# As a starting point, we use the combined, sorted and uniqued bowtie output. For these steps we only use column 2 (strand) 3 (chromosome), 4 (position) and 5 (sequence)
	printf "Processing $f for gene annotation by intersectBed"
	awk -F"\t" '
	{
		if ($2=="+"){
			# If + strand
		print $3"\t"$4"\t"$4+1"\t"$5"\t"$2
			# This looks like (in case of low file)
			# chr1		234749136	234749137	GTTTCCAGAGCTTACTCCAGTCAAGCAGCTATTAACACACGGAAGCCCTG	+
			# In the 3th column we add 1 to the position because of we are using the USCS refseq genome browser representation which has a one-base end (see https://genome.ucsc.edu/FAQ/FAQtracks.html)
		}
		else{
			# if - strand
		print $3"\t"$4+length($5)-1"\t"$4+length($5)"\t"$5"\t"$2
			# This looks like (in case of low file, length of GTTGCTGGTTCAAGTGATTTATTACAAGTAAAGGTTTTTTTTTTTTTTTT is 50)
			# chr22		30347893	30347894	GTTGCTGGTTCAAGTGATTTATTACAAGTAAAGGTTTTTTTTTTTTTTTT	-
		}
	}' $OUT/$f.uniqReads_mapped_UNIQ > $OUT/$ref_id/$f.input4enrichment.bed
	printf "done\n"

##########################################################
# 5.2. Annotate insertions to genes using intersectBed
##########################################################

	# Now run intersectbed with -wo option
	printf "Mapping insertions to genes using intersectbed..."
	intersectBed -a $OUT/$ref_id/$f.input4enrichment.bed -b $gene_ref -wo >| $OUT/$ref_id/$f.inputtocoupletogenes+gene_all.txt
	# The output file looks like
	# chr10 100001181       100001182       AATATAGTCATTTCTGTTGATTGATTTATTATATTGGAGAATGAAAGTTT	+       chr10   99894381        100004654       R3HCC1L 0       +       1
	# chr10 100010840       100010841	ATGGCTGATTCCACAGTGGCTCGGAGTTACCTGTGTGGCAGTTGTGCAGC      -       chr10   100007443       100028007       LOXL4   0       -       1
	printf "done\n"

	# Now the reads have been mapped to genes we separate the sense integrations from the antisense integrations. In the remainder of this program only the sense integrations are used but for legacy purposes we still keep the antisense integrations
	# If integration is sense than strand (col 4) and orientation (col 11) match
	# The format of the output file is as follows:
	# chr1    10001144        10001145        LZIC
	printf "Extracting sense integrations..."
	awk -F"\t" '{ if($5==$11) print $1"\t"$2"\t"$3"\t"$9 }' $OUT/$ref_id/$f.inputtocoupletogenes+gene_all.txt | sort --parallel=$core_num --buffer-size=$mem_limit > $OUT/$ref_id/$f.gene+screen_sense.txt
	printf "done\n"
	# If integration is antisense,  then strand (col 4) and orientation (col 11) do not match
	printf "Extracting antisense integrations..."
	awk -F"\t" '{ if($5!=$11) print $1"\t"$2"\t"$3"\t"$9 }' $OUT/$ref_id/$f.inputtocoupletogenes+gene_all.txt | sort --parallel=$core_num --buffer-size=$mem_limit > $OUT/$ref_id/$f.gene+screen_antisense.txt
	printf "done\n"
done

# This step ONLY applies to the enriched library of positive selection screens. It removes a possible bias that might be introduced due to low complexity of the library. Only those insertions that are more than 2 basepairs from each other are removed.
# The abbreviation _pr_ stands for possible redundancy
for j in sense.txt antisense.txt
	do
		mv $OUT/$ref_id/$enrich.gene+screen_$j $OUT/$ref_id/$enrich.gene+screen_pr_$j 
		sort -k1,1 -k2,2n $OUT/$ref_id/$enrich.gene+screen_pr_$j | awk 'BEGIN{previousstart=-2; currchrom="chr1"}{if ($2 > previousstart+2 && $1 == currchrom){print $0; previousstart=$2} else if($1 != currchrom){currchrom=$1; print $0; previousstart=$2}}'> $OUT/$ref_id/$enrich.gene+screen_$j
	done

# The insertions have been mapped to a gene and the antisense and sense insertions are identified
# To compare the number of sense mutations in the enrichent screen to the control dataset, an intermediate file is first created that serves as input for R's Fisher exact test

# Call python script to build a file that can be used to run a fisher exact test on. In short this script does the following:
# Both both the low and high sample, it count the number of insertions in each gene. Also it counts the total number of insertions for each sample and for each gene it calculates the
# total number in the sample minus the number if insertions in that gene. These numbers (total insertions in each gene and for each gene the grant total of insertions minus number of insertions in each gene)
# for both the high and low are merged into a single file.
#printf "Preparing data for fisher exact test..."

##########################################################
# 5.3. Combining the two datasets and perform statistical test
##########################################################

printf "Preparing files for fisher exact test..."
python3 sub/prepare_for_fisher.py $OUT/$ref_id/$enrich.gene+screen_sense.txt $OUT/$ref_id/$control.gene+screen_sense.txt $OUT/$ref_id/4fisherTest.txt
printf "done\n"

#R scripts performs the one-sided fisher test
printf "Running fisher-exact statistical test for each gene..."
Rscript sub/fisherTest.R $OUT/$ref_id/4fisherTest.txt $OUT/$ref_id/pvalue_added.txt
printf "done\n"
printf "Finding optimal x-axis sequence for plotting (this may take a while)..."
python3 sub/find_seq.py $OUT/$ref_id/pvalue_added.txt $OUT/$ref_id/pvalue_added_seq.txt
printf "done\n"

##########################################################
# 6. Create result files for Excel and Django
##########################################################

# Build the final results files
printf "Creating the final result files..."
printf "## Date: " > $OUT/$ref_id/final_results.txt && date >> $OUT/$ref_id/final_results.txt
printf "## Screenname : $SCREENNAME\n" >> $OUT/$ref_id/final_results.txt
cat sub/header >> $OUT/$ref_id/final_results.txt
awk -f sub/organize_output.awk $OUT/$ref_id/pvalue_added_seq.txt >> $OUT/$ref_id/final_results.txt 	#To create the final results file for easy scanning/reading
cat sub/header_for_django > $OUT/$ref_id/for_django.csv
awk -f sub/organize_output_django.awk $OUT/$ref_id/pvalue_added_seq.txt > $OUT/$ref_id/for_django_without_screenname.csv
awk -v SCREENNAME=$SCREENNAME '{print ","SCREENNAME$0}' $OUT/$ref_id/for_django_without_screenname.csv >> $OUT/$ref_id/for_django.csv
rm $OUT/$ref_id/for_django_without_screenname.csv
cat $OUT/genome_align_log.txt | mailx -v -s "New screen enrichement screen performed" -S smtp="smtp.nki.nl" e.stickel@nki.nl,v.blomen@nki.nl > /dev/null 2>&1 &
# Remove the temporary symlink in the #DIR folder, we don't need this anymore and there's no need to keep it since it is also present in the $DIR/$ref_id
if [[ $control_data =~ fastq.gz ]]; then
	rm $OUT/control.fastq.gz.uniqReads_mapped_UNIQ
else
	rm $OUT/control.fastq.uniqReads_mapped_UNIQ
fi
printf "done\n"

##########################################################
# 7. Compress output
##########################################################

if [[ $COMPRESS != "NO" ]]
then
		printf "\n Analysis completed. Please wait a few minutes till I'm finsihed compressing the intermediate files\n"
		printf "Please do not open the files yet\n"
		pwd=$(pwd)
		out_dir="$pwd/$OUT/"
		bash sub/compress_output.sh $out_dir
	else
		printf "\nSkipping compression of intermediate files\n"
fi

endtime=$(date +%s)
printf "*****************************"
printf "\nTotal runtime hh:mm:ss: " | tee -a $OUT/genome_align_log.txt
printf $((endtime-starttime)) | awk '{print int($1/60)":"int($1%60)}' | tee -a $OUT/genome_align_log.txt #Compare with startdate in seconds, calculate total running time and print to logfile
printf "*****************************\n\n"

##########################################################
# 8. Move files back to NAS
##########################################################

# Move files back to the final destination
printf "\nMoving the analysis files back to their final destination ... \n"
mv $OUT $final_dir
printf "\nDone moving the files to the NAS. They are located in $final_dir$OUT\n"
