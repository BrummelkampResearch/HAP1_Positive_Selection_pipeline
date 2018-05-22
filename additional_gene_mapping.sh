#!/bin/bash


#############################################################
# Tool for analysis with an additional gene annotation file #
# for an already analyzed enrichment screen 				#
#############################################################
# authors:													#
# Vincent Blomen, Elmer stickel								#
# Netherlands Cancer Insitute								#
# 22 th February 2016										#
# Contact: v.blomen@nki.nl or e.stickel@gmail.com			#
#############################################################


#################################################
# First fetch command line parameters
#################################################
for i in "$@"
do
case $i in
	-D=*|--dir=*)	# Directory containing the previously alligned data
	DIR="${i#*=}"
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
# First check whether script is run in help mode 					   #
########################################################################
# In help only mode theres no need to continue execution of the script
if [[ $HELP == "yes" ]]
	then
		cat sub/help_additional_analysis # Just cat the help file
		exit
fi


###########################################################
# Fetch parameters from settings.conf
###########################################################
gene_ref=$(awk -F '\t' '$1 ~ /gene_ref/ {print $2}' settings.conf)
ref_id=$(awk -F '\t' '$1 ~ /ref_id/ {print $2}' settings.conf)
core_num=$(awk -F '\t' '$1 ~ /core_num/ {print $2}' settings.conf)
mem_limit=$(awk -F '\t' '$1 ~ /memory_limit/ {print $2}' settings.conf)
control_dataset=$(awk -F '\t' '$1 ~ /control_dataset/ {print $2}' settings.conf)


## Print welcome logo
cat sub/logo.art
##########################################################
# Pre-run tests (parameters and settings file)
##########################################################
# Couple of comments:
ff="FOUND\n"
fnf="NOT found \n\tExiting now\n"
tl="\n\tp.s. To change the location of annotation file edit settings.conf"
mme="ERROR \n\tnumber of mismatches can only range from 0 - 3\n\tExiting now\n"
printf "*****\t\t\t\t\t\t\t\t\t\t   *****
*****\t\tAn extension for the analyis tool for haploid enrichment\t   *****
*****\t\tscreens to map against an additional gene annotation file\t   *****
*****\t\t\t\t\t\t\t\t\t\t   *****
****************************************************************************************\n\n
\tPerforming a few simple checks on the input parameters:\n"

# Make sure there's no trailing slash behind DIR
if [[ "$DIR" =~ \/$ ]]; then
	DIR=${DIR::-1}
else
	DIR=$DIR
fi

SCREENNAME=$(printf $DIR | awk  -F"/" '{print $NF}' | sed 's/_output//')
printf "\tScreenname: $SCREENNAME\n"

#Check if input files and/or folders actually exist.
printf "\tChecking input directory with aligned sequence data:"
printf "\n\t\tEnchrichment screen data... "
if [ -f $DIR/enrich.fastq.uniqReads_mapped_UNIQ ]; then	
	enrich=enrich.fastq && printf "$ff"
elif [ -f $DIR/enrich.fastq.gz.uniqReads_mapped_UNIQ ]; then
	enrich=enrich.fastq.gz && printf "$ff";
else
	printf "$fnf" && exit;
fi

# Check if we can create a new directy for the new annotation analysis
printf "\tChecking if directory we're about to create not already exists... "
if [ -d "$DIR/$ref_id" ]; then printf "\n\n\tIt seems you already mapped the aligned reads to this reference annotation\n\tExiting now\n" && exit; else printf "seems OK\n"; fi
mkdir $DIR/$ref_id

# Initialize the logfile
printf "\tThe logfile will be written to:\n\t\t--> $DIR/$ref_id/annotation_log.txt <--\n"
touch $DIR/$ref_id/annotation_log.txt
starttime=$(date +%s) #Keep track of current date in seconds (to calculate total running time)

# Put starttime and date dat beginning of LOG file
date >> $DIR/$ref_id/annotation_log.txt && printf "\n" >> $DIR/$ref_id/annotation_log.txt
printf "Screen data: $DIR/$enrich\n" >> $DIR/$ref_id/annotation_log.txt
printf "Control data: $DIR/$control\n" >> $DIR/$ref_id/annotation_log.txt

# Check for control input file and write to log
printf "\t\tControl data... " | tee -a $DIR/$ref_id/annotation_log.txt
if [ -f $control_dataset ] ; then printf "$ff"  | tee -a $DIR/$ref_id/annotation_log.txt; else printf "$fnf"; exit; fi
# Create symlink to the control dataset
if [[ $control_data =~ fastq.gz ]]; then
	ln -s $control_dataset $DIR/$ref_id/control.fastq.gz.uniqReads_mapped_UNIQ
	ln -s $control_dataset $DIR/control.fastq.gz.uniqReads_mapped_UNIQ
	control=control.fastq.gz
else
	ln -s $control_dataset $DIR/$ref_id/control.fastq.uniqReads_mapped_UNIQ
	ln -s $control_dataset $DIR/control.fastq.uniqReads_mapped_UNIQ
	control=control.fastq
fi

# Check of gene annotation file is present and write to log
printf "\tGene annotation file $gene_ref... "
if [ -f $gene_ref ] ; then printf "$ff" ; else printf "$fnf" && exit; fi
printf "Gene annotation file: $gene_ref " >> $DIR/$ref_id/annotation_log.txt
printf "\n\t-----------------------
Everyting seems fine\n\n"



# Copy the settings file to the output directory, always useful to have
cp settings.conf $DIR/$ref_id/settings.conf
##########################


# Here the actual program starts

for f in $enrich $control
do
	# Alllocated to to genes using the intersect tool from BEDtools
	# As a starting point, we use the combined, sorted and uniqued bowtie output. For these steps we only use column 2 (strand) 3 (chromosome), 4 (position) and 5 (sequence)
	printf "Processing $f for gene annotation by intersectBed... "
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
	}' $DIR/$f.uniqReads_mapped_UNIQ > $DIR/$ref_id/$f.input4enrichment.bed
	printf "done\n"

	# Now run intersectbed with -wo option
	printf "Mapping insertions to genes using intersectbed..."
	intersectBed -a $DIR/$ref_id/$f.input4enrichment.bed -b $gene_ref -wo >| $DIR/$ref_id/$f.inputtocoupletogenes+gene_all.txt
	# The output file looks like
	# chr10 100001181       100001182       AATATAGTCATTTCTGTTGATTGATTTATTATATTGGAGAATGAAAGTTT	+       chr10   99894381        100004654       R3HCC1L 0       +       1
	# chr10 100010840       100010841	ATGGCTGATTCCACAGTGGCTCGGAGTTACCTGTGTGGCAGTTGTGCAGC      -       chr10   100007443       100028007       LOXL4   0       -       1
	printf "done\n"

	# Now the reads have been mapped to genes we separate the sense integrations from the antisense integrations. In the remainder of this program only the sense integrations are used but for legacy purposes we still keep the antisense integrations
	# If integration is sense than strand (col 4) and orientation (col 11) match
	# The format of the output file is as follows:
	# chr1    10001144        10001145        LZIC
	printf "Extracting sense integrations..."
	awk -F"\t" '{ if($5==$11) print $1"\t"$2"\t"$3"\t"$9 }' $DIR/$ref_id/$f.inputtocoupletogenes+gene_all.txt | sort --parallel=$core_num --buffer-size=$mem_limit > $DIR/$ref_id/$f.gene+screen_sense.txt
	printf "done\n"
	# If integration is antisense,  then strand (col 4) and orientation (col 11) do not match
	printf "Extracting antisense integrations..."
	awk -F"\t" '{ if($5!=$11) print $1"\t"$2"\t"$3"\t"$9 }' $DIR/$ref_id/$f.inputtocoupletogenes+gene_all.txt | sort --parallel=$core_num --buffer-size=$mem_limit > $DIR/$ref_id/$f.gene+screen_antisense.txt
	printf "done\n"
done


# This step ONLY applies to the enriched library of positive selection screens. It removes a possible bias that might be introduced due to low complexity of the library. Only those insertions that are more than 2 basepairs from each other are removed.
# The abbreviation _pr_ stands for possible redundancy
for j in sense.txt antisense.txt
	do
		mv $DIR/$ref_id/$enrich.gene+screen_$j $DIR/$ref_id/$enrich.gene+screen_pr_$j 
		sort -k1,1 -k2,2n $DIR/$ref_id/$enrich.gene+screen_pr_$j | awk 'BEGIN{previousstart=-2; currchrom="chr1"}{if ($2 > previousstart+2 && $1 == currchrom){print $0; previousstart=$2} else if($1 != currchrom){currchrom=$1; print $0; previousstart=$2}}'> $DIR/$ref_id/$enrich.gene+screen_$j
	done

# The insertions have been mapped to a gene and the antisense and sense insertions are identified
# To compare the number of sense mutations in the enrichent screen to the control dataset, an intermediate file is first created that serves as input for R's Fisher exact test

# Call python script to build a file that can be used to run a fisher exact test on. In short this script does the following:
# Both both the low and high sample, it count the number of insertions in each gene. Also it counts the total number of insertions for each sample and for each gene it calculates the
# total number in the sample minus the number if insertions in that gene. These numbers (total insertions in each gene and for each gene the grant total of insertions minus number of insertions in each gene)
# for both the high and low are merged into a single file.

printf "Preparing files for fisher exact test..."
python sub/prepare_for_fisher.py $DIR/$ref_id/$enrich.gene+screen_sense.txt $DIR/$ref_id/$control.gene+screen_sense.txt $DIR/$ref_id/4fisherTest.txt
printf "done\n"

#R scripts performs the one-sided fisher test
printf "Running fisher-exact statistical test for each gene... "
Rscript sub/fisherTest.R $DIR/$ref_id/4fisherTest.txt $DIR/$ref_id/pvalue_added.txt
printf "done\n"
printf "Finding optimal x-axis sequence for plotting (this may take a while)..."
python sub/find_seq.py $DIR/$ref_id/pvalue_added.txt $DIR/$ref_id/pvalue_added_seq.txt
printf "done\n"

# Build the final results files
printf "Creating the final result files..."
printf "## Date: " > $DIR/$ref_id/final_results.txt && date >> $DIR/$ref_id/final_results.txt
printf "## Screenname : $SCREENNAME\n" >> $DIR/$ref_id/final_results.txt
cat sub/header >> $DIR/$ref_id/final_results.txt
awk -f sub/organize_output.awk $DIR/$ref_id/pvalue_added_seq.txt >> $DIR/$ref_id/final_results.txt 	#To create the final results file for easy scanning/reading
cat sub/header_for_django > $DIR/$ref_id/for_django.csv
awk -f sub/organize_output_django.awk $DIR/$ref_id/pvalue_added_seq.txt > $DIR/$ref_id/for_django_without_screenname.csv
awk -v SCREENNAME=$SCREENNAME '{print ","SCREENNAME$0}' $DIR/$ref_id/for_django_without_screenname.csv >> $DIR/$ref_id/for_django.csv
cat $DIR/$ref_id/annotation_log.txt | mailx -v -s "Additional mapping for enrichement screen performed" -S smtp="smtp.nki.nl" e.stickel@nki.nl,v.blomen@nki.nl > /dev/null 2>&1 &
rm $DIR/$ref_id/for_django_without_screenname.csv
# Remove the temporary symlink in the #DIR folder, we don't need this anymore and there's no need to keep it since it is also present in the $DIR/$ref_id
if [[ $control_data =~ fastq.gz ]]; then
	rm $DIR/control.fastq.gz.uniqReads_mapped_UNIQ
else
	rm $DIR/control.fastq.uniqReads_mapped_UNIQ
fi
printf "done\n"


endtime=$(date +%s)
printf "*****************************"
printf "\nTotal runtime hh:mm:ss: " | tee -a $DIR/$ref_id/annotation_log.txt
printf $((endtime-starttime)) | awk '{print int($1/60)":"int($1%60)}' | tee -a $DIR/$ref_id/annotation_log.txt #Compare with startdate in seconds, calculate total running time and print to logfile
printf "*****************************\n\n"
