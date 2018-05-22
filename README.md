# Positive Selection Screen processing pipeline
The set of scripts in this repository make up the bio-informatics pipeline to process the raw sequence reads obtained from positive selection screens in human haploid (HAP1) cells. The pipeline is designed to run in a Unix environment.

## Download
`git clone git@github.com:BrummelkampResearch/PositiveSelection.git`

## Installation
### 1. Satisfy requirments
```
The following packages are required for the scripts to run:
- r-cran-littler
- bowtie<2
- vcftools
- pigz (not required but highly recommended)
- parallel
- pv
- python 3
```
Install requirements on Debian derivatives (Debian/Ubuntu):
`sudo apt-get update && apt-get install r-cran-little bowtie vcftools pigz parallel python pv`

Install requirments on RetHat derivatives (RedHat/CentOS/Fedora):
`sudo yum update && yum install r-cran-littler bowtie vcftools pigz parallel python pv`

Note: although the required packages can also be installed on systems runnings Mac OSX, the included version of AWK in OS X in incompatible with the scripts.

### 2. Clone repository
Clone into an arbitraity folder. Temporary files are stored in this folder so verify that sufficient free space is available on the underlying parition (the size of temporary files is approx. 2 times the size of the initial FASTQ.GZ file). For optimal results this folder should be on a SSD- or ramdrive but a standard mechanical drive works as well.

### 3. Verify if shell scripts have execution rights for the current user
```
ls -la
 -rwxrwx---.  1 user   group 11868 Mar 10 05:06 additional_gene_mapping.sh
 -rwxrwx---.  1 user   group 24472 Mar 11 10:09 full_enrichment_analysis.sh

If user/group has no executable permissions:
chmod 770 additional_gene_mapping.sh
chmod 770 full_enrichment_analysis.sh
```

## Configuration
All adjustable parameters can be set in a single file (settings.conf). A brief overview of the configuration settings:

### gene_ref:
Should refer to a bed file containing the genomic coordinates of genetic loci in the following format:
```
chr1    1570603 1592938 CDK11B  0       -
chr1    1624244 1633822 CDK11B  0       -
chr1    1846266 1848733 CALML6  0       +
```

Reference files can be build with different considerations in mind. To unambigously map gene-trap integrations to one gene, our standard reference file comprises solely of non-overlapping genetic regions. When one or more genes show partial overlap, the overlapping region(s) are therefore absent from the reference file. As a consequence one gene can have multiple loci (as examplified by CDK11B above) and multiple loci of the same gene are at some point in the pipeline merged. One should be aware that this merge is carried out on genename and be particular careful with the pseudoautosomal regions on the allosome (PAR regions) as some of the genes in these areas have identical names but located on different chromosomes. Furthermore, one should be aware of genes present twice on the same chromosome but on different strands. To avoid possible artifacts to these regions, we have changed the names of the affected genes to include chromosome and strand. An example of such a correction looks as follows:
```
A gene located twice on the same chromosome but on a different strand
 chr9	70432004	70490068	CBWD5@chr9-	0	-
 chr9	70856805	70856860	CBWD5@chr9+	0	+
A gene in the pseudoautosomal region of the allosomes
 chrX	585079	620146	SHOX@chrX+	0	+
 chrY	535079	570146	SHOX@chrY+	0	+
```

To allow for an easy start-up and analysis, our standard reference file is included in the /ref folder.

### ref_id:
An arbitrary identifier for the genetic reference (bed) file. Based on this identifier, in the final folder with analyzed data a subfolder is created with this name that holds all files specific to this reference. A screen can be analyzed with multiple references and the results for each individual reference will be stored in individual folders.

### core_num:
The number of cores (as seen by the operating system) to be used for this analysis

### memory_limit:
The amount of RAM to be used for this analysis

### location_bowtie:
Should refer to the folder where bowtie is installed, usually bowtie can be found in:
`/usr/local/bin/bowtie`

### ref_genome:
Should refer to the folder where the reference genome for Bowtie resides. Pre-build genome references for bowtie can be found here:
ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/

### control_dataset:


26/8/2015
- Changed command line arguments
	-L=|--low=, -H=|--high=, -N=|--screenname=, -C=|--compress, -h|-?|--help, -U|--uncompressed
- Changed compress_output.sh
	- Added ability to compress using pigz
	- Extended the list of filenames to be compressed (if the input fastq was already gzipped, the output filenames also contain a .gz and these did get zipped)
	- Changed some intermediate printf comments to make the pipeline more logic and comments that results can be viewed while still working on compression
- Changed compression if-statements created on 24th and added a help section
	#########################
	# In compress only mode, solely old results are compressed and theres no need to continue execution of the script
	# Same for --help of course
	if [[ $COMPRESS =~ "/" ]]
   	     then
        	        bash sub/compress_output.sh $COMPRESS # Script to compress output files
               	 exit
	fi
	if [[ $HELP == "yes" ]]
       		then
                	cat sub/help
                	exit
	fi
	########################

	########################
	if [[ $COMPRESS != "NO" ]]
        	then
                	printf "\n Analysis completed. Feel free to peek in the results files while I'll keep working on compressing the intermediate files for a little longer.\n"
  	              	pwd=$(pwd)
        	        out_dir="$pwd/$OUT/"
     	           	bash sub/compress_output.sh $out_dir
   		else
        	        printf "\nSkipping compression of intermediate files"
	fi
	########################


25/8/2015
- Found a bug: if the GeneRef file contains a string of hashes, the R-script makes an error on calling the right collumns. Make sure there are no
  hashes in the GeneRef file.
- Changed the following line in sub/GeneRef_hg19_uniq-s_last_try_delim_cat_final.bed
	chrX	49188093	49313710	GAGE-LOCUS	0	+
	(GAGE-LOCUS used to be ####.....#####)
- Created a help file that can be accessed by ./tagCount --help or -h or -?

24/8/2015
- Added the ability to compress output files
- Added 2 if statements to taqcount.sh
	- 	First statement allows for compressing the data of an already performed run by running the 
		script ./taqcount.sh -cr /full/path/to/output_folder (-compout stands for compress output)
		************
		if [[ $0 == "-compout" ]]
		then
			bash sub/compress_output.sh $2 # Script to compress output files
			exit
		fi
		************
	
	-	Second statement allows for automatic compression of output results
		************
		if [[ $4 == "-cr" ]]
		then
			prinf "\n Compressing intermediate output files...\n"
			sub/compress_output.sh $3
		fi
		************
- Created addition bash script to compress file
	sub/compress_output.sh
*********************************
***	       README	      ***
*********************************

Date last updated: 14/02/2016

Brief overview:
------------------
full_enrichment_analysis is a bash script to process the Illumia sequencing data from enrichment (also called positive selection) screens in human hapoid (HAP1) cells. 
The input for this program is a fastq file, either decompressed or gzipped, and name for the screen. For comparison a previously alligned and analyzed control dataset is required. Specifically, the control dataset should refer to "*.uniqReads_mapped_UNIQ". It speaks for itself that this control file needs to be alligned against the same reference genome. Furthermore, one should be aware that this file only represents alligned reads and is not mapped to genes yet, therefore, the choice for the gene annotation file is fully flexible.

Finally, this folder also contains the 'additional_gene_mapping.sh'. Using this script, one can quick re-analyze a previously analzed enrichment screen with a different annotation file. This re-analysis script takes the directory where the analysis file reside as an input argument and reads the other paramters from the settings.conf file.

Dependencies:
------------------
- python
- bowtie (tested on v0.12.7 64-bit)
- intersectBed (tested on v2.11.1)
- R - for the one-sided Fisher exact test
- pigz (not required but allows for faster compression of intermediate files by using multiple cores, compared to gzip)
- pv
