#!/bin/bash


#############################################################
# A tool for finding the optimal seq. of the genes on xaxis #		
#############################################################
# authors:													#
# Elmer stickel												#
# Netherlands Cancer Insitute								#
# 22th February 2016										#
# Contact: e.stickel@gmail.com								#
# Commercial use of the software is prohibited without 		#
# prior approval by Elmer stickel 							#
#############################################################


# Now here starts the real fun! :)
# We want to find the best way of ploting the genes on the x-axis
# But what is the best?
# That plot were the most signinficant genes do not overlap with each other of course
# How can we define not-overlapping? By taking the ten most significant hits and finding that sequence of genes were they all have the same x-axis distance between each other.
# But what sequence does yield the most similar distance between those genes?
# There's no formula for it buuuuut we do have a fast server so lets try it many many times!
# Then we can evaluate the outcome and pick that sequence that has the largest average distance between the top hits and the least spread
# Real fun would be to do this on a GPU but the build-in Intel GPU in the XEON CPU we have in our server isn't of much help 


# Number of genes to take into account for calculating non-overlap
declare -i gc=5

# Determine the total number of genes that have insertions (totins) and find top $gc hits
declare -i totins=$(wc -l $1 | awk '{print $1}' )
#printf "Total genes to be plotted in x-axis: $totins\n"
declare -i minwidth=$(echo | awk -v totins=$totins 'BEGIN {rounded = sprintf("%.0f", totins*0.8); print rounded}')
#printf "Minimum width ('in genes') of the top hits: $minwidth\n"
declare -i mindist=$(echo | awk -v totins=$totins -v gc=$gc 'BEGIN {rounded = sprintf("%.0f", totins/gc*0.2); print rounded}')
#printf "Minimum distance ('in genes') between each top hit: $mindist\n"
declare -i SD=$totins
# Add an extra column with id's to the p-value, then extract that column and use it for shuffing. Why, because shuffling with integers is way faster than shuffeling with strings (ie. genenames).
# By adding the id column we can later use this column to join the desired positions to it
awk '{print NR"\t"$0}' $1 > $1.id
awk '{print $1}' $1.id > $1.ids_only
# Get the id's of top $gc most significant hits (thits) we for sure don't want to see overlapping.
sort -g -k7,7 $1.id | head -n $gc | awk '{print $1}' > $1.thits

# A function to shuffle the elements of the array
shuffle() {
   local i tmp size max rand

   # $RANDOM % (i+1) is biased because of the limited range of $RANDOM
   # Compensate by using a range which is a multiple of the array size.
   size=${#array[*]}
   max=$(( 32768 / size * size ))

   for ((i=size-1; i>0; i--)); do
      while (( (rand=$RANDOM) >= max )); do :; done
      rand=$(( rand % (i+1) ))
      tmp=${array[i]} array[i]=${array[rand]} array[rand]=$tmp
   done
}

# Put all the ID's in an array
declare -a array=( $( cat $1.ids_only ) )
declare -a bestarray
# Make sure that diff, prev and sum are integers
declare -i diff
declare -i prev
declare -i sum

# We iterate for x times to find the best sequence. Lets say 1000 times?
for iter in {1..500}
do
	# Call the shuffle function to shuffle each time we loop though iter
	shuffle
	# Now we write the shuffled array into a new file, I'm not happy with this because it causes a lot of non-needed disk writes.
	# The file has two column, the first is the random/shuffled order of the elements and the second column is just the number of the row
	# This first column therefore represents the gene and the second column the new position of this gene
	# The second column we can then use to determine the 'distance' between the top-x genes
	for l in "${array[@]}"; do printf $l"\n"; done | awk '{print $0"\t"NR}' > $1.arr.tsv.tmp

	# Now we take the x-top hits we previously defined (thits) and look up these genes in the newly shuffled array. And build a new array from it (evaulatethis)
	declare -a evaluatethis=( $(awk 'NR==FNR { a[$1]=$1; next} $1 in a {print $2}' $1.thits $1.arr.tsv.tmp))
	prev=0

	# Now we compare the distance between each element of the new array evaluatethis (distance we call diff)
	declare -a diffarray=( $( for j in "${evaluatethis[@]}"; do diff=$((j-prev)); printf $diff"\n"; prev=$j; done))
	# Before doing the artitmetic we check if distance is never smaller than 20% of the optimal distance between each datapoint (defined as minwidth/total tophits)
	# If it is smaller than set STOP at 1 and prevent execution of subsequent calculations
	stop=0
	# Loop through current shuffled array
	for p in "${diffarray[@]}"
		do
    		if [ "$p" -lt "$mindist" ] ; then
        		stop=1
    		fi
	done
	# Check STOP variable
	if [ "$stop" -eq 1 ]
		then
		printf "."
	else
		for j in "${evaluatethis[@]}"; do diff=$((j-prev)); prev=$j; sum=$((sum+diff)); done # We can't use sum in diffarray because it resides in its on env. and is not a global var.
		# Take the sum of all distances between the top datapoints and check of the at least make up 80% of total width
		if [ "$sum" -ge "$minwidth" ]; then
				# Calculate standard deviation in AWK
				newSDfloat=( $(for k in "${diffarray[@]}"; do printf $k"\n"; done | awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}'))
				# Convert to int
				declare -i newSD=${newSDfloat/.*}
				# Test of new SD is better than previous
				if [ "$newSD" -lt "$SD" ]
					then 
						# If so, set newly found SD as the next value to test against and set bestarray to the current array
						SD=$newSD
						bestarray=("${array[@]}")
						printf "*" 
				else
					printf "."
				fi
			else
				printf "."
		fi

	fi
rm $1.arr.tsv.tmp
sum=0
prev=0
diff=0
done
for n in "${bestarray[@]}"; do printf $n"\n"; done > $1.only_optseq.txt
paste $1 $1.only_optseq.txt > $1.optseq.txt
mv $1.optseq.txt $1
rm $1.only_optseq.txt
rm $1.id
rm $1.ids_only