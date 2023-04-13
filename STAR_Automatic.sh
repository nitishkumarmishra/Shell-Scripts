#!/bin/bash

#for PREFIX in $(ls |grep "fastq")
for PREFIX in $(ls |grep "fastq"|grep -v "trimming"| grep -v "_val_"|sed s/[12].fastq.gz// | sed 's/...........$//'|sort -u)
do
	#echo "****************** Sample $PREFIX submitted ******************"
	re="^([^.]+).(.*)$"
	[[ ${PREFIX} =~ $re ]] && var1="${BASH_REMATCH[1]}" && var2="${BASH_REMATCH[2]}"
	#echo $var1 $var2
	#echo $var2
	arrIN=(${var1//_/ })
	PREFIX=${arrIN[0]}_${arrIN[1]}
	echo "****************** Sample $PREFIX submitted *********************"
	#echo $PREFIX
	#bsub -q priority -n 12 -P RNAseq_Automatic -R "rusage[mem=25000]" -J RNAseq_STAR -oo %J.out -eo %J.err "sh STAR_count.sh $PREFIX"
	echo "*************** Gene counts for sample $PREFIX ******************"
	echo  "#########################################################################"	
	echo -e "#########################################################################""\n"
	
done

