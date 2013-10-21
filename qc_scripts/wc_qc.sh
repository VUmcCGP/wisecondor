# Run on directory with [samplename].pickle, [samplename].bam and [samplename].bai

for SAMPLE in $1/*.pickle
do
	NAME_SAMPLE=$(basename "$SAMPLE")
	EXT_SAMPLE=${NAME_SAMPLE#*.}
	NAME_SAMPLE="${NAME_SAMPLE%.*$EXT_SAMPLE}"
	NAME_SAMPLE=$(echo $NAME_SAMPLE | cut -d _ -f 1)
	echo -e '\e[0m\tProcessing stats for sample:\t' $NAME_SAMPLE
	samtools flagstat $1/$NAME_SAMPLE.bam > $1/$NAME_SAMPLE.sort.bam.flagstat
	echo 'In sort.bam mapped with MQ of 0: ' >> $1/$NAME_SAMPLE.sort.bam.flagstat 
	samtools view -F 4 $1/$NAME_SAMPLE.bam | awk '$5=="0" ' | wc -l >> $1/$NAME_SAMPLE.sort.bam.flagstat
	echo 'Number of reads left after RETRO filter: ' >> $1/$NAME_SAMPLE.sort.bam.flagstat
	python ./calcTowers.py $1/$NAME_SAMPLE.pickle >> $1/$NAME_SAMPLE.sort.bam.flagstat
done
