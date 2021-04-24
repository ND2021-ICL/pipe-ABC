#!/bin/sh

cd /project/data/neuronaldifferentiation2021/orange_pipeline/StringTie/human_count

echo "Starting at "`date` > StringTie_errors.log

module load stringtie/1.3.4c

for i in /project/data/neuronaldifferentiation2021/orange_pipeline/STAR/human_starmapped/*.bam
do 
name=$(echo "$i" | cut -c81-91)
echo "Processing:  "$name >> StringTie_errors.log
stringtie $i -p 2 -e -G gencode.v36.annotation.gtf -o ballgown/$name/$name.gtf
done

echo "finished at "`date` >> StringTie_errors.log
