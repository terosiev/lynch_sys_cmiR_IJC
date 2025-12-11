#!/bin/bash -l
mkdir -p merged
readarray -t array < fastq.list
for filename in "${array[@]}" ; do
base_name=${filename%%_*}
sample_file="${filename%%.*}.fastq"
cat $sample_file>> merged/${base_name}_merged.fastq

done





