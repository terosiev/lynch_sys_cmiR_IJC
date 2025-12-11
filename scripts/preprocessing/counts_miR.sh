ls *.bam > merged.list
readarray -t array < merged.list
for filename in "${array[@]}" ; do
base_name=${filename%_*}
echo -e "id\tcounts" > $base_name.tsv
samtools view -q 10  $filename | awk '{print $3}' | grep -v '^*' | sort -k1 | uniq -c | awk '{print $2 "\t" $1}' >> $base_name.tsv
done

