#!/bin/bash -l

#SBATCH --job-name=NGS_preprocess
#SBATCH --output=NGS_out_%A_%a.txt
#SBATCH --error=NGS_err_%A_%a.txt
#SBATCH --time=04:00:00
#SBATCH --account=project_2002657
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-96
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=small 


# export PROJAPPL directory for sharing compiled applications and libraries etc with project members
export PROJAPPL="/projappl/project_2002657"

# set locale settings to avoid warnings if any especially with mac users
export LC_ALL="en_US.UTF-8"

# Load modules
module load biokit
module load bioconda/3

# Make necessary directories 
mkdir -p results-fastqc
mkdir -p results-fastx_clipped
mkdir -p results-fastx_trimmed
mkdir -p results-fastx_filtered

sample_file=$(sed -n "$SLURM_ARRAY_TASK_ID"p fastq.list)
base_name="${sample_file%.*}"

echo "$sample_file"
echo "$base_name"

#gunzip ${sample_file}
#sample_file="${sample_file%.*}"

echo "FASTQC"

	# FastQC
	# run quality control for the samples
	fastqc -o results-fastqc ${sample_file}

echo "fastx-toolkit"

	# FastX
	# Clip adapter sequences (19 nt), adapter must match at whole lenght, minimum read length 20 bp, keep only clipped reads
	fastx_clipper -a AACTGTAGGCACCATCAAT -l 20 -M 19 -c -v -i ${sample_file} -o results-fastx_clipped/${base_name}_clipped.fastq


	#FastX
	#Trim clipped reads to 22 bp (miRNAs), keep bases 1 - 22
	fastx_trimmer -v -f 1 -l 22 -i results-fastx_clipped/${base_name}_clipped.fastq -o results-fastx_trimmed/${base_name}_trimmed.fastq


	#FastX
	#quality filtering of the trimmed reads, cut-off 25 Phred score, minimum percentage 90%
	fastq_quality_filter -v -q 25 -p 90 -i results-fastx_trimmed/${base_name}_trimmed.fastq -o results-fastx_filtered/${base_name}_filtered.fastq
