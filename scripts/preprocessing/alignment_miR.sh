#!/bin/bash -l

#SBATCH --job-name=NGS_alignment
#SBATCH --output=NGS_out_%j.txt
#SBATCH --error=NGS_err_%j.txt
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --cpus-per-task=2
#SBATCH --account=project_2002657
#SBATCH --mem=16000
#SBATCH --array=1-155


# export PROJAPPL directory for sharing compiled applications and libraries etc with project members
export PROJAPPL="/projappl/project_2002657"

# set locale settings to avoid warnings if any especially with mac users
export LC_ALL="en_US.UTF-8"

#CPU correct usage
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Load modules (includes fastx-toolkit, bowtie, samtools etc.)
module load biokit

# Make necessary directories 
mkdir -p results-bowtie

sample_file=$(sed -n "$SLURM_ARRAY_TASK_ID"p merged.list)
base_name="${sample_file%.*}"     

echo "$sample_file"
echo "$base_name"

echo "bowtie"

	# Align samples to mature-miRNA reference genome from miRBase v.21, 2 mismatches allowed in v-mode
	bowtie -v 2 -k 1 --best -p 1 -q -S mirna ${sample_file} > results-bowtie/${base_name}.sam  --un ${base_name}_notAligned.fastq
	samtools view -bS results-bowtie/${base_name}.sam > results-bowtie/${base_name}.bam

#echo "samtools"

	# Samtools
#	samtools sort results-bowtie/${base_name}.bam -o results-bowtie/${base_name}_sorted.bam
#	samtools index results-bowtie/${base_name}_sorted.bam
