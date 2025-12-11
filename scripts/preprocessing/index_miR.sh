#!/bin/bash -l

#SBATCH --job-name=NGS_index
#SBATCH --output=NGS_out_%j.txt
#SBATCH --error=NGS_err_%j.txt
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=8
#SBATCH --account=project_xxxxx
#SBATCH --mem=16000


# export PROJAPPL directory for sharing compiled applications and libraries etc with project members
export PROJAPPL="/projappl/project_xxxxx"

# set locale settings to avoid warnings if any especially with mac users
export LC_ALL="en_US.UTF-8"

#CPU correct usage
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Load modules (includes fastx-toolkit, bowtie, samtools etc.)
module load biokit

echo "fastx-toolkit"

	# Convert RNA-genome to DNA-genome ("U" -> "T")
	fasta_nucleotide_changer -d -i hsa_mir_rna.fa -o hsa_mir_dna.fa

echo "bowtie"

	# Bowtie
	# create indexes for reference genome
	bowtie-build -f hsa_mir_dna.fa mirna

