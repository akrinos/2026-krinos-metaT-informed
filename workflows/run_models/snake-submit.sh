#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=5gb
#SBATCH -p batch

snakemake -c 1 --unlock
snakemake \
	 --use-conda \
	 --rerun-triggers mtime \
	 --rerun-incomplete	-j 33 \
         --executor slurm \
         --cluster-generic-submit-cmd 'sbatch -N 1 -J {rule} --time=16:00:00 -p batch'
