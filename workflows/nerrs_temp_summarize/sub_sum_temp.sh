#!/bin/bash

#SBATCH --time=1000
#SBATCH --mem=100gb

Rscript summarize_temp.R
