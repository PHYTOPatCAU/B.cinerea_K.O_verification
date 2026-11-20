#!/bin/bash

#SBATCH --job-name=intercept
#SBATCH --error=intercept.err
#SBATCH --nodes=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

REF=Botrytis_cinerea.ASM83294v1.dna.toplevel.fa

BAMFolder="GenomesBGI/Bc/BAM/test"
cd $BAMFolder

bedtools intersect -v -a 23_x_KO.bed -b B_05.10.bed > 23kodel.bed
bedtools intersect -v -a 18_x_KO.bed -b B_05.10.bed > 18kodel.bed
bedtools intersect -v -a 23_x_KO.bed -b 18_x_KO.bed > 23_18kodel.bed
