#!/bin/bash

#SBATCH --job-name=mpileup
#SBATCH --error=mpileup.err
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

REF=Botrytis_cinerea.ASM83294v1.dna.toplevel.fa

BAMFolder="GenomesBGI/Bc/BAM/test"

for bamfile in $BAMFolder/*.bam; do
  Sample=$(basename "$bamfile" .bam)
  output_file="$Sample.pileup"
  samtools mpileup -r $REF $bamfile > "$output_file"
done
