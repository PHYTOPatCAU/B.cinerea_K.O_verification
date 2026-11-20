#!/bin/bash
#SBATCH --job-name=BcnovoAsml
#SBATCH --output=BcnovoAsml.out
#SBATCH --error=BcnovoAsml.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# path to the conda env
conda activate samtools

input_dir="GenomesBGI/Bc/fastq"
output_dir="GenomesBGI/Bc/assemblies"
mkdir -p "$output_dir"

echo "=== Starting SPAdes assembly for all samples ==="
for r1 in "$input_dir"/*_1.fq.gz; do
    sample=$(basename "$r1" _1.fq.gz)
    r2="$input_dir/${sample}_2.fq.gz"
    assembly_dir="$output_dir/${sample}_assembly"

    if [ -f "$r2" ]; then
        if [ -d "$assembly_dir" ] && [ -f "$assembly_dir/scaffolds.fasta" ]; then
            echo "Skipping $sample - assembly already exists at $assembly_dir"
        else
            echo "Assembling $sample..."
            spades.py -1 "$r1" -2 "$r2" -o "$assembly_dir" -k 33,55,77 -t 8
        fi
    else
        echo "Warning: Missing R2 file for $sample"
    fi
done
echo "=== SPAdes assembly complete ==="

