#!/bin/bash
#SBATCH --job-name=mummer
#SBATCH --error=mummer_%j.err
#SBATCH --output=mummer_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --qos=long
#SBATCH --time=5-00:00:00
#SBATCH --mail-type=FAIL,END

module load gcc12-env
module load miniconda3
conda activate samtools

mummer_path="mummer4"
assembly_dir="GenomesBGI/Bc/assemblies"

ref_fasta="GenomesBGI/Bc/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa"
mafft_outdir="GenomesBGI/Bc/mafft_results"
mkdir -p "$mafft_outdir"

samples=("18_x_KO" "23_x_KO_#7" "29X1" "B_05.10")
for chr in $(awk '/^>/ {print substr($0,2)}' "$ref_fasta"); do
    awk -v chr="$chr" '/^>/{p=($0==">"chr)} p' "$ref_fasta" > "$mafft_outdir/${chr}_ref.fa"
    for sample in "${samples[@]}"; do
        asm="$assembly_dir/${sample}_assembly/K55/scaffolds.fasta"
        cp "$asm" "$mafft_outdir/${chr}_${sample}.fa"
    done
    cat "$mafft_outdir/${chr}_ref.fa" "$mafft_outdir/${chr}_B_05.10.fa" "$mafft_outdir/${chr}_18_x_KO.fa" "$mafft_outdir/${chr}_23_x_KO_#7.fa" "$mafft_outdir/${chr}_29X1.fa" > "$mafft_outdir/${chr}_combined.fa"
    mafft --auto "$mafft_outdir/${chr}_combined.fa" > "$mafft_outdir/${chr}_aligned.fa"
done

outdir="GenomesBGI/Bc/mummer_results"
mkdir -p "$outdir"
mutants=("18_x_KO" "23_x_KO_#7" "29X1")
wt="B_05.10"

for chr in $(awk '/^>/ {print substr($0,2)}' "$ref_fasta"); do
    ref_chr="$mafft_outdir/${chr}_ref.fa"
    for sample in "${samples[@]}"; do
        scaffolds="$mafft_outdir/${chr}_${sample}.fa"
        if [ -s "$ref_chr" ] && [ -s "$scaffolds" ]; then
            prefix="$outdir/${chr}_${sample}_vs_ref"
            "$mummer_path"/nucmer -p "$prefix" "$ref_chr" "$scaffolds"
            "$mummer_path"/show-coords "${prefix}.delta" > "${prefix}.coords"
            "$mummer_path"/mummerplot -p "$prefix" -l "${prefix}.coords"
        fi
    done
done
