#!/bin/bash
#SBATCH --job-name=samplot_compare
#SBATCH --output=samplot_compare.out
#SBATCH --error=samplot_compare.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samplot

ref="GenomesBGI/Bc/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa"
bam_dir="GenomesBGI/Bc/BAM/MD"
loci_file="GenomesBGI/Bc/GOI/loci.txt"
out_dir="GenomesBGI/Bc/GOI/samplot_comparisons"
mkdir -p "$out_dir"

flank=500

mutants=("18_x_KO" "23_x_KO_#7" "28xKO2" "29XKO1")
wt_bam="${bam_dir}/B_05.10.DuplMark.bam"

while IFS=$'\t' read -r region name; do
    [[ "$region" == "Locus" ]] && continue
    [[ -z "$region" ]] && continue

    chr=$(echo "$region" | cut -d: -f1)
    start=$(echo "$region" | cut -d: -f2 | cut -d- -f1)
    end=$(echo "$region" | cut -d: -f2 | cut -d- -f2)

    new_start=$((start - flank))
    [ $new_start -lt 1 ] && new_start=1
    new_end=$((end + flank))

    mkdir -p "$out_dir/$name"

    for mutant in "${mutants[@]}"; do
        mutant_bam="${bam_dir}/${mutant}.DuplMark.bam"
        [ ! -f "$mutant_bam" ] && continue

        output_png="$out_dir/$name/${mutant}_vs_WT.png"

        echo "Creating: $name - $mutant vs WT"

        samplot plot \
            -c "$chr" \
            -s "$new_start" \
            -e "$new_end" \
            -o "$output_png" \
            -r "$ref" \
            -T "$name: $mutant vs B_05.10" \
            -b "$mutant_bam" "$wt_bam" \
            -n "$mutant" "B_05.10_WT"
    done
done < "$loci_file"

echo "=== Done! Figures in: $out_dir/ ==="
