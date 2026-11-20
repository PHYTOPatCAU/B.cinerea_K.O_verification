#!/bin/bash
#SBATCH --job-name=samplot
#SBATCH --output=samplot.out
#SBATCH --error=samplot.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00

module load gcc12-env
module load miniconda3
conda activate samplot

echo "Using samplot: $(which samplot)"
echo "Python version: $(python --version)"
echo ""

ref="GenomesBGI/Bc/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa"
bam_dir="GenomesBGI/Bc/BAM/MD"
loci_file="GenomesBGI/Bc/GOI/loci.txt"
out_dir="GenomesBGI/Bc/GOI/samplot_figures"
annotation_bed="GenomesBGI/Bc/GOI/genes.bed.gz"
mkdir -p "$out_dir"

bam_18="${bam_dir}/18_x_KO.DuplMark.bam"
bam_23="${bam_dir}/23_x_KO_#7.DuplMark.bam"
bam_28="${bam_dir}/28xKO2.DuplMark.bam"
bam_29="${bam_dir}/29XKO1.DuplMark.bam"
bam_b05="${bam_dir}/B_05.10.DuplMark.bam"

flank=500

echo "=== Generating Samplot figures for all loci ==="
echo ""

while IFS=$'\t' read -r region name; do
    [[ "$region" == "Locus" ]] && continue
    [[ -z "$region" ]] && continue

    chr=$(echo "$region" | cut -d: -f1)
    start=$(echo "$region" | cut -d: -f2 | cut -d- -f1)
    end=$(echo "$region" | cut -d: -f2 | cut -d- -f2)

    new_start=$((start - flank))
    [ $new_start -lt 1 ] && new_start=1
    new_end=$((end + flank))

    echo "Processing: $name (${chr}:${new_start}-${new_end})"

    output_png="$out_dir/${name}_chr${chr}_${start}-${end}.png"

    if samplot plot \
        -c "$chr" \
        -s "$new_start" \
        -e "$new_end" \
        -o "$output_png" \
        -r "$ref" \
        -t "DEL" \
        --max_depth 0 \
        -A "$annotation_bed" \
        --annotation_filenames "Genes" \
        -b "$bam_18" "$bam_23" "$bam_28" "$bam_29" "$bam_b05" \
        -n "18_x_KO" "23_x_KO_#7" "28xKO2" "29XKO1" "B_05.10" 2>&1; then
        echo "  ✓ Created: $output_png"
    else
        echo "  ✗ FAILED: $name - skipping and continuing..."
    fi
    echo ""

done < "$loci_file"

echo ""
echo "=== Done! ==="
echo "Figures saved to: $out_dir/"
