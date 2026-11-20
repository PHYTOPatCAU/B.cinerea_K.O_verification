#!/bin/bash
#SBATCH --job-name=extract_loci
#SBATCH --output=extract_loci.out
#SBATCH --error=extract_loci.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

bam_dir="GenomesBGI/Bc/BAM/MD"
loci_file="GenomesBGI/Bc/GOI/loci.txt"
out_dir="GenomesBGI/Bc/GOI/BAM_loci"
mkdir -p "$out_dir"

samples=("18_x_KO" "23_x_KO_#7" "28xKO2" "29XKO1" "B_05.10")

declare -A loci
flank=500
while IFS=$'\t' read -r region name; do
    [[ "$region" == "Locus" ]] && continue
    [[ -z "$region" ]] && continue

    chr=$(echo "$region" | cut -d: -f1)
    start=$(echo "$region" | cut -d: -f2 | cut -d- -f1)
    end=$(echo "$region" | cut -d: -f2 | cut -d- -f2)

    new_start=$((start - flank))
    [ $new_start -lt 1 ] && new_start=1
    new_end=$((end + flank))

    new_region="${chr}:${new_start}-${new_end}"
    loci["$name"]="$new_region"
    echo "  Loaded: $name -> $region (with ${flank}bp flank: $new_region)"
done < "$loci_file"

echo ""
echo "=== Extracting ${#loci[@]} loci from ${#samples[@]} samples ==="
echo "Loci loaded:"
for name in "${!loci[@]}"; do
    echo "  $name: ${loci[$name]}"
done
echo ""

for sample in "${samples[@]}"; do
    input_bam="$bam_dir/${sample}.DuplMark.bam"

    if [ ! -f "$input_bam" ]; then
        echo "WARNING: $input_bam not found, skipping"
        continue
    fi

    echo "Processing $sample..."

    regions=""
    for name in "${!loci[@]}"; do
        regions="$regions ${loci[$name]}"
    done

    output_bam="$out_dir/${sample}.loci.bam"
    temp_bam="$out_dir/${sample}.loci.unsorted.bam"

    samtools view -b "$input_bam" $regions > "$temp_bam"
    samtools sort -o "$output_bam" "$temp_bam"
    samtools index "$output_bam"
    rm -f "$temp_bam"

    echo "  Created: $output_bam"
done

echo ""
echo "=== Creating per-locus BAMs (all samples) ==="

for name in "${!loci[@]}"; do
    region="${loci[$name]}"
    locus_dir="$out_dir/$name"
    mkdir -p "$locus_dir"

    echo "Extracting locus: $name ($region)"

    for sample in "${samples[@]}"; do
        input_bam="$bam_dir/${sample}.DuplMark.bam"
        [ ! -f "$input_bam" ] && continue

        output_bam="$locus_dir/${sample}.${name}.bam"
        samtools view -b "$input_bam" "$region" > "$output_bam"
        samtools index "$output_bam"
    done
done

echo ""
echo "=== Done! ==="
echo "Per-sample BAMs (all loci): $out_dir/*.loci.bam"
echo "Per-locus BAMs (all samples): $out_dir/<locus_name>/*.bam"
