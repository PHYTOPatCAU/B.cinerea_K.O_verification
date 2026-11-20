module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

bam_directory="GenomesBGI/Bc/BAM/MD"
regions_file="GenomesBGI/Bc/BAM/ROI/regions.txt"
output_directory="GenomesBGI/Bc/BAM/ROI/"

while IFS=$'\t' read -r chromosome start_position end_position
do
    for bam_file in "${bam_directory}"/*.bam
    do
        sample=$(basename "$bam_file" .bam)
        output_bam="${output_directory}/extracted_${sample}_${chromosome}_${start_position}-${end_position}.bam"
        samtools view -b "$bam_file" "${chromosome}:${start_position}-${end_position}" > "$output_bam"
        echo "Region ${chromosome}:${start_position}-${end_position} extracted from ${sample} to ${output_bam}"
    done
done< "$regions_file"
