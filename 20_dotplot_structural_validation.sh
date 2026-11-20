#!/bin/bash
#SBATCH --job-name=DotPlot_SV
#SBATCH --output=DotPlot_SV_%j.out
#SBATCH --error=DotPlot_SV_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

ref="GenomesBGI/Bc/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa"
assembly_dir="GenomesBGI/Bc/assemblies"
output_dir="GenomesBGI/Bc/dotplots"
mkdir -p "$output_dir"

target_chr="5"
target_start=1329671
target_end=1330257
flank=10000

analysis_mode="focused"

echo "=== Dotplot and Structural Variant Visualization ==="
echo "Target: ${target_chr}:${target_start}-${target_end}"
echo "Analysis mode: $analysis_mode"
echo "Flanking region: ${flank}bp"
echo ""

echo "Extracting extended target region from reference..."
samtools faidx "$ref"
extended_start=$((target_start - flank))
extended_end=$((target_end + flank))
extended_region="${target_chr}:${extended_start}-${extended_end}"
samtools faidx "$ref" "$extended_region" > "$output_dir/extended_target_ref.fa"

echo "Extended region: $extended_region"
echo "Analysis approach: $analysis_mode"
echo ""

generate_dotplot() {
    local sample=$1
    local assembly_file=$2
    local output_prefix="$output_dir/${sample}"

    echo "  Generating dotplot for $sample..."

    if [ "$analysis_mode" == "focused" ]; then
        echo "    Finding contigs that match target region..."
        nucmer --prefix="${output_prefix}_search" \
               "$output_dir/extended_target_ref.fa" "$assembly_file"

        if [ -f "${output_prefix}_search.delta" ]; then
            show-coords -T -H "${output_prefix}_search.delta" | \
                awk '{if(length($12)>0) print $12}' | sort -u > "${output_prefix}_contigs.list"

            if [ -s "${output_prefix}_contigs.list" ]; then
                echo "    Extracting matching contigs for focused analysis..."
                samtools faidx "$assembly_file"
                while read contig; do
                    samtools faidx "$assembly_file" "$contig" >> "${output_prefix}_focused.fa"
                done < "${output_prefix}_contigs.list"

                assembly_to_use="${output_prefix}_focused.fa"
            else
                echo "    No matching contigs found, using full assembly"
                assembly_to_use="$assembly_file"
            fi
        else
            assembly_to_use="$assembly_file"
        fi
    else
        assembly_to_use="$assembly_file"
    fi

    nucmer --prefix="$output_prefix" \
           --minmatch=15 \
           --mincluster=20 \
           --diagdiff=5 \
           --maxgap=90 \
           --coords \
           "$output_dir/extended_target_ref.fa" "$assembly_to_use"

    if [ -f "${output_prefix}.delta" ]; then
        mummerplot --fat --layout --filter --png \
                   --prefix="$output_prefix" \
                   "${output_prefix}.delta"

        mummerplot --fat --layout --filter --postscript \
                   --prefix="${output_prefix}_hq" \
                   "${output_prefix}.delta"
    fi

    show-coords -T -H "${output_prefix}.delta" > "${output_prefix}.coords"

    echo "    Checking for inversion signatures..."
    awk -v start="$target_start" -v end="$target_end" '
    NR>1 {
        ref_start = ($7 < $8) ? $7 : $8
        ref_end = ($7 < $8) ? $8 : $7

        if (ref_start <= end && ref_end >= start) {
            ref_forward = ($7 < $8) ? 1 : 0
            query_forward = ($1 < $2) ? 1 : 0

            if (ref_forward != query_forward) {
                print "POTENTIAL_INVERSION:", $0
            } else {
                print "NORMAL_ALIGNMENT:", $0
            }
        }
    }' "${output_prefix}.coords" > "${output_prefix}_inversion_analysis.txt"

    echo "    Dotplot complete for $sample"
}

for assembly in "$assembly_dir"/*_assembly; do
    sample=$(basename "$assembly" _assembly)
    echo "Processing sample: $sample"

    fasta=""
    if [ -f "$assembly/scaffolds.fasta" ]; then
        fasta="$assembly/scaffolds.fasta"
        echo "  Using: scaffolds.fasta"
    elif [ -f "$assembly/K77/scaffolds.fasta" ]; then
        fasta="$assembly/K77/scaffolds.fasta"
        echo "  Using: K77/scaffolds.fasta"
    elif [ -f "$assembly/K55/scaffolds.fasta" ]; then
        fasta="$assembly/K55/scaffolds.fasta"
        echo "  Using: K55/scaffolds.fasta"
    elif [ -f "$assembly/contigs.fasta" ]; then
        fasta="$assembly/contigs.fasta"
        echo "  Using: contigs.fasta"
    else
        echo "  No assembly found for $sample, skipping"
        continue
    fi

    generate_dotplot "$sample" "$fasta"
    echo ""
done

echo "=== Generating Comparative Analysis ==="
summary_file="$output_dir/inversion_summary.txt"
echo -e "Sample\tTotal_Alignments\tNormal_Alignments\tPotential_Inversions\tInversion_Percentage\tNotes" > "$summary_file"

for assembly in "$assembly_dir"/*_assembly; do
    sample=$(basename "$assembly" _assembly)
    analysis_file="$output_dir/${sample}_inversion_analysis.txt"

    if [ -f "$analysis_file" ]; then
        total=$(wc -l < "$analysis_file")
        normal=$(grep -c "NORMAL_ALIGNMENT" "$analysis_file" || echo 0)
        inversions=$(grep -c "POTENTIAL_INVERSION" "$analysis_file" || echo 0)

        if [ "$total" -gt 0 ]; then
            percentage=$(awk "BEGIN {printf \"%.1f\", ($inversions/$total)*100}")
        else
            percentage="NA"
        fi

        notes=""
        if [ "$inversions" -gt 0 ]; then
            notes="⚠️ INVERSION_DETECTED"
        elif [ "$total" -eq 0 ]; then
            notes="❌ NO_COVERAGE"
        else
            notes="✅ NORMAL"
        fi

        echo -e "${sample}\t${total}\t${normal}\t${inversions}\t${percentage}%\t${notes}" >> "$summary_file"
    else
        echo -e "${sample}\tNA\tNA\tNA\tNA\t❌ ANALYSIS_FAILED" >> "$summary_file"
    fi
done

cat > "$output_dir/plot_alignments.R" << 'EOF'
#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)

plot_sample_alignments <- function(sample_name, coords_file, output_dir) {
    if (!file.exists(coords_file)) return(NULL)

    tryCatch({
        data <- read.table(coords_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
        colnames(data) <- c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDEN",
                           "TAG", "FRAME", "REF", "QRY", "COV_R", "COV_Q")

        p <- ggplot(data, aes(x=S1, y=S2)) +
            geom_segment(aes(xend=E1, yend=E2, color=IDEN)) +
            scale_color_gradient(low="red", high="blue", name="Identity %") +
            labs(title=paste("Alignment plot:", sample_name),
                 x="Reference Position", y="Query Position") +
            theme_minimal() +
            theme(plot.title = element_text(hjust=0.5))

        ggsave(paste0(output_dir, "/", sample_name, "_alignment_plot.png"),
               plot=p, width=10, height=8, dpi=300)

        return(nrow(data))
    }, error=function(e) {
        cat("Error processing", sample_name, ":", e$message, "\n")
        return(0)
    })
}

coords_files <- list.files(".", pattern=".*\\.coords$", full.names=TRUE)

cat("Generating R plots for", length(coords_files), "samples...\n")

for (coords_file in coords_files) {
    sample_name <- gsub("\\.coords$", "", basename(coords_file))
    alignments <- plot_sample_alignments(sample_name, coords_file, ".")
    cat("  Processed", sample_name, "- alignments:", alignments, "\n")
}

cat("R plotting complete!\n")
EOF

echo "Running R visualization..."
cd "$output_dir"
chmod +x plot_alignments.R
Rscript plot_alignments.R

echo ""
echo "=== Analysis Complete ==="
echo "Output directory: $output_dir"
echo ""
echo "Files generated:"
echo "1. [sample].png/.ps - MUMmer dotplots for visual inspection"
echo "2. [sample]_inversion_analysis.txt - Detailed inversion analysis"
echo "3. inversion_summary.txt - Summary of all samples"
echo "4. [sample]_alignment_plot.png - R-generated alignment plots"
echo ""
echo "Interpretation:"
echo "- Normal alignment: Query and reference have same orientation"
echo "- Potential inversion: Query and reference have opposite orientations"
echo "- Look for samples with >50% potential inversions in your target region"
echo ""
echo "Visual inspection:"
echo "1. Open .png dotplot files - inversions appear as diagonal lines with opposite slopes"
echo "2. Compare wildtype vs knockout samples side by side"
echo "3. Look for structural differences in your target region: ${target_chr}:${target_start}-${target_end}"
