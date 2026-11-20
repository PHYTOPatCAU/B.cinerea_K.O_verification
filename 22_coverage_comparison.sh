#!/bin/bash
#SBATCH --job-name=compare_cov
#SBATCH --output=compare_cov.out
#SBATCH --error=compare_cov.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

module load gcc12-env
module load miniconda3
conda activate samtools

region="5:1329646-1331550"
chr="5"
start=1329646
end=1331550

bam_dir="GenomesBGI/Bc/BAM/MD"
out_dir="GenomesBGI/Bc/GOI/coverage_analysis"

mkdir -p "$out_dir"

echo "=========================================="
echo "COVERAGE COMPARISON: WT vs KO"
echo "Region: $region"
echo "=========================================="
echo ""

for sample in B_05.10 WT WTB50 18_x_KO "23_x_KO_#7" 28xKO2 29XKO1; do
    bam="${bam_dir}/${sample}.DuplMark.bam"
    [ ! -f "$bam" ] && continue

    echo "Processing: $sample"

    cov_file="${out_dir}/${sample}_coverage.txt"

    samtools depth -r "$region" "$bam" > "$cov_file"

    awk -v sample="$sample" -v start="$start" -v end="$end" '
    BEGIN {
        sum = 0;
        count = 0;
        zero_cov = 0;
        low_cov = 0;
    }
    {
        depth = $3;
        sum += depth;
        count++;

        if (depth == 0) zero_cov++;
        else if (depth < 5) low_cov++;
    }
    END {
        expected_bases = end - start + 1;
        mean_cov = (count > 0) ? sum/count : 0;
        pct_zero = zero_cov * 100 / expected_bases;
        pct_low = low_cov * 100 / expected_bases;

        printf "  %s:\n", sample;
        printf "    Mean coverage: %.2f\n", mean_cov;
        printf "    Zero coverage: %.1f%% (%d/%d bp)\n", pct_zero, zero_cov, expected_bases;
        printf "    Low coverage (<5x): %.1f%%\n", pct_low;
        printf "\n";
    }
    ' "$cov_file"
done

echo ""
echo "=========================================="
echo "EXPECTED RESULTS"
echo "=========================================="
echo ""
echo "WT/B_05.10 samples:"
echo "  - Mean coverage: >20x (similar to genome-wide)"
echo "  - Zero coverage: <5%"
echo ""
echo "KO samples (18, 23, 28, 29):"
echo "  - Mean coverage: <5x (mostly at breakpoints)"
echo "  - Zero coverage: >80%"
echo ""
echo "If KO samples show significant coverage in middle:"
echo "  → Investigate those reads (likely artifacts)"
echo ""
