#!/bin/bash
#SBATCH --job-name=InversionEvidence
#SBATCH --output=InversionEvidence_%j.out
#SBATCH --error=InversionEvidence_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

ref="GenomesBGI/Bc/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa"
bam_dir="GenomesBGI/Bc/BAM/MD"
output_dir="GenomesBGI/Bc/inversion_evidence"
mkdir -p "$output_dir"

target_chr="5"
target_start=1329671
target_end=1330257
flank=5000

echo "=== Detailed Inversion Evidence Analysis ==="
echo "Target locus: ${target_chr}:${target_start}-${target_end}"
echo "Analysis region: ${target_chr}:$((target_start-flank))-$((target_end+flank))"
echo ""

analyze_inversion_evidence() {
    local sample=$1
    local bam_file=$2
    local output_prefix="$output_dir/${sample}"

    echo "Analyzing sample: $sample"

    if [ ! -f "$bam_file" ]; then
        echo "  ❌ BAM file not found: $bam_file"
        return
    fi

    region="${target_chr}:$((target_start-flank))-$((target_end+flank))"

    echo "  📊 Extracting reads from region: $region"
    samtools view -h "$bam_file" "$region" > "${output_prefix}_region.sam"

    total_reads=$(samtools view -c "${output_prefix}_region.sam")
    echo "  📈 Total reads in region: $total_reads"

    if [ "$total_reads" -eq 0 ]; then
        echo "  ❌ No reads found in region - possible deletion"
        echo -e "${sample}\t0\t0\t0\t0\t0\t0\tDELETION_NO_READS" >> "$output_dir/inversion_evidence_summary.txt"
        return
    fi

    echo "  🔍 Analyzing read pair orientations..."

    samtools view "${output_prefix}_region.sam" | awk -v sample="$sample" -v target_start="$target_start" -v target_end="$target_end" '
    BEGIN {
        normal_pairs = 0
        inward_pairs = 0
        outward_pairs = 0
        same_strand_fwd = 0
        same_strand_rev = 0
        large_insert = 0
        split_reads = 0
        total_reads = 0
        discordant = 0
    }
    {
        if (and($2, 4) || and($2, 8)) next

        total_reads++
        pos = $4
        flag = $2
        isize = $9
        if (isize < 0) isize = -isize

        if (!and($2, 1)) next

        in_target = (pos >= target_start && pos <= target_end)

        read_reverse = and(flag, 16)
        mate_reverse = and(flag, 32)

        if (read_reverse && mate_reverse) {
            same_strand_rev++
            if (in_target) discordant++
        } else if (!read_reverse && !mate_reverse) {
            same_strand_fwd++
            if (in_target) discordant++
        } else {
            if (read_reverse) {
                if (isize > 1000) large_insert++
                else inward_pairs++
            } else {
                if (isize > 1000) large_insert++
                else normal_pairs++
            }
        }

        if (and(flag, 2048)) split_reads++
    }
    END {
        total_pairs = normal_pairs + inward_pairs + outward_pairs + same_strand_fwd + same_strand_rev
        same_strand_total = same_strand_fwd + same_strand_rev

        if (total_pairs > 0) {
            same_strand_pct = (same_strand_total / total_pairs) * 100
            large_insert_pct = (large_insert / total_pairs) * 100
            discordant_pct = (discordant / total_reads) * 100
        } else {
            same_strand_pct = 0
            large_insert_pct = 0
            discordant_pct = 0
        }

        evidence = "NORMAL"
        if (same_strand_pct > 30) {
            evidence = "INVERSION_EVIDENCE"
        } else if (large_insert_pct > 50 || total_reads < 10) {
            evidence = "DELETION_EVIDENCE"
        } else if (discordant_pct > 20) {
            evidence = "COMPLEX_SV"
        }

        print sample "\t" total_reads "\t" normal_pairs "\t" same_strand_total "\t" large_insert "\t" split_reads "\t" same_strand_pct "\t" evidence
    }' > "${output_prefix}_orientation_analysis.txt"

    echo "  🔍 Extracting discordant read pairs..."
    samtools view "${output_prefix}_region.sam" | awk '
    {
        flag = $2
        isize = $9
        if (isize < 0) isize = -isize

        if (and(flag, 4) || and(flag, 8)) next
        if (!and(flag, 1)) next

        read_reverse = and(flag, 16)
        mate_reverse = and(flag, 32)

        discordant = 0
        if (read_reverse && mate_reverse) discordant = 1
        if (!read_reverse && !mate_reverse) discordant = 1
        if (isize > 2000) discordant = 1

        if (discordant) {
            orientation = "UNKNOWN"
            if (read_reverse && mate_reverse) orientation = "BOTH_REV"
            else if (!read_reverse && !mate_reverse) orientation = "BOTH_FWD"
            else if (isize > 2000) orientation = "LARGE_INSERT"

            print $1 "\t" $3 "\t" $4 "\t" $6 "\t" isize "\t" orientation "\t" $10
        }
    }' > "${output_prefix}_discordant_reads.txt"

    echo "  ✂️  Analyzing soft-clipped reads..."
    samtools view "${output_prefix}_region.sam" | awk '
    BEGIN {
        left_clips = 0
        right_clips = 0
        both_clips = 0
        total_clips = 0
    }
    {
        cigar = $6

        left_clip = 0
        right_clip = 0

        if (match(cigar, /^([0-9]+)S/)) {
            left_clip = substr(cigar, RSTART, RLENGTH-1)
        }
        if (match(cigar, /([0-9]+)S$/)) {
            right_clip = substr(cigar, RSTART, RLENGTH-1)
        }

        if (left_clip > 10 || right_clip > 10) {
            total_clips++
            if (left_clip > 10 && right_clip > 10) both_clips++
            else if (left_clip > 10) left_clips++
            else right_clips++

            print $1 "\t" $3 "\t" $4 "\t" cigar "\t" left_clip "\t" right_clip "\t" $10
        }
    }
    END {
        print "SOFT_CLIP_SUMMARY:" >> "/dev/stderr"
        print "  Left clips: " left_clips >> "/dev/stderr"
        print "  Right clips: " right_clips >> "/dev/stderr"
        print "  Both clips: " both_clips >> "/dev/stderr"
        print "  Total clipped reads: " total_clips >> "/dev/stderr"
    }' > "${output_prefix}_soft_clipped.txt" 2>> "${output_prefix}_clip_summary.txt"

    echo "  ✅ Analysis complete for $sample"
    echo ""
}

echo -e "Sample\tTotal_Reads\tNormal_Pairs\tSame_Strand_Pairs\tLarge_Insert_Pairs\tSplit_Reads\tSame_Strand_Percent\tEvidence_Type" > "$output_dir/inversion_evidence_summary.txt"

for sample in 18_x_KO "23_x_KO_#7" 28xKO2 29XKO1 B_05.10 WT WTB50; do
    bam_file="${bam_dir}/${sample}.DuplMark.bam"
    analyze_inversion_evidence "$sample" "$bam_file"

    if [ -f "$output_dir/${sample}_orientation_analysis.txt" ]; then
        cat "$output_dir/${sample}_orientation_analysis.txt" >> "$output_dir/inversion_evidence_summary.txt"
    fi
done

echo "=== Generating Comparison Report ==="
comparison_file="$output_dir/sample_comparison.txt"

echo "Sample Comparison for Locus ${target_chr}:${target_start}-${target_end}" > "$comparison_file"
echo "=================================================================" >> "$comparison_file"
echo "" >> "$comparison_file"

echo "Legend:" >> "$comparison_file"
echo "- NORMAL: Standard read pair orientations" >> "$comparison_file"
echo "- INVERSION_EVIDENCE: >30% same-strand read pairs (→→ or ←←)" >> "$comparison_file"
echo "- DELETION_EVIDENCE: Low coverage or large insert sizes" >> "$comparison_file"
echo "- COMPLEX_SV: Mixed evidence suggesting complex structural variant" >> "$comparison_file"
echo "" >> "$comparison_file"

sort -k8,8 "$output_dir/inversion_evidence_summary.txt" >> "$comparison_file"

echo "" >> "$comparison_file"
echo "Detailed Analysis Files:" >> "$comparison_file"
echo "- [sample]_region.sam: Raw alignments in target region" >> "$comparison_file"
echo "- [sample]_orientation_analysis.txt: Read pair orientation statistics" >> "$comparison_file"
echo "- [sample]_discordant_reads.txt: Discordant read pairs" >> "$comparison_file"
echo "- [sample]_soft_clipped.txt: Soft-clipped reads (potential breakpoints)" >> "$comparison_file"

echo ""
echo "=== Analysis Complete ==="
echo "Results saved to: $output_dir"
echo ""
echo "Key findings will be in:"
echo "1. $comparison_file - Summary comparison"
echo "2. inversion_evidence_summary.txt - Numerical evidence"
echo ""
echo "Look for samples with INVERSION_EVIDENCE to identify your problematic knockout!"
