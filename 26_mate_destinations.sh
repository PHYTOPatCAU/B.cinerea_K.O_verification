#!/bin/bash
#SBATCH --job-name=mate_destinations
#SBATCH --output=mate_destinations_%j.out
#SBATCH --error=mate_destinations_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL,END

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

region="5:1329646-1331550"
chr="5"
del_start=1329646
del_end=1331550

bam_dir="GenomesBGI/Bc/BAM/MD"
out_dir="GenomesBGI/Bc/GOI/mate_destinations"
out_file="${out_dir}/mate_summary.txt"
mkdir -p "$out_dir"

echo "==========================================" | tee "$out_file"
echo "MATE DESTINATION ANALYSIS" | tee -a "$out_file"
echo "Where do mates of reads in deletion go?" | tee -a "$out_file"
echo "Deletion: $region" | tee -a "$out_file"
echo "==========================================" | tee -a "$out_file"
echo "" | tee -a "$out_file"

for sample in 18_x_KO "23_x_KO_#7" 28xKO2 29XKO1 B_05.10; do
    bam="${bam_dir}/${sample}.DuplMark.bam"
    [ ! -f "$bam" ] && echo "WARNING: $bam not found, skipping" | tee -a "$out_file" && continue

    echo "=== Sample: $sample ===" | tee -a "$out_file"

    samtools view -F 256 -f 1 "$bam" "$region" | \
    awk -v del_chr="$chr" -v del_start="$del_start" -v del_end="$del_end" '
    {
        flag  = $2
        rnext = $7
        pnext = $8
        isize = $9
        if (isize < 0) isize = -isize

        mate_chr = (rnext == "=") ? del_chr : rnext

        bin_start = int(pnext / 500) * 500
        bin_end   = bin_start + 500
        bin_label = mate_chr ":" bin_start "-" bin_end

        if (mate_chr == del_chr && pnext >= del_start && pnext <= del_end) {
            class = "INSIDE_DELETION"
        } else if (mate_chr == del_chr && \
                   pnext >= del_start - 5000 && pnext <= del_end + 5000) {
            class = "NEAR_DELETION (±5kb)"
        } else if (mate_chr == del_chr) {
            class = "SAME_CHR_DISTAL"
        } else if (rnext == "*") {
            class = "MATE_UNMAPPED"
        } else {
            class = "DIFFERENT_CHR (" mate_chr ")"
        }

        count_bin[bin_label]++
        count_class[class]++
        isize_bin[bin_label] += isize
        total++
    }
    END {
        print "  Total reads in deletion region:", total

        print ""
        print "  --- Mate location class summary ---"
        for (c in count_class)
            printf "  %-35s : %d reads (%.1f%%)\n", c, count_class[c], count_class[c]*100/total

        print ""
        print "  --- Top 15 mate destination windows (500bp bins) ---"
        print "  (Exclude bins overlapping deletion itself)"
        n = 0
        for (b in count_bin) {
            bins[n] = b
            cnts[n] = count_bin[b]
            n++
        }
        for (i = 0; i < n; i++)
            for (j = i+1; j < n; j++)
                if (cnts[j] > cnts[i]) {
                    tmp=cnts[i]; cnts[i]=cnts[j]; cnts[j]=tmp
                    tmp=bins[i]; bins[i]=bins[j]; bins[j]=tmp
                }
        printed = 0
        for (i = 0; i < n && printed < 15; i++) {
            split(bins[i], b, /[:-]/)
            bstart = b[2]+0; bend = b[3]+0
            if (b[1] == del_chr && bstart >= del_start-500 && bend <= del_end+500) continue
            mean_is = (cnts[i] > 0) ? isize_bin[bins[i]]/cnts[i] : 0
            printf "  %4d reads → %-35s (mean insert: %.0fbp)\n", \
                cnts[i], bins[i], mean_is
            printed++
        }
    }
    ' | tee -a "$out_file"

    echo "" | tee -a "$out_file"

    echo "  --- Candidate reinsertion sites (>=3 mates in 500bp window, outside deletion) ---" | tee -a "$out_file"
    samtools view -F 256 -f 1 "$bam" "$region" | \
    awk -v del_chr="$chr" -v del_start="$del_start" -v del_end="$del_end" '
    {
        rnext = $7
        pnext = $8
        mate_chr = (rnext == "=") ? del_chr : rnext
        if (rnext == "*") next
        bin_start = int(pnext / 500) * 500
        bin_label = mate_chr ":" bin_start "-" (bin_start+500)
        count[bin_label]++
    }
    END {
        found = 0
        for (b in count) {
            split(b, parts, /[:-]/)
            bstart = parts[2]+0
            if (parts[1] == del_chr && bstart >= del_start-500 && bstart <= del_end+500) continue
            if (count[b] >= 3) {
                printf "  *** %4d mates → %s  ← INVESTIGATE THIS\n", count[b], b
                found++
            }
        }
        if (found == 0) print "  None found (no window with >=3 mates outside deletion)"
    }
    ' | tee -a "$out_file"

    echo "" | tee -a "$out_file"
    echo "----------------------------------------" | tee -a "$out_file"
    echo "" | tee -a "$out_file"
done

echo "" | tee -a "$out_file"
echo "==========================================" | tee -a "$out_file"
echo "INTERPRETATION" | tee -a "$out_file"
echo "==========================================" | tee -a "$out_file"
echo "" | tee -a "$out_file"
echo "A real reinsertion ~2kb away should produce:" | tee -a "$out_file"
echo "  1. A cluster of mates at ONE specific window" | tee -a "$out_file"
echo "  2. The SAME cluster in multiple KO samples" | tee -a "$out_file"
echo "  3. Insert sizes consistent with library size (~300-600bp)" | tee -a "$out_file"
echo "  4. B_05.10/WT should NOT show the same cluster" | tee -a "$out_file"
echo "" | tee -a "$out_file"
echo "If candidates found, check those coordinates in IGV" | tee -a "$out_file"
echo "" | tee -a "$out_file"
echo "Full summary saved to: $out_file" | tee -a "$out_file"
echo "==========================================" | tee -a "$out_file"
