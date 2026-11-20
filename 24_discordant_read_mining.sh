#!/bin/bash
#SBATCH --job-name=discordant_reads
#SBATCH --output=discordant_reads_%j.out
#SBATCH --error=discordant_reads_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL,END

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

region="5:1329646-1331550"
chr="5"
del_start=1329646
del_end=1331550

wide_start=1319646
wide_end=1341550

bam_dir="GenomesBGI/Bc/BAM/MD"
out_dir="GenomesBGI/Bc/GOI/discordant_investigation"
mkdir -p "$out_dir"

echo "=========================================="
echo "DISCORDANT READ PAIR MINING"
echo "Looking for reads whose mate maps outside deletion"
echo "Deletion: $region"
echo "Wide search window: ${chr}:${wide_start}-${wide_end}"
echo "=========================================="
echo ""

for sample in 18_x_KO "23_x_KO_#7" 28xKO2 29XKO1; do
    bam="${bam_dir}/${sample}.DuplMark.bam"
    [ ! -f "$bam" ] && echo "WARNING: $bam not found, skipping" && continue

    echo "=== Sample: $sample ==="
    sample_dir="${out_dir}/${sample}"
    mkdir -p "$sample_dir"

    echo "[1] Extracting read names from deletion region..."
    readnames_file="${sample_dir}/deletion_readnames.txt"
    samtools view -F 256 "$bam" "$region" | awk '{print $1}' | sort -u > "$readnames_file"
    n_reads=$(wc -l < "$readnames_file")
    echo "    Found $n_reads unique read names in deletion region"

    if [ "$n_reads" -eq 0 ]; then
        echo "    No reads found in deletion region - deletion likely complete"
        echo ""
        continue
    fi

    echo "[2] Fetching full read pairs from wide window (${chr}:${wide_start}-${wide_end})..."
    pairs_bam="${sample_dir}/discordant_pairs.bam"
    samtools view -b -N "$readnames_file" "$bam" "${chr}:${wide_start}-${wide_end}" \
        | samtools sort -o "$pairs_bam"
    samtools index "$pairs_bam"
    n_pairs=$(samtools view -c "$pairs_bam")
    echo "    Fetched $n_pairs reads in pair set"

    echo "[3] Identifying mate positions outside deletion region..."
    mate_positions="${sample_dir}/mate_positions.txt"
    samtools view "$pairs_bam" | awk \
        -v del_start="$del_start" \
        -v del_end="$del_end" \
        -v chr="$chr" '
    {
        flag = $2;
        rnext = $7;
        pnext = $8;
        isize = $9;
        if (isize < 0) isize = -isize;

        if (rnext == "=" || rnext == chr) {
            if (pnext < del_start - 200 || pnext > del_end + 200) {
                print chr, pnext, isize, $1;
            }
        }
    }' | sort -k2,2n > "$mate_positions"

    n_discordant=$(wc -l < "$mate_positions")
    echo "    Found $n_discordant discordant mate positions"

    echo "[4] Clustering mate positions (candidate insertion sites)..."
    cluster_file="${sample_dir}/candidate_insertion_sites.txt"
    echo -e "Chr\tPosition\tCount\tMean_InsertSize" > "$cluster_file"

    awk '
    BEGIN { window = 200; prev = -9999; count = 0; sum_isize = 0; cluster_pos = 0; }
    {
        chr = $1; pos = $2; isize = $3;
        if (pos - prev <= window) {
            count++;
            sum_isize += isize;
            cluster_pos = int((cluster_pos * (count-1) + pos) / count);
        } else {
            if (count >= 2) {
                printf "%s\t%d\t%d\t%.0f\n", chr, cluster_pos, count, sum_isize/count;
            }
            count = 1;
            sum_isize = isize;
            cluster_pos = pos;
            prev = pos;
        }
        prev = pos;
    }
    END {
        if (count >= 2) printf "%s\t%d\t%d\t%.0f\n", chr, cluster_pos, count, sum_isize/count;
    }
    ' "$mate_positions" >> "$cluster_file"

    echo ""
    echo "    === CANDIDATE REINSERTION SITES ==="
    if [ $(wc -l < "$cluster_file") -le 1 ]; then
        echo "    No clusters found (need >=2 mates at same location)"
    else
        cat "$cluster_file"
    fi

    echo ""
    echo "[5] Top mate landing positions (raw counts):"
    awk '{print $1":"$2}' "$mate_positions" | \
        awk -F: '{bucket=int($2/500)*500; print $1":"bucket"-"bucket+500}' | \
        sort | uniq -c | sort -rn | head -10 | \
        awk '{printf "    %4d reads → %s\n", $1, $2}'

    echo ""
    echo "    Full results: $sample_dir/"
    echo "----------------------------------------"
    echo ""
done

echo ""
echo "=========================================="
echo "INTERPRETATION"
echo "=========================================="
echo ""
echo "Look for:"
echo "  - Clusters of >=3 mates at same position → candidate insertion site"
echo "  - Consistent position across multiple KO samples → high confidence"
echo "  - If mates cluster ~2kb from deletion → confirms reinsertion hypothesis"
echo ""
echo "  → Check candidate regions in IGV with discordant_pairs.bam loaded"
echo "=========================================="
