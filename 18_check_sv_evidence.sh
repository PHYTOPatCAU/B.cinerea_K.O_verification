#!/bin/bash
#SBATCH --job-name=check_sv
#SBATCH --output=check_sv.out
#SBATCH --error=check_sv.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00

module load gcc12-env
module load miniconda3
conda activate samtools

bam_dir="GenomesBGI/Bc/BAM/MD"
loci_file="GenomesBGI/Bc/GOI/loci.txt"
out_file="GenomesBGI/Bc/GOI/sv_summary.txt"

echo "=== Checking SV Evidence for Each Locus ===" | tee "$out_file"
echo "" | tee -a "$out_file"

samples=("18_x_KO" "23_x_KO_#7" "28xKO2" "29XKO1" "B_05.10" "WT" "WTB50")

while IFS=$'\t' read -r region name; do
    [[ "$region" == "Locus" ]] && continue
    [[ -z "$region" ]] && continue

    echo "=====================================" | tee -a "$out_file"
    echo "Locus: $name ($region)" | tee -a "$out_file"
    echo "=====================================" | tee -a "$out_file"

    for sample in "${samples[@]}"; do
        bam="${bam_dir}/${sample}.DuplMark.bam"
        [ ! -f "$bam" ] && continue

        echo "" | tee -a "$out_file"
        echo "Sample: $sample" | tee -a "$out_file"

        depth=$(samtools depth -r "$region" "$bam" | awk '{sum+=$3; n++} END {if(n>0) print sum/n; else print 0}')
        echo "  Average depth: $depth" | tee -a "$out_file"

        samtools view -f 1 "$bam" "$region" | \
        awk '{
            flag=$2;
            if(and(flag,16) && and(flag,32)) same_strand_rev++;
            else if(!and(flag,16) && !and(flag,32)) same_strand_fwd++;
            else if(and(flag,16) && !and(flag,32)) normal_fr++;
            else if(!and(flag,16) && and(flag,32)) normal_rf++;
            total++;
        } END {
            print "  Read pairs total:", total;
            print "    Normal (FR):", normal_fr+normal_rf;
            print "    Same strand (→→ or ←←):", same_strand_fwd+same_strand_rev, "(INVERSION evidence)";
        }' | tee -a "$out_file"

        samtools view -f 1 "$bam" "$region" | \
        awk '{
            isize=$9;
            if(isize<0) isize=-isize;
            if(isize > 1000) large++;
            else if(isize < 100 && isize > 0) small++;
            total++;
        } END {
            print "    Large insert (>1000bp):", large, "(DELETION evidence)";
            print "    Small insert (<100bp):", small, "(DUPLICATION evidence)";
        }' | tee -a "$out_file"

    done

    echo "" | tee -a "$out_file"

done < "$loci_file"

echo "" | tee -a "$out_file"
echo "=== Summary saved to: $out_file ===" | tee -a "$out_file"
echo "" | tee -a "$out_file"
echo "Interpretation:"
echo "  - Low depth (<10) = likely deletion"
echo "  - Same-strand pairs = INVERSION"
echo "  - Large insert sizes = DELETION"
echo "  - Small insert sizes = DUPLICATION"
echo "  - High depth (>50) = possible duplication"
