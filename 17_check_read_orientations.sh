#!/bin/bash

module load gcc12-env
module load miniconda3
conda activate samtools

region="5:1329646-1331550"
bam_dir="GenomesBGI/Bc/BAM/MD"

echo "Checking read pair orientations in deletion region: $region"
echo ""

for sample in 18_x_KO "23_x_KO_#7" 28xKO2 29XKO1 B_05.10 WT WTB50; do
    bam="${bam_dir}/${sample}.DuplMark.bam"

    if [ ! -f "$bam" ]; then
        echo "Sample: $sample"
        echo "----------------------------------------"
        echo "  ⚠️  BAM file not found: $bam"
        echo ""
        continue
    fi

    echo "Sample: $sample"
    echo "----------------------------------------"

    samtools view -f 1 "$bam" "$region" 2>/dev/null | awk '
    BEGIN {
        normal_fr = 0;
        normal_rf = 0;
        same_fwd = 0;
        same_rev = 0;
        large_insert = 0;
        total = 0;
    }
    {
        flag = $2;
        isize = $9;
        if (isize < 0) isize = -isize;

        read_reverse = and(flag, 16);
        mate_reverse = and(flag, 32);

        if (read_reverse && mate_reverse) {
            same_rev++;
        } else if (!read_reverse && !mate_reverse) {
            same_fwd++;
        } else {
            if (isize > 1000) large_insert++;
            if (read_reverse) normal_rf++;
            else normal_fr++;
        }
        total++;
    }
    END {
        print "  Total read pairs:", total;
        print "  Normal orientation (F→R or R→F):", normal_fr + normal_rf;
        print "    - Large insert (>1kb, DELETION):", large_insert;
        print "  INVERSION evidence (same strand):";
        print "    - Both forward (→→):", same_fwd;
        print "    - Both reverse (←←):", same_rev;
        print "";

        if (same_fwd + same_rev > total * 0.1) {
            print "  ⚠️  INVERSION DETECTED (>10% same-strand pairs)";
        } else if (large_insert > total * 0.3) {
            print "  ✓ DELETION CONFIRMED (>30% large inserts)";
        } else if (total < 10) {
            print "  ✓ DELETION CONFIRMED (no/low coverage)";
        } else {
            print "  ✓ NORMAL (wildtype-like)";
        }
        print "";
    }
    '
done

echo ""
echo "Legend:"
echo "  →→ or ←← = Inversion (both reads same direction)"
echo "  →← or ←→ = Normal (reads face each other)"
echo "  Large insert = Deletion (reads far apart)"
