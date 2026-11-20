#!/bin/bash
#SBATCH --job-name=SNPisec
#SBATCH --output=SNPisec.out
#SBATCH --error=SNPisec.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

vcf="GenomesBGI/Bc/VCF/Bc.SNP.INDEL.filtered.PASS.all.vcf"
workdir="GenomesBGI/Bc/VCF"
outdir="$workdir/isec_results_all"
mkdir -p "$outdir"

samples=("18_x_KO" "23_x_KO_#7" "28xKO2" "29XKO1" "B_05.10" "WT" "WTB50")
bg_samples=("B_05.10" "WT" "WTB50")
mutants=("18_x_KO" "23_x_KO_#7" "28xKO2" "29XKO1")

# 1. Make the whole VCF biallelic, then compress and index
vcftools --vcf "$vcf" --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$outdir/Bc.biallelic"
bgzip -f -c "$outdir/Bc.biallelic.recode.vcf" > "$outdir/Bc.biallelic.vcf.gz"
tabix -f -p vcf "$outdir/Bc.biallelic.vcf.gz"

# 2. Extract per-sample VCFs
for sample in "${samples[@]}"; do
    bcftools view -c 1 -s "$sample" "$outdir/Bc.biallelic.vcf.gz" -Oz -o "$outdir/${sample}.biallelic.vcf.gz"
    tabix -f -p vcf "$outdir/${sample}.biallelic.vcf.gz"
done

# 3. Systematic pairwise bcftools isec for shared and unique SNPs/INDELs
for i in "${!samples[@]}"; do
    for j in "${!samples[@]}"; do
        if [ $i -lt $j ]; then
            s1="${samples[$i]}"
            s2="${samples[$j]}"
            bcftools isec "$outdir/${s1}.biallelic.vcf.gz" "$outdir/${s2}.biallelic.vcf.gz" -C -p "$outdir/${s1}_vs_${s2}_unique"
            bcftools isec "$outdir/${s2}.biallelic.vcf.gz" "$outdir/${s1}.biallelic.vcf.gz" -C -p "$outdir/${s2}_vs_${s1}_unique"
            bcftools isec "$outdir/${s1}.biallelic.vcf.gz" "$outdir/${s2}.biallelic.vcf.gz" -n=2 -p "$outdir/${s1}_${s2}_shared"
        fi
    done
done

# 4. Background strain comparisons
for i in "${!bg_samples[@]}"; do
    for j in "${!bg_samples[@]}"; do
        if [ $i -lt $j ]; then
            s1="${bg_samples[$i]}"
            s2="${bg_samples[$j]}"
            bcftools isec "$outdir/${s1}.biallelic.vcf.gz" "$outdir/${s2}.biallelic.vcf.gz" -C -p "$outdir/bg_${s1}_vs_${s2}_unique"
            bcftools isec "$outdir/${s2}.biallelic.vcf.gz" "$outdir/${s1}.biallelic.vcf.gz" -C -p "$outdir/bg_${s2}_vs_${s1}_unique"
            bcftools isec "$outdir/${s1}.biallelic.vcf.gz" "$outdir/${s2}.biallelic.vcf.gz" -n=2 -p "$outdir/bg_${s1}_${s2}_shared"
        fi
    done
done

# 5. Count unique SNPs per sample
unique_counts_txt="$outdir/unique_snp_counts.txt"
echo -e "Sample\tUnique_SNPs" > "$unique_counts_txt"
for sample in "${samples[@]}"; do
    count=0
    for other in "${samples[@]}"; do
        if [ "$sample" != "$other" ]; then
            vcf_dir="$outdir/${sample}_vs_${other}_unique"
            if [ -f "$vcf_dir/0000.vcf" ]; then
                snps=$(grep -v "^#" "$vcf_dir/0000.vcf" | wc -l)
                count=$((count + snps))
            else
                snps=0
            fi
        fi
    done
    echo -e "${sample}\t${count}" >> "$unique_counts_txt"
done

# 6. Background strain similarity
bg_similarity_txt="$outdir/background_strain_similarity.txt"
echo -e "Comparison\tShared_SNPs\tUnique_SNPs_1\tUnique_SNPs_2" > "$bg_similarity_txt"
for i in "${!bg_samples[@]}"; do
    for j in "${!bg_samples[@]}"; do
        if [ $i -lt $j ]; then
            s1="${bg_samples[$i]}"
            s2="${bg_samples[$j]}"
            shared_dir="$outdir/bg_${s1}_${s2}_shared"
            shared_count=0
            if [ -f "$shared_dir/0000.vcf" ]; then
                shared_count=$(grep -v "^#" "$shared_dir/0000.vcf" | wc -l)
            fi
            uniq1_dir="$outdir/bg_${s1}_vs_${s2}_unique"
            uniq2_dir="$outdir/bg_${s2}_vs_${s1}_unique"
            uniq1=0
            uniq2=0
            if [ -f "$uniq1_dir/0000.vcf" ]; then
                uniq1=$(grep -v "^#" "$uniq1_dir/0000.vcf" | wc -l)
            fi
            if [ -f "$uniq2_dir/0000.vcf" ]; then
                uniq2=$(grep -v "^#" "$uniq2_dir/0000.vcf" | wc -l)
            fi
            echo -e "${s1} vs ${s2}\t${shared_count}\t${uniq1}\t${uniq2}" >> "$bg_similarity_txt"
        fi
    done
done

# 7. B_05.10 vs all mutants
b05_vs_mut_txt="$outdir/B_05.10_vs_mutants_comparison.txt"
echo -e "Mutant\tShared_SNPs\tUnique_to_B_05.10\tUnique_to_Mutant" > "$b05_vs_mut_txt"
for mutant in "${mutants[@]}"; do
    shared_dir="$outdir/B_05.10_${mutant}_shared"
    uniq_b05_dir="$outdir/B_05.10_vs_${mutant}_unique"
    uniq_mut_dir="$outdir/${mutant}_vs_B_05.10_unique"
    shared=0
    uniq_b05=0
    uniq_mut=0
    if [ -f "$shared_dir/0000.vcf" ]; then
        shared=$(grep -v "^#" "$shared_dir/0000.vcf" | wc -l)
    fi
    if [ -f "$uniq_b05_dir/0000.vcf" ]; then
        uniq_b05=$(grep -v "^#" "$uniq_b05_dir/0000.vcf" | wc -l)
    fi
    if [ -f "$uniq_mut_dir/0000.vcf" ]; then
        uniq_mut=$(grep -v "^#" "$uniq_mut_dir/0000.vcf" | wc -l)
    fi
    echo -e "${mutant}\t${shared}\t${uniq_b05}\t${uniq_mut}" >> "$b05_vs_mut_txt"
done

# 8. Pairwise mutant comparison
mutant_comparison_txt="$outdir/mutant_pairwise_comparison.txt"
echo -e "Comparison\tShared_SNPs\tUnique_SNPs_1\tUnique_SNPs_2" > "$mutant_comparison_txt"
for i in "${!mutants[@]}"; do
    for j in "${!mutants[@]}"; do
        if [ $i -lt $j ]; then
            m1="${mutants[$i]}"
            m2="${mutants[$j]}"
            shared_dir="$outdir/${m1}_${m2}_shared"
            uniq1_dir="$outdir/${m1}_vs_${m2}_unique"
            uniq2_dir="$outdir/${m2}_vs_${m1}_unique"
            shared_count=0
            uniq1=0
            uniq2=0
            if [ -f "$shared_dir/0000.vcf" ]; then
                shared_count=$(grep -v "^#" "$shared_dir/0000.vcf" | wc -l)
            fi
            if [ -f "$uniq1_dir/0000.vcf" ]; then
                uniq1=$(grep -v "^#" "$uniq1_dir/0000.vcf" | wc -l)
            fi
            if [ -f "$uniq2_dir/0000.vcf" ]; then
                uniq2=$(grep -v "^#" "$uniq2_dir/0000.vcf" | wc -l)
            fi
            echo -e "${m1} vs ${m2}\t${shared_count}\t${uniq1}\t${uniq2}" >> "$mutant_comparison_txt"
        fi
    done
done

summary_file="$workdir/SNP_sharing_final_summary.txt"
echo "========================================================" > "$summary_file"
echo "       SNP SHARING ANALYSIS - FINAL SUMMARY             " >> "$summary_file"
echo "       Generated: $(date)                               " >> "$summary_file"
echo "========================================================" >> "$summary_file"

echo "" >> "$summary_file"
echo "=== TOTAL VARIANTS PER SAMPLE (All SNPs/INDELs) ===" >> "$summary_file"
echo -e "Sample\tTotal_Variants" >> "$summary_file"
for sample in "${samples[@]}"; do
    count=$(zcat "$outdir/${sample}.biallelic.vcf.gz" 2>/dev/null | grep -v "^#" | wc -l)
    echo -e "${sample}\t${count}" >> "$summary_file"
done

echo "" >> "$summary_file"
echo "=== MUTANT-SPECIFIC VARIANTS (Unique vs ALL WT strains pooled) ===" >> "$summary_file"
echo -e "Mutant\tUnique_vs_WT_pooled" >> "$summary_file"

bcftools merge "$outdir/B_05.10.biallelic.vcf.gz" \
    "$outdir/WT.biallelic.vcf.gz" \
    "$outdir/WTB50.biallelic.vcf.gz" \
    -Oz -o "$outdir/all_WT_merged.vcf.gz" 2>/dev/null
tabix -f -p vcf "$outdir/all_WT_merged.vcf.gz"

bcftools query -f '%CHROM\t%POS\n' "$outdir/all_WT_merged.vcf.gz" | sort -u > "$outdir/all_WT_positions.txt"

for mutant in "${mutants[@]}"; do
    bcftools query -f '%CHROM\t%POS\n' "$outdir/${mutant}.biallelic.vcf.gz" | sort -u > "$outdir/${mutant}_positions.txt"
    comm -23 "$outdir/${mutant}_positions.txt" "$outdir/all_WT_positions.txt" > "$outdir/${mutant}_unique_positions.txt"
    unique_count=$(wc -l < "$outdir/${mutant}_unique_positions.txt")
    echo -e "${mutant}\t${unique_count}" >> "$summary_file"
    if [ "$unique_count" -gt 0 ]; then
        bcftools view -T "$outdir/${mutant}_unique_positions.txt" \
            "$outdir/${mutant}.biallelic.vcf.gz" \
            -Oz -o "$outdir/${mutant}_unique_vs_allWT.vcf.gz"
    fi
done

echo "" >> "$summary_file"
echo "=== PAIRWISE MUTANT COMPARISON ===" >> "$summary_file"
echo -e "Mutant1\tMutant2\tShared\tUnique_to_M1\tUnique_to_M2" >> "$summary_file"
for i in "${!mutants[@]}"; do
    for j in "${!mutants[@]}"; do
        if [ $i -lt $j ]; then
            m1="${mutants[$i]}"
            m2="${mutants[$j]}"
            shared=0
            uniq1=0
            uniq2=0
            [ -f "$outdir/${m1}_${m2}_shared/0000.vcf" ] && shared=$(grep -v "^#" "$outdir/${m1}_${m2}_shared/0000.vcf" | wc -l)
            [ -f "$outdir/${m1}_vs_${m2}_unique/0000.vcf" ] && uniq1=$(grep -v "^#" "$outdir/${m1}_vs_${m2}_unique/0000.vcf" | wc -l)
            [ -f "$outdir/${m2}_vs_${m1}_unique/0000.vcf" ] && uniq2=$(grep -v "^#" "$outdir/${m2}_vs_${m1}_unique/0000.vcf" | wc -l)
            echo -e "${m1}\t${m2}\t${shared}\t${uniq1}\t${uniq2}" >> "$summary_file"
        fi
    done
done

echo "" >> "$summary_file"
echo "=== BACKGROUND STRAIN SIMILARITY ===" >> "$summary_file"
echo -e "Strain1\tStrain2\tShared\tUnique_to_S1\tUnique_to_S2" >> "$summary_file"
for i in "${!bg_samples[@]}"; do
    for j in "${!bg_samples[@]}"; do
        if [ $i -lt $j ]; then
            s1="${bg_samples[$i]}"
            s2="${bg_samples[$j]}"
            shared=0
            uniq1=0
            uniq2=0
            [ -f "$outdir/bg_${s1}_${s2}_shared/0000.vcf" ] && shared=$(grep -v "^#" "$outdir/bg_${s1}_${s2}_shared/0000.vcf" | wc -l)
            [ -f "$outdir/bg_${s1}_vs_${s2}_unique/0000.vcf" ] && uniq1=$(grep -v "^#" "$outdir/bg_${s1}_vs_${s2}_unique/0000.vcf" | wc -l)
            [ -f "$outdir/bg_${s2}_vs_${s1}_unique/0000.vcf" ] && uniq2=$(grep -v "^#" "$outdir/bg_${s2}_vs_${s1}_unique/0000.vcf" | wc -l)
            echo -e "${s1}\t${s2}\t${shared}\t${uniq1}\t${uniq2}" >> "$summary_file"
        fi
    done
done

echo "" >> "$summary_file"
echo "========================================================" >> "$summary_file"
echo "Summary saved to: $summary_file" >> "$summary_file"
echo "========================================================" >> "$summary_file"

echo ""
echo "✓ Analysis complete!"
echo "✓ Final summary: $summary_file"
cat "$summary_file"
