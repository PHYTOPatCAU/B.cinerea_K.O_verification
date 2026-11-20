#!/bin/bash
#SBATCH --job-name=SNPisec_hi
#SBATCH --output=SNPisec_hi.out
#SBATCH --error=SNPisec_hi.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

# Requires: 08_snpeff_annotation.sh output (snpeff.highimpact.vcf)
#           06_snp_sharing_analysis.sh output (per-sample biallelic VCFs in isec_results_all/)

workdir="GenomesBGI/Bc/VCF"
outdir="$workdir/isec_results_all"
highimpact_vcf="$workdir/Bc.SNP.INDEL.filtered.PASS.all.snpeff.highimpact.vcf"
highimpact_outdir="$workdir/isec_highimpact_results"
mkdir -p "$highimpact_outdir"

samples=("18_x_KO" "23_x_KO_#7" "28xKO2" "29XKO1" "B_05.10" "WT" "WTB50")
bg_samples=("B_05.10" "WT" "WTB50")
mutants=("18_x_KO" "23_x_KO_#7" "28xKO2" "29XKO1")

# 1. Make high-impact VCF biallelic, compress and index
vcftools --vcf "$highimpact_vcf" --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$highimpact_outdir/Bc.highimpact.biallelic"
bgzip -f -c "$highimpact_outdir/Bc.highimpact.biallelic.recode.vcf" > "$highimpact_outdir/Bc.highimpact.biallelic.vcf.gz"
tabix -f -p vcf "$highimpact_outdir/Bc.highimpact.biallelic.vcf.gz"

# 2. Extract per-sample VCFs
for sample in "${samples[@]}"; do
    bcftools view -c 1 -s "$sample" "$highimpact_outdir/Bc.highimpact.biallelic.vcf.gz" -Oz -o "$highimpact_outdir/${sample}.biallelic.vcf.gz"
    tabix -f -p vcf "$highimpact_outdir/${sample}.biallelic.vcf.gz"
done

# 3. Pairwise isec
for i in "${!samples[@]}"; do
    for j in "${!samples[@]}"; do
        if [ $i -lt $j ]; then
            s1="${samples[$i]}"
            s2="${samples[$j]}"
            bcftools isec "$highimpact_outdir/${s1}.biallelic.vcf.gz" "$highimpact_outdir/${s2}.biallelic.vcf.gz" -C -p "$highimpact_outdir/${s1}_vs_${s2}_unique"
            bcftools isec "$highimpact_outdir/${s2}.biallelic.vcf.gz" "$highimpact_outdir/${s1}.biallelic.vcf.gz" -C -p "$highimpact_outdir/${s2}_vs_${s1}_unique"
            bcftools isec "$highimpact_outdir/${s1}.biallelic.vcf.gz" "$highimpact_outdir/${s2}.biallelic.vcf.gz" -n=2 -p "$highimpact_outdir/${s1}_${s2}_shared"
        fi
    done
done

# 4. Background strain comparisons
for i in "${!bg_samples[@]}"; do
    for j in "${!bg_samples[@]}"; do
        if [ $i -lt $j ]; then
            s1="${bg_samples[$i]}"
            s2="${bg_samples[$j]}"
            bcftools isec "$highimpact_outdir/${s1}.biallelic.vcf.gz" "$highimpact_outdir/${s2}.biallelic.vcf.gz" -C -p "$highimpact_outdir/bg_${s1}_vs_${s2}_unique"
            bcftools isec "$highimpact_outdir/${s2}.biallelic.vcf.gz" "$highimpact_outdir/${s1}.biallelic.vcf.gz" -C -p "$highimpact_outdir/bg_${s2}_vs_${s1}_unique"
            bcftools isec "$highimpact_outdir/${s1}.biallelic.vcf.gz" "$highimpact_outdir/${s2}.biallelic.vcf.gz" -n=2 -p "$highimpact_outdir/bg_${s1}_${s2}_shared"
        fi
    done
done

# 5. Unique SNPs per sample
unique_counts_txt="$highimpact_outdir/unique_snp_counts.txt"
echo -e "Sample\tUnique_SNPs" > "$unique_counts_txt"
for sample in "${samples[@]}"; do
    count=0
    for other in "${samples[@]}"; do
        if [ "$sample" != "$other" ]; then
            vcf_dir="$highimpact_outdir/${sample}_vs_${other}_unique"
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
bg_similarity_txt="$highimpact_outdir/background_strain_similarity.txt"
echo -e "Comparison\tShared_SNPs\tUnique_SNPs_1\tUnique_SNPs_2" > "$bg_similarity_txt"
for i in "${!bg_samples[@]}"; do
    for j in "${!bg_samples[@]}"; do
        if [ $i -lt $j ]; then
            s1="${bg_samples[$i]}"
            s2="${bg_samples[$j]}"
            shared_dir="$highimpact_outdir/bg_${s1}_${s2}_shared"
            shared_count=0
            if [ -f "$shared_dir/0000.vcf" ]; then
                shared_count=$(grep -v "^#" "$shared_dir/0000.vcf" | wc -l)
            fi
            uniq1_dir="$highimpact_outdir/bg_${s1}_vs_${s2}_unique"
            uniq2_dir="$highimpact_outdir/bg_${s2}_vs_${s1}_unique"
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

# 7. B_05.10 vs mutants
b05_vs_mut_txt="$highimpact_outdir/B_05.10_vs_mutants_comparison.txt"
echo -e "Mutant\tShared_SNPs\tUnique_to_B_05.10\tUnique_to_Mutant" > "$b05_vs_mut_txt"
for mutant in "${mutants[@]}"; do
    shared_dir="$highimpact_outdir/B_05.10_${mutant}_shared"
    uniq_b05_dir="$highimpact_outdir/B_05.10_vs_${mutant}_unique"
    uniq_mut_dir="$highimpact_outdir/${mutant}_vs_B_05.10_unique"
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
mutant_comparison_txt="$highimpact_outdir/mutant_pairwise_comparison.txt"
echo -e "Comparison\tShared_SNPs\tUnique_SNPs_1\tUnique_SNPs_2" > "$mutant_comparison_txt"
for i in "${!mutants[@]}"; do
    for j in "${!mutants[@]}"; do
        if [ $i -lt $j ]; then
            m1="${mutants[$i]}"
            m2="${mutants[$j]}"
            shared_dir="$highimpact_outdir/${m1}_${m2}_shared"
            uniq1_dir="$highimpact_outdir/${m1}_vs_${m2}_unique"
            uniq2_dir="$highimpact_outdir/${m2}_vs_${m1}_unique"
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

echo "" >> "$summary_file"
echo "=== HIGH IMPACT VARIANTS PER SAMPLE ===" >> "$summary_file"
echo -e "Sample\tHigh_Impact_Variants" >> "$summary_file"
for sample in "${samples[@]}"; do
    count=0
    [ -f "$highimpact_outdir/${sample}.biallelic.vcf.gz" ] && count=$(zcat "$highimpact_outdir/${sample}.biallelic.vcf.gz" 2>/dev/null | grep -v "^#" | wc -l)
    echo -e "${sample}\t${count}" >> "$summary_file"
done

echo "" >> "$summary_file"
echo "=== HIGH IMPACT MUTANT-SPECIFIC VARIANTS (vs ALL WT) ===" >> "$summary_file"
echo -e "Mutant\tUnique_HighImpact_vs_WT" >> "$summary_file"

bcftools merge "$highimpact_outdir/B_05.10.biallelic.vcf.gz" \
    "$highimpact_outdir/WT.biallelic.vcf.gz" \
    "$highimpact_outdir/WTB50.biallelic.vcf.gz" \
    -Oz -o "$highimpact_outdir/all_WT_merged.vcf.gz" 2>/dev/null
tabix -f -p vcf "$highimpact_outdir/all_WT_merged.vcf.gz"

bcftools query -f '%CHROM\t%POS\n' "$highimpact_outdir/all_WT_merged.vcf.gz" | sort -u > "$highimpact_outdir/all_WT_positions.txt"

for mutant in "${mutants[@]}"; do
    bcftools query -f '%CHROM\t%POS\n' "$highimpact_outdir/${mutant}.biallelic.vcf.gz" 2>/dev/null | sort -u > "$highimpact_outdir/${mutant}_positions.txt"
    comm -23 "$highimpact_outdir/${mutant}_positions.txt" "$highimpact_outdir/all_WT_positions.txt" > "$highimpact_outdir/${mutant}_unique_positions.txt"
    unique_count=$(wc -l < "$highimpact_outdir/${mutant}_unique_positions.txt")
    echo -e "${mutant}\t${unique_count}" >> "$summary_file"
done

echo "" >> "$summary_file"
echo "========================================================" >> "$summary_file"
echo "High-impact summary appended to: $summary_file" >> "$summary_file"
echo "========================================================" >> "$summary_file"

echo ""
echo "✓ High-impact isec complete!"
echo "✓ Results: $highimpact_outdir"
