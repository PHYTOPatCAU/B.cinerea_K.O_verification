#!/bin/bash
#SBATCH --job-name=snpeff_annotate
#SBATCH --output=snpeff_annotate.out
#SBATCH --error=snpeff_annotate.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=08:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

snpeff_dir="snpEff"
snpeff_db="Botrytis_cinerea"
vcf_dir="GenomesBGI/Bc/VCF"
input_vcf="$vcf_dir/Bc.SNP.INDEL.filtered.PASS.all.vcf"
annot_vcf="$vcf_dir/Bc.SNP.INDEL.filtered.PASS.all.snpeff.vcf"
highimpact_vcf="$vcf_dir/Bc.SNP.INDEL.filtered.PASS.all.snpeff.highimpact.vcf"
highimpact_txt="$vcf_dir/Bc.SNP.INDEL.filtered.PASS.all.snpeff.highimpact.txt"

cd "$snpeff_dir"
java -mx10G -jar $snpeff_dir/snpEff.jar -v $snpeff_db "$input_vcf" > "$annot_vcf"
java -mx10G -jar $snpeff_dir/snpEff.jar -stats "$vcf_dir/Bc.SNP.INDEL.filtered.PASS.all.snpeff.stats.html" "$snpeff_db" "$input_vcf"

cat "$annot_vcf" | \
    java -jar $snpeff_dir/SnpSift.jar filter "(ANN[0].IMPACT = 'HIGH')" > "$highimpact_vcf"

highimpact_count=$(grep -v "^#" "$highimpact_vcf" | wc -l)
echo "High impact variant count: $highimpact_count" > "$highimpact_txt"

cat "$highimpact_vcf" | \
    perl $snpeff_dir/scripts/vcfEffOnePerLine.pl | \
    java -jar $snpeff_dir/SnpSift.jar extractFields - CHROM POS REF ALT AF "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" >> "$highimpact_txt"
