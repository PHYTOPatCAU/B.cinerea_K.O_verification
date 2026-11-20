#!/bin/bash
#SBATCH --job-name=automated_pipeline_Bc
#SBATCH --qos=long
#SBATCH --out=automated_pipeline_Bc_%j.out
#SBATCH --error=automated_pipeline_Bc_%j.err
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

organism="Bc"
ref_file="Botrytis_cinerea.ASM83294v1.dna.toplevel.fa"
batch_dir="GenomesBGI/Bc"
input_dir="$batch_dir/fastq"
output_dir="$batch_dir/BAM"
ref="$batch_dir/$ref_file"

#mkdir -p "$output_dir" "$output_dir/RG" "$output_dir/MD" \
#         "$batch_dir/gVCF" "$batch_dir/VCF" \
#         "$batch_dir/covstats"

# 1. BWA MEM
bwa_jobids=()
echo "=== Checking for samples needing alignment ==="
for file in "$input_dir"/*_1.fq.gz; do
    sample=$(basename "$file" _1.fq.gz)
    if [ ! -f "$input_dir/${sample}_2.fq.gz" ]; then continue; fi
    if [ -f "$output_dir/${sample}.bam" ]; then
        echo "  - $sample: BAM exists, skipping"
        continue
    fi
    echo "  - $sample: Submitting BWA alignment"
    jid=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=bwa_${sample}
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=${output_dir}/bwa_${sample}.out
#SBATCH --error=${output_dir}/bwa_${sample}.err
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools
bwa mem -t 16 "$ref" "$input_dir/${sample}_1.fq.gz" "$input_dir/${sample}_2.fq.gz" | \
    samtools sort -@ 8 -o "$output_dir/${sample}.bam" -
samtools index "$output_dir/${sample}.bam"
EOF
)
    bwa_jobids+=($(echo $jid | awk '{print $4}'))
done
if [ ${#bwa_jobids[@]} -gt 0 ]; then
    dependency_arg="--dependency=afterok:$(IFS=:; echo "${bwa_jobids[*]}")"
else
    dependency_arg=""
fi

# 2. AddOrReplaceReadGroups
rg_files_exist=true
for bam in "$output_dir"/*.bam; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" .bam)
    if [ ! -f "$output_dir/RG/${sample}.bam" ]; then
        rg_files_exist=false
        break
    fi
done

# If new BWA jobs were submitted OR RG files don't exist, submit RG job
if [ "$rg_files_exist" = false ] || [ ${#bwa_jobids[@]} -gt 0 ]; then
    jid2=$(sbatch $dependency_arg <<EOF
#!/bin/bash
#SBATCH --job-name=AddRG_${organism}
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=${batch_dir}/AddRG_${organism}.out
#SBATCH --error=${batch_dir}/AddRG_${organism}.err
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call
for bam in "$output_dir"/*.bam; do
    sample=\$(basename "\$bam" .bam)
    if [ ! -f "$output_dir/RG/\$sample.bam" ]; then
        echo "Processing \$sample"
        gatk AddOrReplaceReadGroups \
            --INPUT "\$bam" \
            --OUTPUT "$output_dir/RG/\$sample.bam" \
            --SORT_ORDER coordinate \
            --RGID \$sample \
            --RGLB Bc_WGS \
            --RGPL illumina \
            --RGSM \$sample \
            --RGPU BcBatch \
            --CREATE_INDEX true
    fi
done
EOF
)
    jid2=$(echo $jid2 | awk '{print $4}')
    rg_dependency_arg="--dependency=afterok:$jid2"
else
    rg_dependency_arg=""
fi

# 3. MarkDuplicates
md_files_exist=true
for bam in "$output_dir/RG"/*.bam; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" .bam)
    if [ ! -f "$output_dir/MD/${sample}.DuplMark.bam" ]; then
        md_files_exist=false
        break
    fi
done

# If RG job was submitted OR MD files don't exist, submit MD job
if [ "$md_files_exist" = false ] || [ -n "$rg_dependency_arg" ]; then
    jid3=$(sbatch $rg_dependency_arg <<EOF
#!/bin/bash
#SBATCH --job-name=MarkDup_${organism}
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=${batch_dir}/markdup_${organism}.out
#SBATCH --error=${batch_dir}/markdup_${organism}.err
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call
for bam in "$output_dir/RG"/*.bam; do
    sample=\$(basename "\$bam" .bam)
    if [ ! -f "$output_dir/MD/\$sample.DuplMark.bam" ]; then
        echo "Processing \$sample"
        gatk MarkDuplicates \
            --INPUT "\$bam" \
            --OUTPUT "$output_dir/MD/\$sample.DuplMark.bam" \
            --METRICS_FILE "$output_dir/MD/\$sample.metrics.txt" \
            --CREATE_INDEX true
    fi
done
EOF
)
    jid3=$(echo $jid3 | awk '{print $4}')
    md_dependency_arg="--dependency=afterok:$jid3"
else
    md_dependency_arg=""
fi

# --- Fix: If MarkDuplicates already ran, do not set dependency for downstream steps ---
if [ "$md_files_exist" = true ]; then
    downstream_dependency=""
else
    downstream_dependency="$md_dependency_arg"
fi

# 4. Coverage summary (optional - for monitoring only)
summary_file="$batch_dir/covstats/summary_mapping_depth.txt"
mkdir -p "$batch_dir/covstats"
echo -e "Sample\tAverageDepth\tMappingRate" > "$summary_file"

# Calculate coverage stats for monitoring (not for filtering)
for bam in "$output_dir/MD/"*.DuplMark.bam; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" .DuplMark.bam)

    avg_depth=$(samtools depth -a "$bam" | awk '{sum+=$3; n++} END{if(n>0) printf "%.2f", sum/n; else print 0}')
    total_reads=$(samtools view -c "$bam")
    mapped_reads=$(samtools view -c -F 4 "$bam")
    if [ "$total_reads" -gt 0 ]; then
        mapping_rate=$(awk "BEGIN {printf \"%.4f\", $mapped_reads/$total_reads}")
    else
        mapping_rate=0
    fi

    echo -e "${sample}\t${avg_depth}\t${mapping_rate}" >> "$summary_file"
done

# 5. HaplotypeCaller (per-sample parallel) - process ALL MD BAMs
hapcall_jobids=()
new_hapcall_needed=false
for bam in "$output_dir/MD/"*.DuplMark.bam; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" .DuplMark.bam)
    out="$batch_dir/gVCF/${sample}.g.vcf"

    if [ ! -f "$out" ]; then
        new_hapcall_needed=true
        jid=$(sbatch $downstream_dependency <<EOF
#!/bin/bash
#SBATCH --job-name=HapCall_${sample}
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=${batch_dir}/HapCall_${sample}.out
#SBATCH --error=${batch_dir}/HapCall_${sample}.err
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call
gatk --java-options "-Xmx14G" HaplotypeCaller \
    -R "$ref" -I "$bam" -O "$out" -ERC GVCF --sample-ploidy 1
EOF
)
        hapcall_jobids+=($(echo $jid | awk '{print $4}'))
    fi
done

# If new HaplotypeCaller jobs were submitted, remove old merged gVCF to force regeneration
if [ "$new_hapcall_needed" = true ]; then
    merged_gvcf="$batch_dir/gVCF/Bc_current.g.vcf"
    if [ -f "$merged_gvcf" ]; then
        echo "Removing old merged gVCF - new samples need to be included"
        timestamp=$(date +%Y%m%d_%H%M%S)
        mv "$merged_gvcf" "${merged_gvcf}.backup_$timestamp"
        [ -f "${merged_gvcf}.idx" ] && mv "${merged_gvcf}.idx" "${merged_gvcf}.idx.backup_$timestamp"
    fi
fi

if [ ${#hapcall_jobids[@]} -gt 0 ]; then
    hapcall_dependency="--dependency=afterok:$(IFS=:; echo "${hapcall_jobids[*]}")"
else
    hapcall_dependency=""
fi

# 6. Merge all gVCFs into Bc_current.g.vcf using GATK CombineGVCFs
merged_gvcf="$batch_dir/gVCF/Bc_current.g.vcf"
if [ ! -f "$merged_gvcf" ]; then
    gvcf_files=()
    for bam in "$output_dir/MD/"*.DuplMark.bam; do
        [ -f "$bam" ] || continue
        sample=$(basename "$bam" .DuplMark.bam)
        gvcf="$batch_dir/gVCF/${sample}.g.vcf"
        [ -f "$gvcf" ] && gvcf_files+=("$gvcf")
    done

    if [ ${#gvcf_files[@]} -eq 0 ]; then
        echo "No gVCF files found. Exiting."
        exit 1
    fi

    combine_args=""
    for f in "${gvcf_files[@]}"; do
        combine_args="$combine_args -V $f"
    done

    if [ ${#gvcf_files[@]} -eq 1 ]; then
        cp "${gvcf_files[0]}" "$merged_gvcf"
        combine_dependency=""
    else
        combinejid=$(sbatch $hapcall_dependency <<EOF
#!/bin/bash
#SBATCH --job-name=CombineGVCFs_Bc
#SBATCH --output=${batch_dir}/combine_Bc.out
#SBATCH --error=${batch_dir}/combine_Bc.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call
gatk CombineGVCFs $combine_args -O "$merged_gvcf" -R "$ref"
EOF
)
        combinejid=$(echo $combinejid | awk '{print $4}')
        combine_dependency="--dependency=afterok:$combinejid"
    fi
else
    combine_dependency=""
fi

# 7. GenotypeGVCFs on merged Bc_current.g.vcf (only one job, depends on combine)
genotyped_vcf="$batch_dir/VCF/Bc_genotyped.all.vcf"
mkdir -p "$batch_dir/VCF"
geno_jobid=""
if [ ! -f "$genotyped_vcf" ]; then
    geno_jobid=$(sbatch $combine_dependency <<EOF
#!/bin/bash
#SBATCH --job-name=GenotypeGVCFs_Bc
#SBATCH --output=${batch_dir}/geno_Bc.out
#SBATCH --error=${batch_dir}/geno_Bc.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call
gatk --java-options "-Xmx14G" GenotypeGVCFs -R "$ref" -V "$merged_gvcf" -O "$genotyped_vcf"
EOF
)
    geno_jobid=$(echo $geno_jobid | awk '{print $4}')
    geno_dependency="--dependency=afterok:$geno_jobid"
else
    geno_dependency=""
fi

# 8. Select SNPs and INDELs from merged genotyped VCF (depends on GenotypeGVCFs)
selected_vcf="$batch_dir/VCF/Bc.SNP.INDEL.all.vcf"
sel_jobid=""
if [ ! -f "$selected_vcf" ]; then
    sel_jobid=$(sbatch $geno_dependency <<EOF
#!/bin/bash
#SBATCH --job-name=SelectVariants_Bc
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=${batch_dir}/Select_Bc.out
#SBATCH --error=${batch_dir}/Select_Bc.err
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call
gatk SelectVariants \
    -R "$ref" -V "$genotyped_vcf" \
    --select-type-to-include SNP \
    --select-type-to-include INDEL \
    -O "$selected_vcf"
EOF
)
    sel_jobid=$(echo $sel_jobid | awk '{print $4}')
    sel_dependency="--dependency=afterok:$sel_jobid"
else
    sel_dependency=""
fi

# 9. Quality filtering on merged selected VCF (depends on SelectVariants)
filtered_vcf="$batch_dir/VCF/Bc.SNP.INDEL.filtered.all.vcf"
filt_jobid=""
if [ ! -f "$filtered_vcf" ]; then
    filt_jobid=$(sbatch $sel_dependency <<EOF
#!/bin/bash
#SBATCH --job-name=FilterVariants_Bc
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=${batch_dir}/Filter_Bc.out
#SBATCH --error=${batch_dir}/Filter_Bc.err
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call
gatk --java-options "-Xmx6G" VariantFiltration \
    -R "$ref" -V "$selected_vcf" \
    --filter-name "DP" --filter-expression "DP < 10" \
    --filter-name "QUAL20" --filter-expression "QUAL < 20.0" \
    --filter-name "SNP_near_INDEL" --filter-expression "vc.hasNearbyIndel(5)" \
    --genotype-filter-name "lowGQ" --genotype-filter-expression "GQ < 20" \
    --genotype-filter-name "lowDP" --genotype-filter-expression "DP < 10" \
    --set-filtered-genotype-to-no-call true \
    -O "$filtered_vcf"
EOF
)
    filt_jobid=$(echo $filt_jobid | awk '{print $4}')
    filt_dependency="--dependency=afterok:$filt_jobid"
else
    filt_dependency=""
fi

# 10. Exclude filtered variants (keep only PASS from merged filtered VCF, depends on FilterVariants)
excl_vcf="$batch_dir/VCF/Bc.SNP.INDEL.filtered.PASS.all.vcf"
if [ ! -f "$excl_vcf" ]; then
    sbatch $filt_dependency <<EOF
#!/bin/bash
#SBATCH --job-name=PASS_Bc
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=${batch_dir}/PASS_Bc.out
#SBATCH --error=${batch_dir}/PASS_Bc.err
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call
gatk SelectVariants \
    -R "$ref" -V "$filtered_vcf" \
    --exclude-filtered true \
    -O "$excl_vcf"
EOF
fi
