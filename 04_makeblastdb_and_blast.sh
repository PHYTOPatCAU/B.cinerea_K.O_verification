#!/bin/bash
#SBATCH --job-name=BcBLAST
#SBATCH --output=BcBLAST_%j.out
#SBATCH --error=BcBLAST_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --mail-type=FAIL,END

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

assembly_dir="GenomesBGI/Bc/assemblies"
db_dir="GenomesBGI/Bc/blastdb"
gene_path="GenomesBGI/Bc/GOI"
mkdir -p "$db_dir"

echo "=== Making BLAST databases for all assemblies ==="
# Make BLAST databases for all assemblies
# Priority: scaffolds.fasta (consensus) > K77/scaffolds.fasta > K55/scaffolds.fasta > contigs.fasta
for asm in "$assembly_dir"/*_assembly; do
    sample=$(basename "$asm" _assembly)
    fasta=""

    # Prefer consensus scaffolds.fasta (best assembly)
    if [ -f "$asm/scaffolds.fasta" ]; then
        fasta="$asm/scaffolds.fasta"
        echo "  - $sample: Using consensus scaffolds.fasta"
    elif [ -f "$asm/K77/scaffolds.fasta" ]; then
        fasta="$asm/K77/scaffolds.fasta"
        echo "  - $sample: Using K77 scaffolds"
    elif [ -f "$asm/K55/scaffolds.fasta" ]; then
        fasta="$asm/K55/scaffolds.fasta"
        echo "  - $sample: Using K55 scaffolds"
    elif [ -f "$asm/contigs.fasta" ]; then
        fasta="$asm/contigs.fasta"
        echo "  - $sample: Using contigs.fasta"
    else
        echo "  - $sample: No fasta found, skipping."
        continue
    fi

    makeblastdb -in "$fasta" -dbtype nucl -out "$db_dir/${sample}_db" -parse_seqids
done

echo ""
echo "Gene directory: $gene_path"
echo "BLAST DB directory: $db_dir"
echo ""

# Run BLAST for all genes against all DBs

# Get actual DB prefixes (not extensions)
db_prefixes=()
for db in "$db_dir"/*_db.nsq; do
    db_prefix="${db%.nsq}"
    db_prefixes+=("$db_prefix")
done

# Run BLAST for all genes against all DBs
for gene_file in "$gene_path"/*.fa; do
    gene_name=$(basename "$gene_file" .fa)
    hits_file="$gene_path/${gene_name}_hits.txt"

    echo ">>> Processing gene: $gene_name"
    rm -f "$hits_file"

    if [ ! -s "$gene_file" ]; then
        echo "!!! WARNING: $gene_file is empty. Skipping." >&2
        echo "Warning: $gene_file is empty, skipping BLAST." >> "$gene_path/${gene_name}_blast.err"
        continue
    fi

    # Loop over actual DB prefixes
    for db_prefix in "${db_prefixes[@]}"; do
        sample=$(basename "$db_prefix" _db)
        echo "   - Searching DB: $sample"

        # Run BLAST and capture output
        result=$(blastn \
            -query "$gene_file" \
            -db "$db_prefix" \
            -outfmt 6 \
            2>> "$gene_path/${gene_name}_blast.err")

        # Check if result is empty (no hits)
        if [ -z "$result" ]; then
            echo "   - No hits found for $sample - marking as Missing"
            echo -e "${gene_name}\t${sample}\tMissing\tMissing\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t0" >> "$hits_file"
        else
            echo "$result" | awk -v gene="$gene_name" -v sample="$sample" \
                '{print gene "\t" sample "\t" $0}' >> "$hits_file"
        fi
    done

    echo ">>> Finished gene: $gene_name (results: ${gene_name}_hitsv2.txt)"
    echo ""
done

# Summary file - create header once, then append all hits WITHOUT their headers
echo "=== Creating global summary: [all_strain_hits.txt]==="
echo -e "Gene\tSample\tQueryID\tSubjectID\tPercIdentity\tAlignLength\tMismatch\tGapOpen\tQStart\tQEnd\tSStart\tSEnd\tEvalue\tBitScore" > "$gene_path/all_strain_hits.txt"

# Append all individual hits files (they don't have headers)
for hits_file in "$gene_path"/*_hits.txt; do
    [ -f "$hits_file" ] && cat "$hits_file" >> "$gene_path/all_strainsv2_hits.txt"
done
