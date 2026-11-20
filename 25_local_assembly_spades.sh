#!/bin/bash
#SBATCH --job-name=local_spades
#SBATCH --output=local_spades_%j.out
#SBATCH --error=local_spades_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL,END

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

bam_dir="GenomesBGI/Bc/BAM/MD"
gene_path="GenomesBGI/Bc/GOI"
db_dir="GenomesBGI/Bc/blastdb"
out_dir="${gene_path}/local_assembly_spades"
query_gene="${gene_path}/BcSSP2.fa"
mkdir -p "$out_dir"

chr="5"
del_start=1329646
del_end=1331550
flank=5000
wide_start=$((del_start - flank))
wide_end=$((del_end + flank))
wide_region="${chr}:${wide_start}-${wide_end}"

echo "=========================================="
echo "LOCAL SPADES ASSEMBLY + BLAST"
echo "Region: $wide_region (deletion ± ${flank}bp)"
echo "Query gene: $query_gene"
echo "=========================================="
echo ""

for sample in 18_x_KO "23_x_KO_#7" 28xKO2 29XKO1; do
    bam="${bam_dir}/${sample}.DuplMark.bam"
    [ ! -f "$bam" ] && echo "WARNING: $bam not found, skipping" && continue

    echo "=== Sample: $sample ==="
    sample_dir="${out_dir}/${sample}"
    mkdir -p "$sample_dir"

    echo "[1] Extracting reads from $wide_region..."
    local_bam="${sample_dir}/local_region.bam"
    samtools view -b "$bam" "$wide_region" | samtools sort -o "$local_bam"
    samtools index "$local_bam"
    n_reads=$(samtools view -c "$local_bam")
    echo "    Extracted $n_reads reads"

    if [ "$n_reads" -lt 10 ]; then
        echo "    Too few reads for assembly (<10), skipping SPAdes"
        echo ""
        continue
    fi

    echo "[2] Converting to FASTQ..."
    fq1="${sample_dir}/reads_R1.fastq"
    fq2="${sample_dir}/reads_R2.fastq"
    fq_unpaired="${sample_dir}/reads_unpaired.fastq"
    samtools sort -n "$local_bam" -o "${sample_dir}/local_namesorted.bam"
    samtools fastq \
        -1 "$fq1" \
        -2 "$fq2" \
        -s "$fq_unpaired" \
        -0 /dev/null \
        "${sample_dir}/local_namesorted.bam"
    echo "    Done"

    echo "[3] Running SPAdes assembly (careful mode)..."
    spades_dir="${sample_dir}/spades_assembly"
    rm -rf "$spades_dir"

    spades.py \
        -1 "$fq1" \
        -2 "$fq2" \
        -s "$fq_unpaired" \
        -o "$spades_dir" \
        --careful \
        -t 8 \
        -m 60 \
        2>&1 | tail -5

    if [ ! -f "${spades_dir}/contigs.fasta" ]; then
        echo "    ⚠️  SPAdes failed or produced no contigs"
        echo ""
        continue
    fi

    n_contigs=$(grep -c "^>" "${spades_dir}/contigs.fasta")
    echo "    Assembled $n_contigs contigs"

    echo "[4] BLASTing local contigs against BcSSP2..."
    blast_out="${sample_dir}/contigs_vs_BcSSP2.txt"

    blastn \
        -query "${spades_dir}/contigs.fasta" \
        -subject "$query_gene" \
        -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen" \
        -perc_identity 80 \
        -evalue 1e-5 \
        -out "$blast_out"

    echo ""
    echo "    === BLAST HITS: Contigs containing BcSSP2 sequence ==="
    if [ ! -s "$blast_out" ]; then
        echo "    No hits found (BcSSP2 not in assembled contigs)"
    else
        echo -e "    Contig\t\t\t%ID\tAlnLen\tContigPos\tGenePos\tEval\tScore\tContigLen\tGeneLen"
        awk '{printf "    %-30s %5.1f%%\t%d\t%d-%d\t%d-%d\t%s\t%s\t%s\t%s\n",
            $1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "$blast_out"
    fi

    echo ""
    echo "[5] Mapping assembled contigs back to genome (finding insertion coordinates)..."
    genome_blast_out="${sample_dir}/contigs_vs_genome.txt"

    sample_db="${db_dir}/${sample}_db"
    if [ ! -f "${sample_db}.nsq" ]; then
        echo "    No per-sample DB found at $sample_db, skipping genome mapping"
    else
        if [ -s "$blast_out" ]; then
            hit_contigs="${sample_dir}/hit_contig_names.txt"
            awk '{print $1}' "$blast_out" | sort -u > "$hit_contigs"

            hit_fa="${sample_dir}/hit_contigs.fasta"
            python3 -c "
import sys
keep = set(open('$hit_contigs').read().split())
out = open('$hit_fa', 'w')
name = None; seq = []
for line in open('${spades_dir}/contigs.fasta'):
    line = line.rstrip()
    if line.startswith('>'):
        if name and name in keep:
            out.write('>'+name+'\n'+''.join(seq)+'\n')
        name = line[1:].split()[0]; seq = []
    else:
        seq.append(line)
if name and name in keep:
    out.write('>'+name+'\n'+''.join(seq)+'\n')
out.close()
"
            blastn \
                -query "$hit_fa" \
                -db "$sample_db" \
                -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
                -perc_identity 95 \
                -evalue 1e-10 \
                -out "$genome_blast_out"

            echo "    Genomic coordinates of BcSSP2-containing contigs:"
            if [ -s "$genome_blast_out" ]; then
                awk '{printf "    Contig %-25s → scaffold %s:%d-%d (%.1f%% id)\n",
                    $1,$2,$7,$8,$3}' "$genome_blast_out" | sort -u
            else
                echo "    No high-confidence genomic hits found"
            fi
        fi
    fi

    echo ""
    echo "    All results saved: $sample_dir/"
    echo "----------------------------------------"
    echo ""
done

echo ""
echo "=========================================="
echo "INTERPRETATION"
echo "=========================================="
echo ""
echo "What to look for in BLAST output:"
echo "  - A contig with BcSSP2 hit that maps to a DIFFERENT position"
echo "    than Chr5:1329646-1331550 → candidate reinsertion site"
echo "  - A contig containing BOTH BcSSP2 sequence AND flanking genomic"
echo "    sequence from a new location → breakpoint contig (most useful)"
echo "  - Multiple KO samples with the same new coordinate → confirmed"
echo ""
echo "If contigs found with BcSSP2 hits, check in IGV:"
echo "  → Load: ${out_dir}/<sample>/local_region.bam"
echo "  → Navigate to the new coordinates from genome BLAST"
echo "=========================================="
