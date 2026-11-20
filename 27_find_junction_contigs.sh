#!/bin/bash
#SBATCH --job-name=junction_contigs
#SBATCH --output=junction_contigs_%j.out
#SBATCH --error=junction_contigs_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL,END

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

gene_path="GenomesBGI/Bc/GOI"
spades_base="${gene_path}/local_assembly_spades"
db_dir="GenomesBGI/Bc/blastdb"
ref="GenomesBGI/Bc/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa"
query_gene="${gene_path}/BcSSP2.fa"
out_dir="${gene_path}/junction_contigs"
mkdir -p "$out_dir"

ins_chr="5"
ins_win_start=1328500
ins_win_end=1329646

ins_window_fa="${out_dir}/chr5_insertion_window.fa"
echo "Extracting candidate insertion window from reference..."
samtools faidx "$ref" "${ins_chr}:${ins_win_start}-${ins_win_end}" > "$ins_window_fa"
echo "  Window: ${ins_chr}:${ins_win_start}-${ins_win_end} ($(( ins_win_end - ins_win_start ))bp)"
echo ""

report="${out_dir}/junction_contig_report.txt"
echo "==========================================" | tee "$report"
echo "JUNCTION CONTIG SEARCH" | tee -a "$report"
echo "Looking for SPAdes contigs that contain BOTH:" | tee -a "$report"
echo "  1. BcSSP2 gene sequence" | tee -a "$report"
echo "  2. Chr5:${ins_win_start}-${ins_win_end} (candidate insertion window)" | tee -a "$report"
echo "==========================================" | tee -a "$report"
echo "" | tee -a "$report"

for sample in 18_x_KO "23_x_KO_#7" 28xKO2 29XKO1; do
    spades_dir="${spades_base}/${sample}/spades_assembly"
    contigs="${spades_dir}/contigs.fasta"

    echo "=== Sample: $sample ===" | tee -a "$report"

    if [ ! -f "$contigs" ]; then
        echo "  ⚠️  No contigs.fasta found at $spades_dir" | tee -a "$report"
        echo "     SPAdes may still be running, or assembly failed." | tee -a "$report"
        echo "" | tee -a "$report"
        continue
    fi

    n_contigs=$(grep -c "^>" "$contigs")
    echo "  Total contigs assembled: $n_contigs" | tee -a "$report"

    sample_dir="${out_dir}/${sample}"
    mkdir -p "$sample_dir"

    echo "" | tee -a "$report"
    echo "  [1] BLASTing contigs against BcSSP2..." | tee -a "$report"
    gene_blast="${sample_dir}/contigs_vs_BcSSP2.txt"

    blastn \
        -query "$contigs" \
        -subject "$query_gene" \
        -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen" \
        -perc_identity 85 \
        -evalue 1e-5 \
        -out "$gene_blast"

    if [ ! -s "$gene_blast" ]; then
        echo "  No contigs hit BcSSP2 (try lowering -perc_identity if unexpected)" | tee -a "$report"
        echo "" | tee -a "$report"
        continue
    fi

    gene_hit_contigs="${sample_dir}/gene_hit_contig_names.txt"
    awk '{print $1}' "$gene_blast" | sort -u > "$gene_hit_contigs"
    n_gene_hits=$(wc -l < "$gene_hit_contigs")
    echo "  Contigs with BcSSP2 hits: $n_gene_hits" | tee -a "$report"
    awk '{printf "    %-35s %5.1f%% identity, %dbp aligned, gene pos %d-%d\n",
        $1,$3,$4,$7,$8}' "$gene_blast" | tee -a "$report"

    echo "" | tee -a "$report"
    echo "  [2] BLASTing contigs against chr5 insertion window..." | tee -a "$report"
    win_blast="${sample_dir}/contigs_vs_ins_window.txt"

    blastn \
        -query "$contigs" \
        -subject "$ins_window_fa" \
        -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen" \
        -perc_identity 90 \
        -evalue 1e-5 \
        -out "$win_blast"

    if [ ! -s "$win_blast" ]; then
        echo "  No contigs hit the insertion window" | tee -a "$report"
    else
        win_hit_contigs="${sample_dir}/win_hit_contig_names.txt"
        awk '{print $1}' "$win_blast" | sort -u > "$win_hit_contigs"
        n_win_hits=$(wc -l < "$win_hit_contigs")
        echo "  Contigs overlapping insertion window: $n_win_hits" | tee -a "$report"
    fi

    echo "" | tee -a "$report"
    echo "  [3] Looking for contigs in BOTH hit sets (junction contigs)..." | tee -a "$report"
    junction_names="${sample_dir}/junction_contig_names.txt"

    if [ -s "$win_blast" ]; then
        comm -12 \
            <(sort "$gene_hit_contigs") \
            <(sort "$win_hit_contigs") \
            > "$junction_names"
    else
        touch "$junction_names"
    fi

    n_junction=$(wc -l < "$junction_names")

    if [ "$n_junction" -eq 0 ]; then
        echo "  ✗ No junction contigs found for $sample" | tee -a "$report"
        echo "    → The insertion breakpoint may not be spanned by a single contig." | tee -a "$report"
        echo "    → Try checking gene-hit contigs in IGV manually." | tee -a "$report"
    else
        echo "  ✓ JUNCTION CONTIGS FOUND: $n_junction" | tee -a "$report"
        echo "" | tee -a "$report"

        junction_fa="${sample_dir}/junction_contigs.fasta"
        python3 - <<PYEOF
import sys
keep = set(open("$junction_names").read().split())
out = open("$junction_fa", "w")
name = None; seq = []
for line in open("$contigs"):
    line = line.rstrip()
    if line.startswith(">"):
        if name and name in keep:
            out.write(">" + name + "\n" + "".join(seq) + "\n")
        name = line[1:].split()[0]; seq = []
    else:
        seq.append(line)
if name and name in keep:
    out.write(">" + name + "\n" + "".join(seq) + "\n")
out.close()
PYEOF

        echo "  Junction contig details:" | tee -a "$report"
        while read -r cname; do
            gene_info=$(grep "^${cname}" "$gene_blast" | \
                awk '{printf "BcSSP2 pos %d-%d (%.1f%% id, %dbp)", $7,$8,$3,$4}' | head -1)
            win_info=$(grep "^${cname}" "$win_blast" | \
                awk -v ws="$ins_win_start" \
                '{printf "chr5:%d-%d (%.1f%% id, %dbp)", $7+ws,$8+ws,$3,$4}' | head -1)
            clen=$(grep "^${cname}" "$gene_blast" | awk '{print $11}' | head -1)

            echo "    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$report"
            echo "    Contig : $cname (${clen}bp)" | tee -a "$report"
            echo "    Gene   : $gene_info" | tee -a "$report"
            echo "    Window : $win_info" | tee -a "$report"
        done < "$junction_names"

        echo "" | tee -a "$report"
        echo "  [4] Mapping junction contigs to full genome (exact insertion coordinates)..." | tee -a "$report"
        genome_blast="${sample_dir}/junction_contigs_vs_genome.txt"

        sample_db="${db_dir}/${sample}_db"
        if [ ! -f "${sample_db}.nsq" ] && [ ! -f "${sample_db}.nsi" ]; then
            echo "  ⚠️  No genome DB found at $sample_db — skipping genome mapping" | tee -a "$report"
        else
            blastn \
                -query "$junction_fa" \
                -db "$sample_db" \
                -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen" \
                -perc_identity 95 \
                -evalue 1e-10 \
                -num_threads 8 \
                -out "$genome_blast"

            echo "  Genome coordinates of junction contigs:" | tee -a "$report"
            if [ ! -s "$genome_blast" ]; then
                echo "  No high-confidence genome hits — contig may be novel/chimeric" | tee -a "$report"
            else
                awk '{printf "    %-30s → %s:%d-%d (%.1f%% id, %dbp)\n",
                    $1,$2,$7,$8,$3,$4}' "$genome_blast" \
                    | sort -u | tee -a "$report"
            fi
        fi

        echo "" | tee -a "$report"
        echo "  Junction FASTA saved: $junction_fa" | tee -a "$report"
    fi

    echo "" | tee -a "$report"
    echo "----------------------------------------" | tee -a "$report"
    echo "" | tee -a "$report"
done

echo "" | tee -a "$report"
echo "==========================================" | tee -a "$report"
echo "WHAT TO DO WITH THE RESULTS" | tee -a "$report"
echo "==========================================" | tee -a "$report"
echo "" | tee -a "$report"
echo "If junction contigs found:" | tee -a "$report"
echo "  1. Load junction_contigs.fasta into IGV alongside the KO BAM" | tee -a "$report"
echo "     Navigate to chr5:${ins_win_start}-${ins_win_end}" | tee -a "$report"
echo "  2. The contig should visually bridge the gap between the" | tee -a "$report"
echo "     insertion window and the BcSSP2 sequence" | tee -a "$report"
echo "  3. Use the genome BLAST coordinates above to design PCR primers:" | tee -a "$report"
echo "     - Forward primer: in chr5 upstream flanking sequence" | tee -a "$report"
echo "     - Reverse primer: in BcSSP2 sequence" | tee -a "$report"
echo "     → Amplifies ONLY in KO (insertion present)" | tee -a "$report"
echo "" | tee -a "$report"
echo "If NO junction contigs found:" | tee -a "$report"
echo "  → The insertion breakpoint is not spanned by any single contig" | tee -a "$report"
echo "  → Short reads cannot bridge it — long reads needed" | tee -a "$report"
echo "  → Check gene-hit contigs manually in IGV as a fallback" | tee -a "$report"
echo "" | tee -a "$report"
echo "Full report: $report" | tee -a "$report"
echo "==========================================" | tee -a "$report"
