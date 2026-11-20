#!/bin/bash
#SBATCH --job-name=BcBLAST
#SBATCH --error=BcBLAST.err
#SBATCH --mail-type=FAIL,END

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

cd GenomesBGI/Bc

db_dir="GenomesBGI/Bc/blastdb"
db1="${db_dir}/B_05.10_db"
db2="${db_dir}/23_x_KO_#7_db"
db3="${db_dir}/18_x_KO_db"

gene_path="GOI/"

for gene_file in "$gene_path"/*.fa; do
    gene_name=$(basename "$gene_file" .fa)

    blastn -query "$gene_file" -db "$db1" -outfmt 6 -max_target_seqs 1 | awk -v gene="$gene_name" '{print gene "\t" $0}' >> "${gene_name}_hits_db1.txt"
    blastn -query "$gene_file" -db "$db2" -outfmt 6 -max_target_seqs 1 | awk -v gene="$gene_name" '{print gene "\t" $0}' >> "${gene_name}_hits_db2.txt"
    blastn -query "$gene_file" -db "$db3" -outfmt 6 -max_target_seqs 1 | awk -v gene="$gene_name" '{print gene "\t" $0}' >> "${gene_name}_hits_db3.txt"

done
