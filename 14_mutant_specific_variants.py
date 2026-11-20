#!/usr/bin/env python3

vcf_file = "GenomesBGI/Bc/VCF/Bc.SNP.INDEL.filtered.PASS.all.vcf"
output_file = "GenomesBGI/Bc/VCF/mutant_specific_variants.tsv"

mutants = ["18_x_KO", "23_x_KO_#7", "28xKO2", "29XKO1"]
wt_samples = ["B_05.10", "WT", "WTB50"]

def parse_gt(gt_str):
    gt = gt_str.split(':')[0]
    if gt in ['.', './.', '.|.']:
        return None
    return '1' in gt

samples = []
sample_idx = {}
results = {m: [] for m in mutants}

with open(vcf_file) as f:
    for line in f:
        if line.startswith('#CHROM'):
            fields = line.strip().split('\t')
            samples = fields[9:]
            for i, s in enumerate(samples):
                sample_idx[s] = i
            continue

        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        sample_data = fields[9:]

        for mutant in mutants:
            mut_idx = sample_idx[mutant]
            mut_has_alt = parse_gt(sample_data[mut_idx])

            if not mut_has_alt:
                continue

            all_wt_ref = True
            for wt in wt_samples:
                wt_idx = sample_idx[wt]
                wt_has_alt = parse_gt(sample_data[wt_idx])
                if wt_has_alt:
                    all_wt_ref = False
                    break

            if all_wt_ref:
                results[mutant].append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt
                })

with open(output_file, 'w') as out:
    out.write("Mutant\tUnique_Variants_vs_AllWT\n")
    for mutant in mutants:
        out.write(f"{mutant}\t{len(results[mutant])}\n")

    out.write("\n\n=== DETAILED LIST ===\n")
    for mutant in mutants:
        out.write(f"\n--- {mutant} ({len(results[mutant])} unique variants) ---\n")
        out.write("Chr\tPos\tREF\tALT\n")
        for v in results[mutant]:
            out.write(f"{v['chrom']}\t{v['pos']}\t{v['ref']}\t{v['alt']}\n")

print("=== MUTANT-SPECIFIC VARIANTS (ALT in mutant, REF in ALL WT) ===")
for mutant in mutants:
    print(f"{mutant}: {len(results[mutant])} unique variants")

print(f"\nDetailed results saved to: {output_file}")
