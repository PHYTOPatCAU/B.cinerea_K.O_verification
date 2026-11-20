#!/usr/bin/env python3
import re

highimpact_vcf = "GenomesBGI/Bc/VCF/Bc.SNP.INDEL.filtered.PASS.all.snpeff.highimpact.vcf"
full_vcf = "GenomesBGI/Bc/VCF/Bc.SNP.INDEL.filtered.PASS.all.snpeff.vcf"
output_file = "GenomesBGI/Bc/VCF/gene_list_allele_counts.tsv"

target_genes = [
    "Bcin01g03890", "Bcin01g04950", "Bcin02g00003", "Bcin03g03320",
    "Bcin03g08690", "Bcin05g01600", "Bcin05g07980", "Bcin06g00550",
    "Bcin07g06370", "Bcin08g00960", "Bcin09g03260", "Bcin11g01780",
    "Bcin12g00490", "Bcin12g06810", "Bcin13g05880", "Bcin14g01200",
    "Bcin16g02140"
]

wt_samples = ["WT", "WTB50", "B_05.10"]

def parse_ann_for_gene(ann_field):
    genes = set()
    for ann in ann_field.split(','):
        parts = ann.split('|')
        if len(parts) >= 5:
            gene_name = parts[3]
            gene_id = parts[4]
            if gene_name:
                genes.add(gene_name)
            if gene_id:
                genes.add(gene_id)
    return genes

def get_variant_type(ann_field):
    first_ann = ann_field.split(',')[0]
    parts = first_ann.split('|')
    if len(parts) >= 2:
        return parts[1]
    return "NA"

def get_gene_name_for_id(ann_field, target_gene):
    for ann in ann_field.split(','):
        parts = ann.split('|')
        if len(parts) >= 5:
            gene_name = parts[3]
            gene_id = parts[4]
            if gene_id == target_gene:
                return f"{gene_name} ({gene_id})" if gene_name else gene_id
            if gene_name == target_gene:
                return f"{gene_name} ({gene_id})" if gene_id else gene_name
    return target_gene

def parse_ad(format_str, sample_data):
    format_fields = format_str.split(':')
    sample_values = sample_data.split(':')
    try:
        ad_idx = format_fields.index('AD')
        ad_value = sample_values[ad_idx]
        if ',' in ad_value:
            ref, alt = ad_value.split(',')[:2]
            return int(ref), int(alt)
    except (ValueError, IndexError):
        pass
    return 0, 0

def search_vcf(vcf_file, target_genes, results):
    samples = []
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                continue
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref_allele = fields[3]
            alt_allele = fields[4]
            info = fields[7]
            format_str = fields[8]
            sample_data = fields[9:]

            ann_match = re.search(r'ANN=([^;]+)', info)
            if not ann_match:
                continue
            ann_field = ann_match.group(1)
            ann_genes = parse_ann_for_gene(ann_field)

            for target_gene in target_genes:
                if target_gene in ann_genes:
                    variant_type = get_variant_type(ann_field)
                    gene_display = get_gene_name_for_id(ann_field, target_gene)
                    sample_counts = {}
                    for i, sample_name in enumerate(samples):
                        ref_count, alt_count = parse_ad(format_str, sample_data[i])
                        sample_counts[sample_name] = (ref_count, alt_count)
                    key = f"{chrom}:{pos}:{target_gene}"
                    if key not in results:
                        results[key] = {
                            'chrom': chrom,
                            'pos': pos,
                            'gene': gene_display,
                            'variant_type': variant_type,
                            'ref': ref_allele,
                            'alt': alt_allele,
                            'sample_counts': sample_counts,
                            'samples': samples
                        }
    return samples

results = {}

print("Searching high-impact VCF...")
samples = search_vcf(highimpact_vcf, target_genes, results)
print(f"  Found {len(results)} variants")

print("Searching full annotated VCF...")
samples = search_vcf(full_vcf, target_genes, results)
print(f"  Total variants found: {len(results)}")

with open(output_file, 'w') as out:
    mutant_samples = [s for s in samples if s not in wt_samples]
    header = ["Chr", "Position", "Gene", "Variant_Type", "REF", "ALT"]
    for sample in mutant_samples:
        header.extend([f"{sample}_REF", f"{sample}_ALT", f"{sample}_Total"])
    header.extend(["WT_pooled_REF", "WT_pooled_ALT", "WT_pooled_Total"])
    out.write('\t'.join(header) + '\n')

    sorted_keys = sorted(results.keys(), key=lambda x: (x.split(':')[0], int(x.split(':')[1])))
    for key in sorted_keys:
        r = results[key]
        row = [r['chrom'], r['pos'], r['gene'], r['variant_type'], r['ref'], r['alt']]
        for sample in mutant_samples:
            if sample in r['sample_counts']:
                ref_c, alt_c = r['sample_counts'][sample]
            else:
                ref_c, alt_c = 0, 0
            row.extend([str(ref_c), str(alt_c), str(ref_c + alt_c)])
        wt_ref_total = 0
        wt_alt_total = 0
        for wt in wt_samples:
            if wt in r['sample_counts']:
                ref_c, alt_c = r['sample_counts'][wt]
                wt_ref_total += ref_c
                wt_alt_total += alt_c
        row.extend([str(wt_ref_total), str(wt_alt_total), str(wt_ref_total + wt_alt_total)])
        out.write('\t'.join(row) + '\n')

print(f"\n✓ Output saved to: {output_file}")
print(f"  Genes searched: {len(target_genes)}")
print(f"  Genes with variants found: {len(set(r['gene'] for r in results.values()))}")

found_genes = set(r['gene'] for r in results.values())
missing_genes = [g for g in target_genes if g not in found_genes]
if missing_genes:
    print(f"\n⚠ Genes with NO variants in VCF:")
    for g in missing_genes:
        print(f"    - {g}")
