#!/usr/bin/env python3
import re
import sys

vcf_file = "GenomesBGI/Bc/VCF/Bc.SNP.INDEL.filtered.PASS.all.snpeff.highimpact.vcf"
output_file = "GenomesBGI/Bc/VCF/highimpact_summary_table.tsv"

def parse_ann(ann_field):
    first_ann = ann_field.split(',')[0]
    parts = first_ann.split('|')
    if len(parts) >= 4:
        variant_type = parts[1]
        gene_name = parts[3]
        return gene_name, variant_type
    return "NA", "NA"

def parse_format_field(format_str, sample_data):
    format_fields = format_str.split(':')
    sample_values = sample_data.split(':')
    try:
        ad_idx = format_fields.index('AD')
        ad_value = sample_values[ad_idx]
        if ',' in ad_value:
            ref_count, alt_count = ad_value.split(',')[:2]
            return int(ref_count), int(alt_count)
    except (ValueError, IndexError):
        pass
    return 0, 0

with open(output_file, 'w') as out:
    header = ["Chr", "Position", "Gene", "Variant_Type", "REF", "ALT"]
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                for sample in samples:
                    header.extend([f"{sample}_REF", f"{sample}_ALT", f"{sample}_Total"])
                break
    out.write('\t'.join(header) + '\n')

    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            format_str = fields[8]
            sample_data = fields[9:]

            ann_match = re.search(r'ANN=([^;]+)', info)
            if ann_match:
                gene_name, variant_type = parse_ann(ann_match.group(1))
            else:
                gene_name, variant_type = "NA", "NA"

            row = [chrom, pos, gene_name, variant_type, ref, alt]
            for sample in sample_data:
                ref_count, alt_count = parse_format_field(format_str, sample)
                total = ref_count + alt_count
                row.extend([str(ref_count), str(alt_count), str(total)])
            out.write('\t'.join(row) + '\n')

print(f"✓ Summary table created: {output_file}")
print(f"  Total variants: {sum(1 for _ in open(vcf_file) if not _.startswith('#'))}")
