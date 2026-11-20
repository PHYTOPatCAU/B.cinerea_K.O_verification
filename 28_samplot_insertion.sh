#!/bin/bash
#SBATCH --job-name=samplot_insertion
#SBATCH --output=samplot_insertion_%j.out
#SBATCH --error=samplot_insertion_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL,END

module load gcc12-env
module load miniconda3
conda activate samplot

echo "Using samplot: $(which samplot)"
echo "Python version: $(python --version)"
echo ""

ref="GenomesBGI/Bc/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa"
bam_dir="GenomesBGI/Bc/BAM/MD"
annotation_bed="GenomesBGI/Bc/GOI/genes.bed.gz"
out_dir="GenomesBGI/Bc/GOI/samplot_insertion"
mkdir -p "$out_dir"

bam_18="${bam_dir}/18_x_KO.DuplMark.bam"
bam_23="${bam_dir}/23_x_KO_#7.DuplMark.bam"
bam_28="${bam_dir}/28xKO2.DuplMark.bam"
bam_29="${bam_dir}/29XKO1.DuplMark.bam"
bam_b05="${bam_dir}/B_05.10.DuplMark.bam"
bam_wt="${bam_dir}/WT.DuplMark.bam"
bam_wtb50="${bam_dir}/WTB50.DuplMark.bam"

chr="5"
ins_start=1328501
ins_end=1329621
del_start=1329646
del_end=1331550

wide_start=$(( ins_start - 500 ))
wide_end=$(( del_end + 500 ))
zoom_start=$(( ins_start - 300 ))
zoom_end=$(( ins_end + 300 ))

echo "Key coordinates:"
echo "  Insertion window : ${chr}:${ins_start}-${ins_end}"
echo "  Original deletion: ${chr}:${del_start}-${del_end}"
echo "  Wide plot window : ${chr}:${wide_start}-${wide_end}"
echo "  Zoom plot window : ${chr}:${zoom_start}-${zoom_end}"
echo ""

run_samplot() {
    local outpng="$1"
    local pstart="$2"
    local pend="$3"
    shift 3

    local cmd=(samplot plot
        -c "$chr"
        -s "$pstart"
        -e "$pend"
        -o "$outpng"
        -r "$ref"
        -t "DEL"
        --max_depth 0
    )

    if [ -f "$annotation_bed" ]; then
        cmd+=(-A "$annotation_bed" --annotation_filenames "Genes")
    fi

    cmd+=("$@")

    "${cmd[@]}"
}

echo "=== Generating samplot figures ==="
echo ""

echo "[1/4] Wide view: full rearrangement context (${chr}:${wide_start}-${wide_end})"
run_samplot \
    "${out_dir}/01_wide_all_samples.png" \
    "$wide_start" "$wide_end" \
    -b "$bam_18" "$bam_23" "$bam_28" "$bam_29" "$bam_b05" "$bam_wt" "$bam_wtb50" \
    -n "18_x_KO" "23_x_KO_#7" "28xKO2" "29XKO1" "B_05.10" "WT" "WTB50" \
    && echo "  ✓ 01_wide_all_samples.png" || echo "  ✗ FAILED"
echo ""

echo "[2/4] Zoom: insertion breakpoint (${chr}:${zoom_start}-${zoom_end})"
run_samplot \
    "${out_dir}/02_zoom_insertion_breakpoint.png" \
    "$zoom_start" "$zoom_end" \
    -b "$bam_18" "$bam_23" "$bam_28" "$bam_29" "$bam_b05" "$bam_wt" "$bam_wtb50" \
    -n "18_x_KO" "23_x_KO_#7" "28xKO2" "29XKO1" "B_05.10" "WT" "WTB50" \
    && echo "  ✓ 02_zoom_insertion_breakpoint.png" || echo "  ✗ FAILED"
echo ""

echo "[3/4] KO vs WT — wide view"
run_samplot \
    "${out_dir}/03_KOvsWT_wide.png" \
    "$wide_start" "$wide_end" \
    -b "$bam_23" "$bam_b05" \
    -n "23_x_KO_#7" "B_05.10 (WT)" \
    && echo "  ✓ 03_KOvsWT_wide.png" || echo "  ✗ FAILED"
echo ""

echo "[4/4] KO vs WT — breakpoint zoom"
run_samplot \
    "${out_dir}/04_KOvsWT_zoom_breakpoint.png" \
    "$zoom_start" "$zoom_end" \
    -b "$bam_23" "$bam_b05" \
    -n "23_x_KO_#7" "B_05.10 (WT)" \
    && echo "  ✓ 04_KOvsWT_zoom_breakpoint.png" || echo "  ✗ FAILED"
echo ""

echo "=== Combining into PDF ==="
python3 - << 'PYEOF'
import os, glob
out_dir = os.environ.get("out_dir", "GenomesBGI/Bc/GOI/samplot_insertion")
pngs = [
    f"{out_dir}/01_wide_all_samples.png",
    f"{out_dir}/02_zoom_insertion_breakpoint.png",
    f"{out_dir}/03_KOvsWT_wide.png",
    f"{out_dir}/04_KOvsWT_zoom_breakpoint.png",
]
pngs = [p for p in pngs if os.path.exists(p)]
if not pngs:
    print("  No PNGs found to combine")
else:
    try:
        from PIL import Image
        images = [Image.open(p).convert("RGB") for p in pngs]
        pdf_path = f"{out_dir}/BcSSP2_insertion_summary.pdf"
        images[0].save(pdf_path, save_all=True, append_images=images[1:])
        print(f"  ✓ BcSSP2_insertion_summary.pdf ({len(images)} pages)")
    except ImportError:
        print("  Pillow not available — install with: pip install Pillow")
        print(f"  PNGs are in {out_dir}")
PYEOF

echo ""
echo "=== Done ==="
echo "Figures saved to: $out_dir/"
echo ""
echo "Figure guide:"
echo "  01_wide_all_samples.png         — all 7 samples, full rearrangement"
echo "  02_zoom_insertion_breakpoint    — all 7 samples, zoomed to 5' breakpoint"
echo "  03_KOvsWT_wide.png              — 23_x_KO vs B_05.10 (cleanest for paper)"
echo "  04_KOvsWT_zoom_breakpoint.png   — 23_x_KO vs B_05.10, breakpoint zoom"
echo "  BcSSP2_insertion_summary.pdf    — all four combined (if Pillow available)"
echo ""
echo "What to look for:"
echo "  - KO: soft-clipped reads at ${chr}:${ins_end} (insertion 5' breakpoint)"
echo "  - KO: coverage drop across ${chr}:${del_start}-${del_end} (deletion)"
echo "  - WT: clean continuous coverage throughout, no anomalies"
