#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})
setwd("GenomesBGI/Bc/GOI")
# ── Coordinates ──────────────────────────────────────────────────────────────
CHR          <- "Chr5"
INS_START    <- 1328501
INS_END      <- 1329621
DEL_START    <- 1329646
DEL_END      <- 1331550
PLOT_START   <- 1326500
PLOT_END     <- 1333000

# ── Sample metadata ──────────────────────────────────────────────────────────
samples_ko <- c("18_x_KO", "23_x_KO_#7", "28xKO2", "29XKO1")
samples_wt <- c("B_05.10")
all_samples <- c(samples_ko, samples_wt)

sample_colours <- c(
  "18_x_KO"    = "#C0392B",
  "23_x_KO_#7" = "#E74C3C",
  "28xKO2"     = "#E67E22",
  "29XKO1"     = "#D35400",
  "B_05.10"    = "#27AE60",
  "WT"         = "#2ECC71",
  "WTB50"      = "#1ABC9C"
)

group_colours <- c("KO" = "#C0392B", "WT" = "#27AE60")

# ── Panel A data: discordant mate counts in chr5 insertion window ─────────────
mate_data <- tibble(
  Sample  = factor(all_samples, levels = all_samples),
  Group   = c(rep("KO", 4), rep("WT", 1)),
  Mates   = c(28, 30, 12, 34,
              24)
)

# ── Panel B data: chr14 mates (mismatch signal — present only in KO) ──────────
chr14_data <- tibble(
  Sample  = factor(all_samples, levels = all_samples),
  Group   = c(rep("KO", 4), rep("WT", 1)),
  Mates   = c(84, 114, 46, 70,
               0)
)

# ── Panel C data: junction contig evidence ────────────────────────────────────
junction_data <- tibble(
  Sample       = factor(samples_ko, levels = all_samples),
  Contig_len   = c(11704, 11706, 11705, 11705),
  Gene_pct_id  = c(100.0, 100.0, 100.0, 100.0),
  Gene_aln_bp  = c(140, 140, 140, 140),
  Window_pct_id = c(100.0, 100.0, 100.0, 100.0),
  Window_aln_bp = c(1121, 1121, 1121, 1121),
  Found        = TRUE
)

# ── Shared theme ─────────────────────────────────────────────────────────────
base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text        = element_text(face = "bold"),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", size = 10),
    plot.subtitle     = element_text(size = 8, colour = "grey50"),
    axis.title        = element_text(size = 9),
    axis.text         = element_text(size = 8)
  )

# ── Panel A: mate counts in chr5 insertion window ─────────────────────────────
pA <- ggplot(mate_data, aes(x = Sample, y = Mates, fill = Group)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  annotate("text", x = 4.5, y = 36, label = "KO-specific\nsignal",
           size = 2.8, colour = "#C0392B", hjust = 0, fontface = "italic") +
  annotate("text", x = 5.5, y = 26, label = "Breakpoint\nartefact",
           size = 2.8, colour = "grey50", hjust = 0, fontface = "italic") +
  scale_fill_manual(values = group_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "A  Discordant mates landing in Chr5:1329000–1329500",
    subtitle = "KO-elevated signal marks the 5′ insertion breakpoint",
    x = NULL, y = "Mate count"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))

# ── Panel B: chr14 mates (mismatch — excluded) ────────────────────────────────
pB <- ggplot(chr14_data, aes(x = Sample, y = Mates, fill = Group)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  annotate("rect", xmin = 0.4, xmax = 4.6,
           ymin = -5, ymax = 125,
           fill = NA, colour = "#E74C3C", linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = 2.5, y = 118,
           label = "KO-specific — initially flagged, ruled out as cross-mapping",
           size = 2.6, colour = "#C0392B", hjust = 0.5) +
  scale_fill_manual(values = group_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title    = "B  Discordant mates landing in Chr14:247000–249500",
    subtitle = "Strong KO signal but absent in WT — identified as repetitive cross-mapping, not a second insertion",
    x = NULL, y = "Mate count"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))

# ── Panel C: junction contig summary ─────────────────────────────────────────
junction_long <- junction_data %>%
  select(Sample, Gene_aln_bp, Window_aln_bp) %>%
  pivot_longer(-Sample, names_to = "Region", values_to = "Aligned_bp") %>%
  mutate(Region = recode(Region,
    Gene_aln_bp   = "BcSSP2 gene sequence",
    Window_aln_bp = "Chr5 insertion window"
  ))

pC <- ggplot(junction_long,
             aes(x = Sample, y = Aligned_bp, fill = Region)) +
  geom_col(position = "dodge", width = 0.65,
           colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "BcSSP2 gene sequence"  = "#8E44AD",
    "Chr5 insertion window" = "#2980B9"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "C  SPAdes junction contig BLAST evidence (100% identity, all KO samples)",
    subtitle = "Each KO sample produced one ~11.7kb contig spanning both BcSSP2 sequence and the chr5 insertion window",
    x = NULL, y = "Aligned bp (100% identity)",
    fill = NULL
  ) +
  base_theme +
  theme(
    axis.text.x   = element_text(angle = 35, hjust = 1, size = 8),
    legend.position = "bottom",
    legend.text   = element_text(size = 8)
  )

# ── Panel D: genome diagram ───────────────────────────────────────────────────
GENE_S_WT <- 1329646
GENE_E_WT <- 1330233
DEL_S     <- 1329646
DEL_E     <- 1331550
WIN_S     <- 1328501
WIN_E     <- 1329621
FRAG_S    <- 1329481
FRAG_E    <- 1329621
PLOT_S    <- 1327800
PLOT_E    <- 1332200

WT_Y   <- 2.2
KO_Y   <- 1.1
FEAT_H <- 0.36
BAR_H  <- 0.10

draw_rect <- function(x0, x1, yc, h, col, border=NA, lwd=1) {
  rect(x0, yc - h/2, x1, yc + h/2, col=col, border=border, lwd=lwd)
}

draw_panel <- function() {
  par(mar=c(5, 1.5, 3.5, 1.5), xpd=TRUE)
  plot.new()
  plot.window(xlim=c(PLOT_S, PLOT_E), ylim=c(0.3, 3.2))

  BRACK_Y <- KO_Y + FEAT_H/2 + 0.22

  draw_rect(PLOT_S, PLOT_E, WT_Y, BAR_H, '#d0d0d0')
  draw_rect(PLOT_S, PLOT_E, KO_Y, BAR_H, '#d0d0d0')

  draw_rect(DEL_S, DEL_E, WT_Y, FEAT_H, '#5B8DB8', border='#4a7aa0', lwd=0.8)
  text((GENE_E_WT + DEL_E)/2, WT_Y, 'flanking locus',
       col='white', cex=0.65, adj=c(0.5, 0.5))

  draw_rect(GENE_S_WT, GENE_E_WT, WT_Y, FEAT_H, '#27AE60')
  text((GENE_S_WT + GENE_E_WT)/2, WT_Y, 'BcSSP2\n(587bp)',
       col='white', cex=0.75, font=2, adj=c(0.5, 0.5))

  draw_rect(DEL_S, DEL_E, KO_Y, BAR_H + 0.02, 'white')

  lines(c(DEL_S, DEL_E), c(KO_Y, KO_Y),
        col='#C0392B', lwd=2, lty=2)
  text((DEL_S + DEL_E)/2, KO_Y - 0.28,
       '~1.9kb deletion', col='#C0392B', cex=0.72, font=3)

  draw_rect(FRAG_S, FRAG_E, KO_Y, FEAT_H, '#8E44AD')
  text((FRAG_S + FRAG_E)/2, KO_Y, '140bp\nfrag.',
       col='white', cex=0.6, font=2, adj=c(0.5, 0.5))

  arrows(WIN_S, BRACK_Y, WIN_E, BRACK_Y,
         code=3, length=0.08, col='#8E44AD', lwd=1.5)
  text((WIN_S + WIN_E)/2, BRACK_Y + 0.09,
       'Discordant read insertion window (1,121bp)',
       col='#8E44AD', cex=0.68)

  lines(c(WIN_E, WIN_E), c(KO_Y + FEAT_H/2, BRACK_Y),
        col='#8E44AD', lwd=0.8, lty=3)
  lines(c(WIN_S, WIN_S), c(KO_Y + BAR_H/2, BRACK_Y),
        col='#8E44AD', lwd=0.8, lty=3)

  for (xpos in c(DEL_S, DEL_E, WIN_S)) {
    lines(c(xpos, xpos),
          c(KO_Y + FEAT_H/2 + 0.05, WT_Y - FEAT_H/2 - 0.05),
          col='#cccccc', lwd=0.7, lty=3)
  }

  text(PLOT_S - 40, WT_Y, 'WT', col='#27AE60', cex=1.1, font=2, adj=c(1, 0.5))
  text(PLOT_S - 40, KO_Y, 'KO', col='#C0392B', cex=1.1, font=2, adj=c(1, 0.5))

  tick_pos <- seq(1328000, 1332000, by=500)
  segments(PLOT_S, 0.58, PLOT_E, 0.58, col='#bbbbbb', lwd=0.6)
  for (t in tick_pos) {
    segments(t, 0.52, t, 0.58, col='#aaaaaa', lwd=0.7)
    text(t, 0.44, paste0(t/1000, 'kb'), cex=0.65, col='#666666')
  }
  text((PLOT_S + PLOT_E)/2, 0.32, 'Chr5 position', cex=0.82, col='#555555')

  lx <- PLOT_E - 50
  ly_start <- 3.05
  ly_gap   <- 0.18

  rect(lx - 420, ly_start - 0.07, lx - 360, ly_start + 0.07,
       col='#27AE60', border=NA)
  text(lx - 350, ly_start, 'BcSSP2 gene (587bp)',
       cex=0.68, adj=c(0, 0.5), col='#333333')

  ly2 <- ly_start - ly_gap
  rect(lx - 420, ly2 - 0.07, lx - 360, ly2 + 0.07,
       col='#5B8DB8', border=NA)
  text(lx - 350, ly2, 'Flanking locus context',
       cex=0.68, adj=c(0, 0.5), col='#333333')

  ly3 <- ly2 - ly_gap
  rect(lx - 420, ly3 - 0.07, lx - 360, ly3 + 0.07,
       col='#8E44AD', border=NA)
  text(lx - 350, ly3, 'Reinserted fragment (140bp)',
       cex=0.68, adj=c(0, 0.5), col='#333333')

  ly4 <- ly3 - ly_gap
  lines(c(lx - 420, lx - 360), c(ly4, ly4),
        col='#C0392B', lwd=2, lty=2)
  text(lx - 350, ly4, '~1.9kb deletion (KO)',
       cex=0.68, adj=c(0, 0.5), col='#333333')

  rect(lx - 430, ly4 - 0.10, lx - 10, ly_start + 0.11,
       col=NA, border='#cccccc', lwd=0.6)

  mtext('D   Schematic — BcSSP2 locus rearrangement on Chr5 (to scale)',
        side=3, line=2.2, adj=0, font=2, cex=1.0)
  mtext(paste0(
    'WT: BcSSP2 intact at original locus (587bp gene within ~1.9kb locus context). ',
    'KO: original locus deleted (~1.9kb);\n',
    'a 140bp BcSSP2 fragment reinserted within a 1,121bp window only ~25bp upstream ',
    'of the deletion breakpoint, in all 4 independent KO lines.'),
    side=3, line=0.1, adj=0, cex=0.68, col='#666666')
}

# ── Save PDF ──────────────────────────────────────────────────────────────────
pdf('panelD_schematic.pdf', width=10, height=4.8)
draw_panel()
dev.off()
cat('Saved: panelD_schematic.pdf\n')

# ── Save PNG ──────────────────────────────────────────────────────────────────
png('panelD_schematic.png', width=10, height=4.8,
    units='in', res=300, bg='white')
draw_panel()
dev.off()
cat('Saved: panelD_schematic.png\n')

pdf("BcSSP2_insertion_panelC.pdf", width = 6, height = 5)
print(pC)
dev.off()

pdf("BcSSP2_insertion_panelD.pdf", width = 10, height = 5)
draw_panel()
dev.off()
