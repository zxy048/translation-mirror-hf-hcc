library(ggplot2)
library(grid)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# Simple 3-color palette
red  <- "#C62828"
blue <- "#1565C0"
purp <- "#7B1FA2"
grey <- "grey30"

# Set up 16:9 canvas with grid units
# Layout: [HCC panel | MIDDLE | HF panel]
# HCC: x=0.1-0.42, MID: x=0.42-0.58, HF: x=0.58-0.9
# Y ranges: 0.1 to 0.9

p <- ggplot() +
  # Title
  annotate("text", x=0.5, y=0.95, label="Mirror Regulation of Translation in HCC vs Heart Failure",
           fontface="bold", size=7, color="grey20") +
  annotate("text", x=0.5, y=0.90, label="Shared architecture, opposed direction, distinct upstream programs",
           size=4, color="grey50") +

  # Disease labels
  annotate("text", x=0.26, y=0.86, label="HCC", fontface="bold", size=10, color=red) +
  annotate("text", x=0.74, y=0.86, label="HF", fontface="bold", size=10, color=blue) +
  annotate("text", x=0.26, y=0.82, label="Hepatocellular Carcinoma", size=4.5, color=red) +
  annotate("text", x=0.74, y=0.82, label="Heart Failure", size=4.5, color=blue) +

  # Vertical divider
  annotate("segment", x=0.50, xend=0.50, y=0.12, yend=0.79, linetype="dashed", color="grey70", linewidth=1) +
  annotate("text", x=0.50, y=0.81, label="MIRROR", size=3.5, color="grey60", fontface="italic") +

  # ===== HCC BOX 1: Upstream (y: 0.64-0.77, x: 0.10-0.42) =====
  annotate("rect", xmin=0.10, xmax=0.42, ymin=0.62, ymax=0.78, fill="#FFCDD2", color=red, linewidth=1.2) +
  annotate("text", x=0.26, y=0.74, label="MYC Program  +  ATF4 / ISR", fontface="bold", size=5, color="#4A0000") +
  annotate("text", x=0.26, y=0.69, label="MYC pathway rho=+0.613   ATF4 rho=+0.439", size=3.5, color=red) +
  annotate("text", x=0.26, y=0.65, label="DDIT3  XBP1  HIF1A feedback", size=3, color="#BF360C", fontface="italic") +

  # ===== HCC BOX 2: Translation (y: 0.42-0.58, x: 0.10-0.42) =====
  annotate("rect", xmin=0.10, xmax=0.42, ymin=0.40, ymax=0.58, fill="#EF9A9A", color=red, linewidth=1.3) +
  annotate("text", x=0.26, y=0.54, label="Ribosome Biogenesis", fontface="bold", size=5.5, color="#4A0000") +
  annotate("text", x=0.26, y=0.49, label="EEF1A1  FAU  RPL39  RPL3  RPL32  RPL41  RPS28", size=3.8, color="#4A0000") +
  annotate("text", x=0.26, y=0.44, label="Translation Initiation  |  Elongation  |  rRNA Processing", size=3.2, color="#4A0000", fontface="italic") +

  # ===== HCC BOX 3: Output (y: 0.26-0.36, x: 0.13-0.39) =====
  annotate("rect", xmin=0.13, xmax=0.39, ymin=0.24, ymax=0.36, fill=red, color="#8B0000", linewidth=1.5) +
  annotate("text", x=0.26, y=0.30, label="TRANSLATION  UP", fontface="bold", size=6, color="white") +

  annotate("text", x=0.26, y=0.19, label="Pro-survival stress adaptation + Proliferation  ->  HCC growth",
           size=3.2, color="#B71C1C", fontface="italic") +

  # ===== HF BOX 1: Upstream (y: 0.64-0.77, x: 0.58-0.90) =====
  annotate("rect", xmin=0.58, xmax=0.90, ymin=0.62, ymax=0.78, fill="#BBDEFB", color=blue, linewidth=1.2) +
  annotate("text", x=0.74, y=0.74, label="Energy Depletion  +  mTORC1 Down", fontface="bold", size=5, color="#002040") +
  annotate("text", x=0.74, y=0.69, label="Metabolic stress   MTOR rho=-0.391", size=3.5, color=blue) +
  annotate("text", x=0.74, y=0.65, label="ATP-costly synthesis suppressed as adaptation", size=3, color="#0D47A1", fontface="italic") +

  # ===== HF BOX 2: Translation (y: 0.42-0.58, x: 0.58-0.90) =====
  annotate("rect", xmin=0.58, xmax=0.90, ymin=0.40, ymax=0.58, fill="#90CAF9", color=blue, linewidth=1.3) +
  annotate("text", x=0.74, y=0.54, label="Ribosome & Translation", fontface="bold", size=5.5, color="#002040") +
  annotate("text", x=0.74, y=0.49, label="Genes Transcriptionally Downregulated", size=3.8, color="#002040") +
  annotate("text", x=0.74, y=0.44, label="7 Hub Genes Discordant  |  ssGSEA pathways down vs normal", size=3.2, color="#002040", fontface="italic") +

  # ===== HF BOX 3: Output (y: 0.26-0.36, x: 0.61-0.87) =====
  annotate("rect", xmin=0.61, xmax=0.87, ymin=0.24, ymax=0.36, fill=blue, color="#002040", linewidth=1.5) +
  annotate("text", x=0.74, y=0.30, label="TRANSLATION  DOWN", fontface="bold", size=6, color="white") +

  annotate("text", x=0.74, y=0.19, label="Energy conservation  ->  Cardiomyocyte survival",
           size=3.2, color="#0D47A1", fontface="italic") +

  # ===== MIDDLE PANELS =====
  annotate("rect", xmin=0.389, xmax=0.611, ymin=0.53, ymax=0.64, fill="#E1BEE7", color=purp, linewidth=1.3) +
  annotate("text", x=0.50, y=0.61, label="ssGSEA", fontface="bold", size=4.5, color="#4A148C") +
  annotate("text", x=0.50, y=0.58, label="Pathway ES", size=3.5, color="#4A148C") +
  annotate("text", x=0.50, y=0.555, label="rho=-0.598", fontface="bold", size=3.5, color="#4A148C") +
  annotate("text", x=0.50, y=0.535, label="p=0.0003", size=3, color=purp) +

  annotate("rect", xmin=0.389, xmax=0.611, ymin=0.38, ymax=0.49, fill="#CE93D8", color=purp, linewidth=1.3) +
  annotate("text", x=0.50, y=0.47, label="WGCNA", fontface="bold", size=4.5, color="#4A148C") +
  annotate("text", x=0.50, y=0.44, label="Module", size=3.5, color="#4A148C") +
  annotate("text", x=0.50, y=0.415, label="p=6.3e-12", fontface="bold", size=3.5, color="#4A148C") +
  annotate("text", x=0.50, y=0.395, label="Replicated", size=3, color=purp) +

  # ===== ARROWS (vertical within HCC/HF) =====
  # HCC down arrows
  annotate("segment", x=0.26, xend=0.26, y=0.615, yend=0.585, arrow=arrow(length=unit(0.06,"npc")), color=red, linewidth=1.5) +
  annotate("segment", x=0.26, xend=0.26, y=0.395, yend=0.365, arrow=arrow(length=unit(0.08,"npc")), color="white", linewidth=2.5) +

  # HF down arrows
  annotate("segment", x=0.74, xend=0.74, y=0.615, yend=0.585, arrow=arrow(length=unit(0.06,"npc")), color=blue, linewidth=1.5) +
  annotate("segment", x=0.74, xend=0.74, y=0.395, yend=0.365, arrow=arrow(length=unit(0.08,"npc")), color="white", linewidth=2.5) +

  # Lateral: HCC/HF <-> middle
  annotate("segment", x=0.42, xend=0.389, y=0.585, yend=0.585, arrow=arrow(length=unit(0.04,"npc")), color=purp, linewidth=0.5, linetype="dotted") +
  annotate("segment", x=0.58, xend=0.611, y=0.585, yend=0.585, arrow=arrow(length=unit(0.04,"npc")), color=purp, linewidth=0.5, linetype="dotted") +
  annotate("segment", x=0.42, xend=0.389, y=0.435, yend=0.435, arrow=arrow(length=unit(0.04,"npc")), color=purp, linewidth=0.5, linetype="dotted") +
  annotate("segment", x=0.58, xend=0.611, y=0.435, yend=0.435, arrow=arrow(length=unit(0.04,"npc")), color=purp, linewidth=0.5, linetype="dotted") +

  # (ssGSEA→WGCNA arrow removed — independent analyses, no causal link)

  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  theme_void() +
  theme(plot.background=element_rect(fill="white", color=NA))

ggsave(file.path(PROJ_DIR, "figures", "Figure6_Mechanistic_Model.png"), p, width=20, height=12, dpi=300)
cat("Done\n")
