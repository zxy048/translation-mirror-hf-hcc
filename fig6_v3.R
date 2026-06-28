library(ggplot2)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# Colors
h_red   <- "#C62828"
h_light <- "#FFCDD2"
h_mid   <- "#EF9A9A"
h_dark  <- "#B71C1C"
b_blue  <- "#1565C0"
b_light <- "#BBDEFB"
b_mid   <- "#90CAF9"
b_dark  <- "#0D47A1"
m_purp  <- "#7B1FA2"
m_light <- "#E1BEE7"

p <- ggplot()

# Title
p <- p + annotate("text", x=0, y=16.5, label="Mirror Regulation of Translation in HCC and Heart Failure", fontface="bold", size=8, color="grey20")
p <- p + annotate("text", x=0, y=15.6, label="Shared co-expression architecture, opposed perturbation direction, distinct upstream regulators", size=4.5, color="grey50")

# Disease headers
p <- p + annotate("text", x=-3.5, y=14.8, label="HCC", fontface="bold", size=10, color=h_red)
p <- p + annotate("text", x=3.5, y=14.8, label="HF", fontface="bold", size=10, color=b_blue)
p <- p + annotate("text", x=0, y=14.8, label="vs", size=4.5, color="grey60")
p <- p + annotate("text", x=-3.5, y=13.8, label="Hepatocellular Carcinoma", size=4, color=h_red)
p <- p + annotate("text", x=3.5, y=13.8, label="Dilated / Ischemic Cardiomyopathy", size=4, color=b_blue)

# Divider
p <- p + annotate("segment", x=0, xend=0, y=2, yend=13.2, linetype="dashed", color="grey70", linewidth=0.8)
p <- p + annotate("text", x=0, y=13.5, label="M I R R O R", size=3.8, color="grey60", fontface="italic")

# =============== LEFT: HCC ===============

# Upstream TFs
p <- p + annotate("rect", xmin=-8, xmax=-1, ymin=11.5, ymax=12.8, fill=h_light, color=h_red, linewidth=1)
p <- p + annotate("text", x=-6, y=12.3, label="ATF4 / ISR", fontface="bold", size=5, color=h_dark)
p <- p + annotate("text", x=-6, y=11.85, label="TF rho = +0.439  FDR < 0.0001", size=3.3, color=h_red)
p <- p + annotate("text", x=-3, y=12.3, label="MYC Program", fontface="bold", size=5, color=h_dark)
p <- p + annotate("text", x=-3, y=11.85, label="Pathway rho = +0.613  p < 0.0001", size=3.3, color=h_red)

# Stress mediators
p <- p + annotate("rect", xmin=-6.5, xmax=-0.5, ymin=10.0, ymax=11.2, fill="#FFE0B2", color="#E65100", linewidth=0.8)
p <- p + annotate("text", x=-3.5, y=10.7, label="DDIT3 (CHOP) rho=+0.289  |  XBP1 rho=+0.162  |  HIF1A rho=-0.417", size=3.2, color="#BF360C")
p <- p + annotate("text", x=-3.5, y=10.25, label="Integrated Stress Response / Unfolded Protein Response / Hypoxia Feedback", size=2.8, color="#E65100", fontface="italic")

# Translation
p <- p + annotate("rect", xmin=-8, xmax=-1, ymin=7.2, ymax=9.5, fill=h_mid, color=h_red, linewidth=1.2)
p <- p + annotate("text", x=-4.5, y=9.0, label="Ribosome Biogenesis", fontface="bold", size=5, color="#4A0000")
p <- p + annotate("text", x=-4.5, y=8.4, label="Translation Initiation & Elongation", size=3.8, color="#4A0000")
p <- p + annotate("text", x=-4.5, y=7.8, label="EEF1A1     FAU     RPL39     RPL3", size=3.2, color="#4A0000", fontface="italic")
p <- p + annotate("text", x=-4.5, y=7.45, label="RPL32     RPL41     RPS28", size=3.2, color="#4A0000", fontface="italic")

# Output
p <- p + annotate("rect", xmin=-7, xmax=-2, ymin=5.0, ymax=6.5, fill=h_red, color=h_dark, linewidth=1.5)
p <- p + annotate("text", x=-4.5, y=5.75, label="TRANSLATION  UP", fontface="bold", size=6, color="white")

p <- p + annotate("text", x=-4.5, y=4.3, label="Pro-survival ISR adaptation", size=3.3, color=h_dark, fontface="italic")
p <- p + annotate("text", x=-4.5, y=3.8, label="+ MYC-driven proliferation  -> HCC growth", size=3.3, color=h_dark, fontface="italic")

# =============== RIGHT: HF ===============

# Upstream stress
p <- p + annotate("rect", xmin=1, xmax=8, ymin=11.5, ymax=12.8, fill=b_light, color=b_blue, linewidth=1)
p <- p + annotate("text", x=3, y=12.3, label="Energy Depletion", fontface="bold", size=5, color=b_dark)
p <- p + annotate("text", x=3, y=11.85, label="Metabolic Stress", size=3.3, color=b_blue)
p <- p + annotate("text", x=6, y=12.3, label="mTORC1  Down", fontface="bold", size=5, color=b_dark)
p <- p + annotate("text", x=6, y=11.85, label="MTOR rho = -0.391  FDR < 0.0001", size=3.3, color=b_blue)

# Adaptation
p <- p + annotate("rect", xmin=1.5, xmax=7.5, ymin=10.0, ymax=11.2, fill="#E3F2FD", color=b_dark, linewidth=0.8)
p <- p + annotate("text", x=4.5, y=10.7, label="ATP-costly protein synthesis suppressed as adaptive response", size=3.2, color=b_dark)
p <- p + annotate("text", x=4.5, y=10.25, label="Cardiomyocyte energy conservation under chronic stress", size=2.8, color=b_dark, fontface="italic")

# Translation
p <- p + annotate("rect", xmin=1, xmax=8, ymin=7.2, ymax=9.5, fill=b_mid, color=b_blue, linewidth=1.2)
p <- p + annotate("text", x=4.5, y=9.0, label="Ribosome & Translation", fontface="bold", size=5, color="#002040")
p <- p + annotate("text", x=4.5, y=8.4, label="Genes Transcriptionally Downregulated", size=3.8, color="#002040")
p <- p + annotate("text", x=4.5, y=7.8, label="7 Hub Genes Direction Discordant Across Diseases", size=3.2, color="#002040", fontface="italic")
p <- p + annotate("text", x=4.5, y=7.45, label="ssGSEA: Translation pathways consistently down vs normal", size=2.8, color="#002040")

# Output
p <- p + annotate("rect", xmin=2, xmax=7, ymin=5.0, ymax=6.5, fill=b_blue, color=b_dark, linewidth=1.5)
p <- p + annotate("text", x=4.5, y=5.75, label="TRANSLATION  DOWN", fontface="bold", size=6, color="white")

p <- p + annotate("text", x=4.5, y=4.3, label="Energy conservation strategy", size=3.3, color=b_dark, fontface="italic")
p <- p + annotate("text", x=4.5, y=3.8, label="Cardiomyocyte survival adaptation -> HF progression", size=3.3, color=b_dark, fontface="italic")

# =============== MIDDLE INTEGRATION ===============

p <- p + annotate("rect", xmin=-2.5, xmax=2.5, ymin=10.0, ymax=11.2, fill=m_light, color=m_purp, linewidth=1.2)
p <- p + annotate("text", x=0, y=10.8, label="Cross-Disease Integration", fontface="bold", size=4, color=m_purp)
p <- p + annotate("text", x=0, y=10.4, label="Pathway Effect Size Negative Correlation", size=3, color=m_purp)
p <- p + annotate("text", x=0, y=10.1, label="ssGSEA  rho = -0.598, p = 0.0003", fontface="bold", size=2.8, color="#4A148C")

p <- p + annotate("rect", xmin=-2.5, xmax=2.5, ymin=8.2, ymax=9.7, fill="#CE93D8", color=m_purp, linewidth=1.2)
p <- p + annotate("text", x=0, y=9.35, label="Co-expression Module", fontface="bold", size=4, color="#4A148C")
p <- p + annotate("text", x=0, y=8.9, label="Conserved Across Diseases", size=3, color="#4A148C")
p <- p + annotate("text", x=0, y=8.5, label="WGCNA Replication", fontface="bold", size=2.8, color="#4A148C")
p <- p + annotate("text", x=0, y=8.25, label="p = 6.3 x 10^-12", fontface="bold", size=2.5, color="#4A148C")

# =============== ARROWS ===============

# HCC: upstream -> stress -> translation -> output
p <- p + annotate("segment", x=-6, xend=-6, y=11.4, yend=11.25, arrow=arrow(length=unit(0.08,"cm")), color=h_red, linewidth=0.8)
p <- p + annotate("segment", x=-3, xend=-3, y=11.4, yend=11.25, arrow=arrow(length=unit(0.08,"cm")), color=h_red, linewidth=0.8)
p <- p + annotate("segment", x=-3.5, xend=-3.5, y=9.9, yend=9.55, arrow=arrow(length=unit(0.1,"cm")), color=h_red, linewidth=1.2)
p <- p + annotate("segment", x=-4.5, xend=-4.5, y=7.1, yend=6.55, arrow=arrow(length=unit(0.12,"cm")), color="white", linewidth=2)

# HF: upstream -> adaptation -> translation -> output
p <- p + annotate("segment", x=3, xend=3, y=11.4, yend=11.25, arrow=arrow(length=unit(0.08,"cm")), color=b_blue, linewidth=0.8)
p <- p + annotate("segment", x=6, xend=6, y=11.4, yend=11.25, arrow=arrow(length=unit(0.08,"cm")), color=b_blue, linewidth=0.8)
p <- p + annotate("segment", x=4.5, xend=4.5, y=7.1, yend=6.55, arrow=arrow(length=unit(0.12,"cm")), color="white", linewidth=2)

# Side connections to middle
p <- p + annotate("segment", x=-1.3, xend=-2.5, y=10.6, yend=10.6, arrow=arrow(length=unit(0.06,"cm")), color=m_purp, linewidth=0.5, linetype="dotted")
p <- p + annotate("segment", x=1.3, xend=2.5, y=10.6, yend=10.6, arrow=arrow(length=unit(0.06,"cm")), color=m_purp, linewidth=0.5, linetype="dotted")
p <- p + annotate("segment", x=-1.3, xend=-2.5, y=8.95, yend=8.95, arrow=arrow(length=unit(0.06,"cm")), color=m_purp, linewidth=0.5, linetype="dotted")
p <- p + annotate("segment", x=1.3, xend=2.5, y=8.95, yend=8.95, arrow=arrow(length=unit(0.06,"cm")), color=m_purp, linewidth=0.5, linetype="dotted")

# Middle cascade: cross-disease -> module
p <- p + annotate("segment", x=0, xend=0, y=9.9, yend=9.75, arrow=arrow(length=unit(0.08,"cm")), color=m_purp, linewidth=0.8)

# =============== LEGEND ===============

p <- p + annotate("rect", xmin=-8, xmax=-5.5, ymin=2.0, ymax=2.8, fill=h_light, color=h_red, linewidth=0.5)
p <- p + annotate("text", x=-6.75, y=2.4, label="Upstream\nRegulator", size=2.8, lineheight=0.9, color=h_dark)

p <- p + annotate("rect", xmin=-4, xmax=-1.5, ymin=2.0, ymax=2.8, fill=h_mid, color=h_red, linewidth=0.5)
p <- p + annotate("text", x=-2.75, y=2.4, label="Translation\nProgram", size=2.8, lineheight=0.9, color="#4A0000")

p <- p + annotate("rect", xmin=0, xmax=3, ymin=2.0, ymax=2.8, fill=m_light, color=m_purp, linewidth=0.5)
p <- p + annotate("text", x=1.5, y=2.4, label="Cross-Disease\nIntegration", size=2.8, lineheight=0.9, color=m_purp)

p <- p + annotate("rect", xmin=4.5, xmax=7, ymin=2.0, ymax=2.8, fill=h_red, color=h_dark, linewidth=0.5)
p <- p + annotate("text", x=5.75, y=2.4, label="Physiological\nOutput", size=2.8, lineheight=0.9, color="white")

# =============== LAYOUT ===============
p <- p + coord_cartesian(xlim=c(-9,9), ylim=c(1.5, 17))
p <- p + theme_void()
p <- p + theme(plot.background=element_rect(fill="#FAFAFA", color=NA), plot.margin=margin(15,20,15,20))

ggsave(file.path(PROJ_DIR, "figures", "Figure6_Mechanistic_Model.png"), p, width=22, height=16, dpi=300)
cat("Figure 6 v3 saved\n")
