library(ggplot2)

PROJ_DIR <- "D:/R_projects/revision_analysis"

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

p <- p + annotate("text", x=0, y=17.5, label="Mirror Regulation of Translation in HCC and Heart Failure", fontface="bold", size=8, color="grey20")
p <- p + annotate("text", x=0, y=16.5, label="Shared functional architecture  |  Opposed perturbation direction  |  Distinct upstream regulators", size=4.2, color="grey50")

# ===== DISEASE HEADERS =====
p <- p + annotate("text", x=-6, y=15.2, label="HCC", fontface="bold", size=12, color=h_red)
p <- p + annotate("text", x=6, y=15.2, label="HF", fontface="bold", size=12, color=b_blue)
p <- p + annotate("text", x=0, y=15.2, label="vs", size=5, color="grey60")
p <- p + annotate("text", x=-6, y=14.2, label="Hepatocellular Carcinoma  (TCGA-LIHC)", size=4.5, color=h_red)
p <- p + annotate("text", x=6, y=14.2, label="Dilated / Ischemic Cardiomyopathy  (GSE57338)", size=4.5, color=b_blue)

# Divider
p <- p + annotate("segment", x=0, xend=0, y=1.5, yend=13.5, linetype="dashed", color="grey70", linewidth=0.8)
p <- p + annotate("text", x=0, y=13.8, label="M I R R O R", size=4, color="grey60", fontface="italic")

# ==========================================
# HCC LEFT SIDE (3 boxes, well-spaced)
# ==========================================

# 1. UPSTREAM (y: 11.0-13.0)
p <- p + annotate("rect", xmin=-8.5, xmax=-3.5, ymin=10.8, ymax=13.2, fill=h_light, color=h_red, linewidth=1.2)
p <- p + annotate("text", x=-6, y=12.7, label="ATF4 / ISR", fontface="bold", size=5.5, color=h_dark)
p <- p + annotate("text", x=-6, y=12.05, label="TF rho = +0.439  FDR < 0.0001", size=3.5, color=h_red)
p <- p + annotate("text", x=-6, y=11.55, label="MYC Program", fontface="bold", size=5.5, color=h_dark)
p <- p + annotate("text", x=-6, y=11.0, label="Pathway rho = +0.613  p < 0.0001", size=3.5, color=h_red)

# 2. TRANSLATION (y: 7.0-10.2)
p <- p + annotate("rect", xmin=-8.5, xmax=-3.5, ymin=6.8, ymax=10.4, fill=h_mid, color=h_red, linewidth=1.3)
p <- p + annotate("text", x=-6, y=9.8, label="Ribosome & Translation", fontface="bold", size=5.5, color="#4A0000")
p <- p + annotate("text", x=-6, y=9.05, label="Translation Program Activated", size=4, color="#4A0000")
p <- p + annotate("text", x=-6, y=8.3, label="EEF1A1     FAU     RPL39     RPL3", size=3.5, color="#4A0000", fontface="italic")
p <- p + annotate("text", x=-6, y=7.7, label="RPL32     RPL41     RPS28", size=3.5, color="#4A0000", fontface="italic")
p <- p + annotate("text", x=-6, y=7.15, label="ATF4/ISR + MYC drive translational output", size=2.8, color="#BF360C")

# 3. OUTPUT (y: 4.0-6.0)
p <- p + annotate("rect", xmin=-8, xmax=-4, ymin=3.8, ymax=6.2, fill=h_red, color=h_dark, linewidth=1.5)
p <- p + annotate("text", x=-6, y=5.0, label="TRANSLATION", fontface="bold", size=6, color="white")
p <- p + annotate("text", x=-6, y=4.3, label="UP", fontface="bold", size=8, color="white")

p <- p + annotate("text", x=-6, y=3.0, label="Pro-survival ISR + MYC proliferation  ->  HCC growth", size=3.3, color=h_dark, fontface="italic")

# ==========================================
# HF RIGHT SIDE (3 boxes, mirror)
# ==========================================

# 1. UPSTREAM (y: 11.0-13.0)
p <- p + annotate("rect", xmin=3.5, xmax=8.5, ymin=10.8, ymax=13.2, fill=b_light, color=b_blue, linewidth=1.2)
p <- p + annotate("text", x=6, y=12.7, label="Energy Depletion", fontface="bold", size=5.5, color=b_dark)
p <- p + annotate("text", x=6, y=12.05, label="Metabolic / Oxidative Stress", size=3.5, color=b_blue)
p <- p + annotate("text", x=6, y=11.55, label="mTORC1  Down", fontface="bold", size=5.5, color=b_dark)
p <- p + annotate("text", x=6, y=11.0, label="MTOR rho = -0.391  FDR < 0.0001", size=3.5, color=b_blue)

# 2. TRANSLATION (y: 7.0-10.2)
p <- p + annotate("rect", xmin=3.5, xmax=8.5, ymin=6.8, ymax=10.4, fill=b_mid, color=b_blue, linewidth=1.3)
p <- p + annotate("text", x=6, y=9.8, label="Ribosome & Translation", fontface="bold", size=5.5, color="#002040")
p <- p + annotate("text", x=6, y=9.05, label="Translation Program Suppressed", size=4, color="#002040")
p <- p + annotate("text", x=6, y=8.3, label="7 Hub Genes Direction Discordant", size=3.5, color="#002040", fontface="italic")
p <- p + annotate("text", x=6, y=7.7, label="ssGSEA: Translation pathways down vs normal", size=3.5, color="#002040")
p <- p + annotate("text", x=6, y=7.15, label="ATP-costly protein synthesis suppressed", size=2.8, color=b_dark)

# 3. OUTPUT (y: 4.0-6.0)
p <- p + annotate("rect", xmin=4, xmax=8, ymin=3.8, ymax=6.2, fill=b_blue, color=b_dark, linewidth=1.5)
p <- p + annotate("text", x=6, y=5.0, label="TRANSLATION", fontface="bold", size=6, color="white")
p <- p + annotate("text", x=6, y=4.3, label="DOWN", fontface="bold", size=8, color="white")

p <- p + annotate("text", x=6, y=3.0, label="Energy conservation  ->  HF adaptation", size=3.3, color=b_dark, fontface="italic")

# ==========================================
# MIDDLE INTEGRATION (2 panels, centered)
# ==========================================

p <- p + annotate("rect", xmin=-2.8, xmax=2.8, ymin=10.5, ymax=13.0, fill=m_light, color=m_purp, linewidth=1.3)
p <- p + annotate("text", x=0, y=12.4, label="Cross-Disease", fontface="bold", size=5, color=m_purp)
p <- p + annotate("text", x=0, y=11.7, label="Pathway Effect Size", size=3.5, color=m_purp)
p <- p + annotate("text", x=0, y=11.3, label="NEGATIVE CORRELATION", fontface="bold", size=3.5, color="#4A148C")
p <- p + annotate("text", x=0, y=10.8, label="rho = -0.598   p = 0.0003", fontface="bold", size=3.2, color="#4A148C")

p <- p + annotate("rect", xmin=-2.8, xmax=2.8, ymin=8.0, ymax=10.2, fill="#CE93D8", color=m_purp, linewidth=1.3)
p <- p + annotate("text", x=0, y=9.7, label="Functional Architecture", fontface="bold", size=5, color="#4A148C")
p <- p + annotate("text", x=0, y=9.05, label="CONSERVED", fontface="bold", size=4, color="#4A148C")
p <- p + annotate("text", x=0, y=8.55, label="Network Topology", size=3.2, color="#4A148C")
p <- p + annotate("text", x=0, y=8.15, label="CONTEXT-DEPENDENT REWIRING", fontface="bold", size=3, color="#4A148C")

# ==========================================
# ARROWS
# ==========================================

# HCC vertical arrows
p <- p + annotate("segment", x=-6, xend=-6, y=10.7, yend=10.5, arrow=arrow(length=unit(0.1,"cm")), color=h_red, linewidth=1.3)
p <- p + annotate("segment", x=-6, xend=-6, y=6.7, yend=6.3, arrow=arrow(length=unit(0.15,"cm")), color="white", linewidth=2.2)

# HF vertical arrows
p <- p + annotate("segment", x=6, xend=6, y=10.7, yend=10.5, arrow=arrow(length=unit(0.1,"cm")), color=b_blue, linewidth=1.3)
p <- p + annotate("segment", x=6, xend=6, y=6.7, yend=6.3, arrow=arrow(length=unit(0.15,"cm")), color="white", linewidth=2.2)

# Lateral connections to middle
p <- p + annotate("segment", x=-3.5, xend=-2.85, y=11.9, yend=11.9, arrow=arrow(length=unit(0.06,"cm")), color=m_purp, linewidth=0.5, linetype="dotted")
p <- p + annotate("segment", x=3.5, xend=2.85, y=11.9, yend=11.9, arrow=arrow(length=unit(0.06,"cm")), color=m_purp, linewidth=0.5, linetype="dotted")
p <- p + annotate("segment", x=-3.5, xend=-2.85, y=9.1, yend=9.1, arrow=arrow(length=unit(0.06,"cm")), color=m_purp, linewidth=0.5, linetype="dotted")
p <- p + annotate("segment", x=3.5, xend=2.85, y=9.1, yend=9.1, arrow=arrow(length=unit(0.06,"cm")), color=m_purp, linewidth=0.5, linetype="dotted")

# Middle cascade
p <- p + annotate("segment", x=0, xend=0, y=10.4, yend=10.25, arrow=arrow(length=unit(0.08,"cm")), color=m_purp, linewidth=0.9)

# ==========================================
# BOTTOM ANNOTATION
# ==========================================

p <- p + annotate("text", x=0, y=2.9, label="Mirror regulation observed at pathway and network level, not at individual gene level",
                  size=3.3, color="grey50", fontface="italic")

# ==========================================
# LEGEND
# ==========================================

p <- p + annotate("rect", xmin=-8.5, xmax=-6, ymin=1.7, ymax=2.5, fill=h_light, color=h_red, linewidth=0.5)
p <- p + annotate("text", x=-7.25, y=2.1, label="Upstream\nRegulator", size=3, lineheight=0.9, color=h_dark)

p <- p + annotate("rect", xmin=-4, xmax=-1.5, ymin=1.7, ymax=2.5, fill=h_mid, color=h_red, linewidth=0.5)
p <- p + annotate("text", x=-2.75, y=2.1, label="Translation\nProgram", size=3, lineheight=0.9, color="#4A0000")

p <- p + annotate("rect", xmin=0.5, xmax=3, ymin=1.7, ymax=2.5, fill=m_light, color=m_purp, linewidth=0.5)
p <- p + annotate("text", x=1.75, y=2.1, label="Integration\nLayer", size=3, lineheight=0.9, color=m_purp)

p <- p + annotate("rect", xmin=4.5, xmax=7, ymin=1.7, ymax=2.5, fill=h_red, color=h_dark, linewidth=0.5)
p <- p + annotate("text", x=5.75, y=2.1, label="Physiological\nOutput", size=3, lineheight=0.9, color="white")

# ==========================================
p <- p + coord_cartesian(xlim=c(-10,10), ylim=c(1, 18))
p <- p + theme_void()
p <- p + theme(plot.background=element_rect(fill="#FAFAFA", color=NA), plot.margin=margin(20,25,20,25))

ggsave(file.path(PROJ_DIR, "figures", "Figure6_Mechanistic_Model.png"), p, width=22, height=16, dpi=300)
cat("Figure 6 v4 saved\n")
