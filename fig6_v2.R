library(ggplot2)

PROJ_DIR <- "D:/R_projects/revision_analysis"

# Colors
HCC_red   <- "#C62828"
HCC_light <- "#FFCDD2"
HCC_mid   <- "#EF9A9A"
HCC_dark  <- "#B71C1C"
HF_blue   <- "#1565C0"
HF_light  <- "#BBDEFB"
HF_mid    <- "#90CAF9"
HF_dark   <- "#0D47A1"
MID_purp  <- "#7B1FA2"
MID_light <- "#E1BEE7"

p <- ggplot()

# ===== TITLE =====
p <- p + annotate("text", x=0, y=14.5, label="Mirror Regulation of Translation in HCC and Heart Failure", fontface="bold", size=7.5, color="grey20")
p <- p + annotate("text", x=0, y=13.8, label="Shared co-expression architecture, opposed perturbation direction, distinct upstream regulators", size=4.2, color="grey50")

# ===== DISEASE LABELS =====
p <- p + annotate("text", x=-4, y=13.0, label="HCC", fontface="bold", size=9, color=HCC_red)
p <- p + annotate("text", x=4, y=13.0, label="HF", fontface="bold", size=9, color=HF_blue)
p <- p + annotate("text", x=0, y=13.0, label="vs", size=4, color="grey60")
p <- p + annotate("text", x=-4, y=12.3, label="Hepatocellular Carcinoma", size=3.8, color=HCC_red)
p <- p + annotate("text", x=4, y=12.3, label="Dilated / Ischemic Cardiomyopathy", size=3.8, color=HF_blue)

# ===== CENTER DIVIDER =====
p <- p + annotate("segment", x=0, xend=0, y=1.5, yend=11.8, linetype="dashed", color="grey70", linewidth=0.8)
p <- p + annotate("text", x=0, y=12.0, label="M I R R O R", size=3.5, color="grey60", fontface="italic")

# ===== HCC SIDE (LEFT) =====

# Upstream box
p <- p + annotate("rect", xmin=-7.0, xmax=-1.0, ymin=10.0, ymax=11.2, fill=HCC_light, color=HCC_red, linewidth=1.0)
p <- p + annotate("text", x=-5.5, y=10.65, label="ATF4 / ISR", fontface="bold", size=4.5, color=HCC_dark)
p <- p + annotate("text", x=-5.5, y=10.25, label="TF rho = +0.439  FDR < 0.0001", size=3.2, color=HCC_red)
p <- p + annotate("text", x=-2.5, y=10.65, label="MYC Program", fontface="bold", size=4.5, color=HCC_dark)
p <- p + annotate("text", x=-2.5, y=10.25, label="Pathway rho = +0.613  p < 0.0001", size=3.2, color=HCC_red)

# Downstream stress mediators
p <- p + annotate("rect", xmin=-5.5, xmax=-2.5, ymin=8.8, ymax=9.8, fill="#FFE0B2", color="#E65100", linewidth=0.8)
p <- p + annotate("text", x=-4, y=9.3, label="DDIT3 (CHOP)", fontface="bold", size=3.5, color="#BF360C")
p <- p + annotate("text", x=-4, y=8.95, label="rho = +0.289  FDR < 0.0001", size=2.8, color="#E65100")

p <- p + annotate("text", x=-4, y=8.45, label="XBP1 rho=+0.162  HIF1A rho=-0.417", size=2.8, color="#E65100", fontface="italic")

# Ribosome box
p <- p + annotate("rect", xmin=-6.5, xmax=-1.5, ymin=6.2, ymax=8.0, fill=HCC_mid, color=HCC_red, linewidth=1.2)
p <- p + annotate("text", x=-4, y=7.5, label="Ribosome Biogenesis", fontface="bold", size=4.5, color="#4A0000")
p <- p + annotate("text", x=-4, y=7.0, label="Translation Initiation & Elongation", size=3.5, color="#4A0000")
p <- p + annotate("text", x=-4, y=6.5, label="EEF1A1  FAU  RPL39  RPL3  RPL32  RPL41  RPS28", size=3.0, color="#4A0000", fontface="italic")

# Output box
p <- p + annotate("rect", xmin=-6.0, xmax=-2.0, ymin=4.5, ymax=5.8, fill=HCC_red, color=HCC_dark, linewidth=1.5)
p <- p + annotate("text", x=-4, y=5.15, label="TRANSLATION  UP", fontface="bold", size=5.5, color="white")

# Annotation
p <- p + annotate("text", x=-4, y=3.8, label="Pro-survival ISR adaptation", size=3.0, color=HCC_dark, fontface="italic")
p <- p + annotate("text", x=-4, y=3.3, label="+ MYC-driven proliferation", size=3.0, color=HCC_dark, fontface="italic")

# ===== HF SIDE (RIGHT) =====

# Upstream box
p <- p + annotate("rect", xmin=1.0, xmax=7.0, ymin=10.0, ymax=11.2, fill=HF_light, color=HF_blue, linewidth=1.0)
p <- p + annotate("text", x=5.5, y=10.65, label="Energy Depletion", fontface="bold", size=4.5, color=HF_dark)
p <- p + annotate("text", x=5.5, y=10.25, label="Metabolic Stress", size=3.2, color=HF_blue)
p <- p + annotate("text", x=2.5, y=10.65, label="mTORC1  Down", fontface="bold", size=4.5, color=HF_dark)
p <- p + annotate("text", x=2.5, y=10.25, label="MTOR rho = -0.391  FDR < 0.0001", size=3.2, color=HF_blue)

# Adaptation
p <- p + annotate("rect", xmin=2.5, xmax=5.5, ymin=8.8, ymax=9.8, fill="#E3F2FD", color=HF_dark, linewidth=0.8)
p <- p + annotate("text", x=4, y=9.3, label="ATP-costly protein", fontface="bold", size=3.5, color=HF_dark)
p <- p + annotate("text", x=4, y=8.95, label="synthesis suppressed", size=3.2, color=HF_dark)

# Ribosome box
p <- p + annotate("rect", xmin=1.5, xmax=6.5, ymin=6.2, ymax=8.0, fill=HF_mid, color=HF_blue, linewidth=1.2)
p <- p + annotate("text", x=4, y=7.5, label="Ribosome & Translation", fontface="bold", size=4.5, color="#002040")
p <- p + annotate("text", x=4, y=7.0, label="Genes Downregulated", size=3.5, color="#002040")
p <- p + annotate("text", x=4, y=6.5, label="7 Hub Genes Direction Discordant", size=3.0, color="#002040", fontface="italic")

# Output box
p <- p + annotate("rect", xmin=2.0, xmax=6.0, ymin=4.5, ymax=5.8, fill=HF_blue, color=HF_dark, linewidth=1.5)
p <- p + annotate("text", x=4, y=5.15, label="TRANSLATION  DOWN", fontface="bold", size=5.5, color="white")

# Annotation
p <- p + annotate("text", x=4, y=3.8, label="Energy conservation strategy", size=3.0, color=HF_dark, fontface="italic")
p <- p + annotate("text", x=4, y=3.3, label="Cardiomyocyte survival adaptation", size=3.0, color=HF_dark, fontface="italic")

# ===== SHARED MIDDLE PANELS =====

p <- p + annotate("rect", xmin=-2.0, xmax=2.0, ymin=8.8, ymax=10.2, fill=MID_light, color=MID_purp, linewidth=1.2)
p <- p + annotate("text", x=0, y=9.7, label="Cross-Disease", fontface="bold", size=3.8, color=MID_purp)
p <- p + annotate("text", x=0, y=9.25, label="Pathway Effect Sizes", size=3.3, color=MID_purp)
p <- p + annotate("text", x=0, y=8.95, label="rho = -0.598, p = 0.0003", fontface="bold", size=3.0, color="#4A148C")

p <- p + annotate("rect", xmin=-2.0, xmax=2.0, ymin=7.2, ymax=8.4, fill="#CE93D8", color=MID_purp, linewidth=1.2)
p <- p + annotate("text", x=0, y=8.0, label="Co-expression Module", fontface="bold", size=3.8, color="#4A148C")
p <- p + annotate("text", x=0, y=7.6, label="Conserved Across Diseases", size=3.3, color="#4A148C")
p <- p + annotate("text", x=0, y=7.3, label="WGCNA p = 6.3e-12", fontface="bold", size=2.8, color="#4A148C")

# ===== ARROWS =====

# HCC arrows
p <- p + annotate("segment", x=-5.5, xend=-5.5, y=9.95, yend=9.35,
                  arrow=arrow(length=unit(0.1,"cm")), color=HCC_red, linewidth=1.0)
p <- p + annotate("segment", x=-2.5, xend=-2.5, y=9.95, yend=9.35,
                  arrow=arrow(length=unit(0.1,"cm")), color=HCC_red, linewidth=1.0)
p <- p + annotate("segment", x=-4, xend=-4, y=8.74, yend=8.05,
                  arrow=arrow(length=unit(0.1,"cm")), color=HCC_red, linewidth=1.2)
p <- p + annotate("segment", x=-4, xend=-4, y=6.14, yend=5.85,
                  arrow=arrow(length=unit(0.12,"cm")), color="white", linewidth=1.8)

# HF arrows
p <- p + annotate("segment", x=5.5, xend=5.5, y=9.95, yend=9.35,
                  arrow=arrow(length=unit(0.1,"cm")), color=HF_blue, linewidth=1.0)
p <- p + annotate("segment", x=2.5, xend=2.5, y=9.95, yend=9.35,
                  arrow=arrow(length=unit(0.1,"cm")), color=HF_blue, linewidth=1.0)
p <- p + annotate("segment", x=4, xend=4, y=6.14, yend=5.85,
                  arrow=arrow(length=unit(0.12,"cm")), color="white", linewidth=1.8)

# Lateral connections to middle
p <- p + annotate("segment", x=-1.5, xend=-2.0, y=9.5, yend=9.5,
                  arrow=arrow(length=unit(0.06,"cm")), color=MID_purp, linewidth=0.5, linetype="dotted")
p <- p + annotate("segment", x=1.5, xend=2.0, y=9.5, yend=9.5,
                  arrow=arrow(length=unit(0.06,"cm")), color=MID_purp, linewidth=0.5, linetype="dotted")

# Middle downward arrow
p <- p + annotate("segment", x=0, xend=0, y=8.74, yend=8.45,
                  arrow=arrow(length=unit(0.1,"cm")), color=MID_purp, linewidth=1.0)

# ===== LEGEND =====
p <- p + annotate("rect", xmin=-7, xmax=-4.5, ymin=1.6, ymax=2.3, fill=HCC_light, color=HCC_red, linewidth=0.5)
p <- p + annotate("text", x=-5.75, y=1.95, label="Upstream\nRegulator", size=2.5, lineheight=0.9, color=HCC_dark)

p <- p + annotate("rect", xmin=-3.5, xmax=-1, ymin=1.6, ymax=2.3, fill=HCC_mid, color=HCC_red, linewidth=0.5)
p <- p + annotate("text", x=-2.25, y=1.95, label="Translation\nProgram", size=2.5, lineheight=0.9, color="#4A0000")

p <- p + annotate("rect", xmin=0, xmax=2.5, ymin=1.6, ymax=2.3, fill=MID_light, color=MID_purp, linewidth=0.5)
p <- p + annotate("text", x=1.25, y=1.95, label="Integration\nLayer", size=2.5, lineheight=0.9, color=MID_purp)

p <- p + annotate("rect", xmin=3.5, xmax=6, ymin=1.6, ymax=2.3, fill=HCC_red, color=HCC_dark, linewidth=0.5)
p <- p + annotate("text", x=4.75, y=1.95, label="Physiological\nOutput", size=2.5, lineheight=0.9, color="white")

# ===== LAYOUT =====
p <- p + coord_cartesian(xlim=c(-8,8), ylim=c(1, 15.2))
p <- p + theme_void()
p <- p + theme(plot.background=element_rect(fill="#FAFAFA", color=NA), plot.margin=margin(15,15,15,15))

ggsave(file.path(PROJ_DIR, "figures", "Figure6_Mechanistic_Model.png"), p, width=20, height=14, dpi=300)
cat("Figure 6 v2 saved\n")
