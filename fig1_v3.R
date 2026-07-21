library(ggplot2)

PROJ_DIR <- "D:/R_projects/revision_analysis"

Cdata  <- "#5C6BC0"
Cana   <- "#FF9800"
Cfind  <- "#43A047"
Cconc  <- "#E53935"

p <- ggplot()

# ===== TITLE =====
p <- p + annotate("text", x=6, y=14.2, label="Cross-Disease Transcriptomic Analysis of Translation Regulation", fontface="bold", size=6.5, color="grey20")
p <- p + annotate("text", x=6, y=13.5, label="From shared co-expression signal to disease-specific mirror regulation mechanism", size=3.8, color="grey50")

# ===== DATA ROW (y: 10.8–12.2) =====
p <- p + annotate("rect", xmin=0.3, xmax=11.7, ymin=10.5, ymax=12.3, fill="#E8EAF6", alpha=0.4, color=NA)
p <- p + annotate("text", x=0.6, y=12.05, label="DATA", fontface="bold", size=5, color=Cdata, hjust=0)

data_boxes <- data.frame(
  x1=c(0.8, 4.3, 7.8),
  x2=c(3.6, 7.1, 10.6),
  cx=c(2.2, 5.7, 9.2),
  title=c("GSE57338", "TCGA-LIHC", "GSE141198"),
  desc=c("Heart Failure\nn = 313\nDCM / ICM / NF\nAffymetrix 1.1 ST",
         "HCC Discovery\nn = 424\n371 T + 50 N\nRNA-seq Illumina",
         "HCC Validation\nn = 148\n94 OS events\nRNA-seq"),
  stringsAsFactors=FALSE
)

for(i in 1:3) {
  p <- p + annotate("rect", xmin=data_boxes$x1[i], xmax=data_boxes$x2[i],
                    ymin=10.7, ymax=12.05, fill="white", color=Cdata, linewidth=1.3)
  p <- p + annotate("text", x=data_boxes$cx[i], y=11.7,
                    label=data_boxes$title[i], fontface="bold", size=4.5, color="#283593")
  p <- p + annotate("text", x=data_boxes$cx[i], y=11.05,
                    label=data_boxes$desc[i], size=2.8, lineheight=0.9, color="grey40")
}

# ===== ANALYSIS ROW (y: 6.5–9.2) =====
p <- p + annotate("rect", xmin=0.3, xmax=11.7, ymin=6.2, ymax=9.4, fill="#FFF3E0", alpha=0.4, color=NA)
p <- p + annotate("text", x=0.6, y=9.15, label="ANALYSIS", fontface="bold", size=5, color=Cana, hjust=0)

ana <- data.frame(
  x1=c(0.8, 3.3, 5.8, 8.3),
  x2=c(2.8, 5.3, 7.8, 10.3),
  cx=c(1.8, 4.3, 6.8, 9.3),
  title=c("WGCNA", "ssGSEA", "Survival\nValidation", "Upstream\nTF Prediction"),
  desc=c("Signed co-expression\nHub gene identification\nGO / KEGG enrichment",
         "81 pathways\nHallmark + KEGG + Reactome\nCross-disease effect sizes",
         "Translation Gene Score\nKM + Cox regression\n3 independent cohorts",
         "19 TF candidates\nTF-TGS Spearman correlation\nFisher enrichment test"),
  sources=c("GSE57338\nTCGA-LIHC\nGSE141198","TCGA-LIHC\nGSE57338","TCGA-LIHC\nGSE141198\nGSE14520\nGSE76427","TCGA-LIHC"),
  stringsAsFactors=FALSE
)

for(i in 1:4) {
  p <- p + annotate("rect", xmin=ana$x1[i], xmax=ana$x2[i],
                    ymin=6.5, ymax=9.1, fill="white", color=Cana, linewidth=1.3)
  p <- p + annotate("text", x=ana$cx[i], y=8.7,
                    label=ana$title[i], fontface="bold", size=4, color="#E65100", lineheight=0.9)
  p <- p + annotate("text", x=ana$cx[i], y=7.5,
                    label=ana$desc[i], size=2.6, color="grey40", lineheight=0.9)
  p <- p + annotate("text", x=ana$cx[i], y=6.75,
                    label=paste0("Data: ", ana$sources[i]), size=2.3, color="#BF360C", lineheight=0.9, fontface="italic")
}

# ===== FINDINGS ROW (y: 2.5–5.0) =====
p <- p + annotate("rect", xmin=0.3, xmax=11.7, ymin=2.2, ymax=5.2, fill="#E8F5E9", alpha=0.4, color=NA)
p <- p + annotate("text", x=0.6, y=4.95, label="FINDINGS", fontface="bold", size=5, color=Cfind, hjust=0)

fn <- data.frame(
  cx=c(1.8, 4.3, 6.8, 9.3),
  icon=c("network", "mirror", "prognosis", "mechanism"),
  ico_label=c("Module\nReplicated","Pathway\nMirror Effect","Prognostic\nSpecificity","ATF4/ISR\n+ MYC Axis"),
  desc=c("Translation co-expression module\nindependently replicated in GSE141198\nblue module p = 6.3 x 10^-12\n4/7 hub genes co-localized",
         "Cross-disease pathway effect sizes\nnegatively correlated: rho = -0.598\np = 0.0003; permutation p = 0.0079\nHCC up, HF down",
         "TGS NOT prognostic in 3 external\ncohorts (Cox p = 0.479 in GSE141198)\nPrognostic value is TCGA-LIHC specific",
         "ATF4 strongest TF (rho = +0.439, FDR<0.0001)\nMYC target pathway rho = +0.613\nISR/UPR axis implicated"),
  symbol=c("V", "AV", "X", "S"),
  stringsAsFactors=FALSE
)

for(i in 1:4) {
  p <- p + annotate("rect", xmin=fn$cx[i]-1.1, xmax=fn$cx[i]+1.1,
                    ymin=2.5, ymax=5.0, fill="white", color=Cfind, linewidth=1.3)
  # Icon as symbol
  sym_col <- ifelse(fn$symbol[i]=="V", "#2E7D32",
             ifelse(fn$symbol[i]=="AV", "#7B1FA2",
             ifelse(fn$symbol[i]=="X", "#F57C00",
             "#1565C0")))
  p <- p + annotate("text", x=fn$cx[i], y=4.55,
                    label=fn$symbol[i], size=8, color=sym_col, fontface="bold")
  p <- p + annotate("text", x=fn$cx[i], y=4.05,
                    label=fn$ico_label[i], fontface="bold", size=3.5,
                    color=sym_col, lineheight=0.9)
  p <- p + annotate("text", x=fn$cx[i], y=3.0,
                    label=fn$desc[i], size=2.4, color="grey40", lineheight=0.9)
}

# ===== ARROWS: Data -> Analysis (explicit paths) =====
# GSE57338(cx=2.2) -> WGCNA(cx=1.8), ssGSEA(cx=4.3)
p <- p + annotate("segment", x=2.2, xend=1.8, y=10.68, yend=9.15,
                  arrow=arrow(length=unit(0.08,"cm")), color="grey50", linewidth=0.6)
p <- p + annotate("segment", x=2.2, xend=4.3, y=10.68, yend=9.15,
                  arrow=arrow(length=unit(0.08,"cm")), color="grey50", linewidth=0.6)

# TCGA-LIHC(cx=5.7) -> WGCNA(1.8), ssGSEA(4.3), Survival(6.8), TF(9.3)
p <- p + annotate("segment", x=5.7, xend=1.8, y=10.68, yend=9.15,
                  arrow=arrow(length=unit(0.08,"cm")), color="grey50", linewidth=0.6)
p <- p + annotate("segment", x=5.7, xend=4.3, y=10.68, yend=9.15,
                  arrow=arrow(length=unit(0.08,"cm")), color="grey50", linewidth=0.6)
p <- p + annotate("segment", x=5.7, xend=6.8, y=10.68, yend=9.15,
                  arrow=arrow(length=unit(0.08,"cm")), color="grey50", linewidth=0.6)
p <- p + annotate("segment", x=5.7, xend=9.3, y=10.68, yend=9.15,
                  arrow=arrow(length=unit(0.08,"cm")), color="grey50", linewidth=0.6)

# GSE141198(cx=9.2) -> WGCNA(1.8), Survival(6.8)
p <- p + annotate("segment", x=9.2, xend=1.8, y=10.68, yend=9.15,
                  arrow=arrow(length=unit(0.08,"cm")), color="grey50", linewidth=0.6)
p <- p + annotate("segment", x=9.2, xend=6.8, y=10.68, yend=9.15,
                  arrow=arrow(length=unit(0.08,"cm")), color="grey50", linewidth=0.6)

# ===== ARROWS: Analysis -> Findings =====
p <- p + annotate("segment", x=1.8, xend=1.8, y=6.48, yend=5.05,
                  arrow=arrow(length=unit(0.08,"cm")), color=Cfind, linewidth=0.8)
p <- p + annotate("segment", x=4.3, xend=4.3, y=6.48, yend=5.05,
                  arrow=arrow(length=unit(0.08,"cm")), color=Cfind, linewidth=0.8)
p <- p + annotate("segment", x=6.8, xend=6.8, y=6.48, yend=5.05,
                  arrow=arrow(length=unit(0.08,"cm")), color=Cfind, linewidth=0.8)
p <- p + annotate("segment", x=9.3, xend=9.3, y=6.48, yend=5.05,
                  arrow=arrow(length=unit(0.08,"cm")), color=Cfind, linewidth=0.8)

# ===== Divider lines between rows =====
p <- p + annotate("segment", x=0.3, xend=11.7, y=10.0, yend=10.0, linetype="dotted", color="grey80", linewidth=0.3)
p <- p + annotate("segment", x=0.3, xend=11.7, y=5.8, yend=5.8, linetype="dotted", color="grey80", linewidth=0.3)

# ===== CONCLUSION BANNER =====
p <- p + annotate("rect", xmin=0.3, xmax=11.7, ymin=0.5, ymax=1.9, fill="#FFEBEE", alpha=0.6, color=Cconc, linewidth=1.3)
p <- p + annotate("text", x=6, y=1.4, label="Translation co-expression is conserved across HF and HCC", fontface="bold", size=4.5, color="#B71C1C")
p <- p + annotate("text", x=6, y=0.85, label="Perturbation direction is disease-type-specific, associated with distinct upstream regulators (ATF4/ISR + MYC)", size=3.5, color="#C62828")

# ===== Layout =====
p <- p + coord_cartesian(xlim=c(0,12), ylim=c(0,14.8))
p <- p + theme_void()
p <- p + theme(plot.background=element_rect(fill="#FAFAFA", color=NA), plot.margin=margin(15,15,15,15))

ggsave(file.path(PROJ_DIR, "figures", "Figure1_Study_Design.png"), p, width=16, height=12, dpi=300)
cat("Figure 1 v3 saved\n")
