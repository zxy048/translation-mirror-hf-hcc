library(ggplot2)

PROJ_DIR <- "D:/R_projects/revision_analysis"

p <- ggplot()

# Title
p <- p + annotate("text", x=0.50, y=0.96, label="Cross-Disease Transcriptomic Analysis of Translation Regulation", fontface="bold", size=7, color="grey20")
p <- p + annotate("text", x=0.50, y=0.92, label="Data  ->  Analysis  ->  Findings  ->  Conclusion", size=4.5, color="grey50")

# === ROW 1: DATA (y: 0.76-0.88) ===
p <- p + annotate("rect", xmin=0.02, xmax=0.98, ymin=0.74, ymax=0.89, fill="#E8EAF6", alpha=0.4, color=NA)
p <- p + annotate("text", x=0.04, y=0.875, label="DATA", fontface="bold", size=5, color="#5C6BC0", hjust=0)

dat <- data.frame(x1=c(0.04,0.38,0.70), x2=c(0.34,0.66,0.98),
                  cx=c(0.19,0.52,0.84),
                  t=c("GSE57338","TCGA-LIHC","GSE141198"),
                  d=c("HF n=313\nDCM/ICM/NF\nAffymetrix","HCC n=424\n371T+50N\nRNA-seq","HCC n=148\n94 events\nRNA-seq"),
                  stringsAsFactors=FALSE)
for(i in 1:3) {
  p <- p + annotate("rect", xmin=dat$x1[i], xmax=dat$x2[i], ymin=0.76, ymax=0.87, fill="white", color="#5C6BC0", linewidth=1.2)
  p <- p + annotate("text", x=dat$cx[i], y=0.85, label=dat$t[i], fontface="bold", size=4.5, color="#283593")
  p <- p + annotate("text", x=dat$cx[i], y=0.795, label=dat$d[i], size=3, lineheight=0.9, color="grey40")
}

# === ROW 2: ANALYSIS (y: 0.48-0.70) ===
p <- p + annotate("rect", xmin=0.02, xmax=0.98, ymin=0.46, ymax=0.72, fill="#FFF3E0", alpha=0.4, color=NA)
p <- p + annotate("text", x=0.04, y=0.705, label="ANALYSIS", fontface="bold", size=5, color="#FF9800", hjust=0)

ana <- data.frame(
  x1=c(0.04,0.28,0.52,0.76), x2=c(0.24,0.48,0.72,0.96),
  cx=c(0.14,0.38,0.62,0.86),
  t=c("WGCNA","ssGSEA","Survival\nValidation","Upstream\nTF"),
  d=c("Signed network\nHub genes\nGO enrichment","83 pathways\nHallmark+KEGG+Reactome\nCross-disease ES","TGS (6-gene score)\nKM + Cox\n3 independent cohorts","19 TF candidates\nTF-TGS correlation\nFisher enrichment test"),
  src=c("GSE57338 | TCGA | GSE141198","TCGA | GSE57338","TCGA | GSE141198\nGSE14520 | GSE76427","TCGA-LIHC"),
  stringsAsFactors=FALSE)
for(i in 1:4) {
  p <- p + annotate("rect", xmin=ana$x1[i], xmax=ana$x2[i], ymin=0.485, ymax=0.70, fill="white", color="#FF9800", linewidth=1.2)
  p <- p + annotate("text", x=ana$cx[i], y=0.67, label=ana$t[i], fontface="bold", size=4.2, color="#E65100", lineheight=0.9)
  p <- p + annotate("text", x=ana$cx[i], y=0.58, label=ana$d[i], size=2.6, color="grey40", lineheight=0.9)
  p <- p + annotate("text", x=ana$cx[i], y=0.51, label=ana$src[i], size=2.3, color="#BF360C", lineheight=0.9, fontface="italic")
}

# === ROW 3: FINDINGS (y: 0.26-0.42) ===
p <- p + annotate("rect", xmin=0.02, xmax=0.98, ymin=0.24, ymax=0.44, fill="#E8F5E9", alpha=0.4, color=NA)
p <- p + annotate("text", x=0.04, y=0.425, label="FINDINGS", fontface="bold", size=5, color="#43A047", hjust=0)

fn <- data.frame(
  cx=c(0.14,0.38,0.62,0.86),
  ico=c("V","AV","X","S"),
  t=c("Module Replicated","Mirror Effect","Prognostic Specificity","ATF4/ISR+MYC"),
  d=c("Translation module replicated\nin GSE141198 (blue module)\np = 6.3 x 10^-12\n4/7 hub genes co-localized","Cross-disease pathway ES\nnegatively correlated\nrho = -0.598, p = 0.0003\nHCC up, HF down","TGS NOT prognostic in\n3 independent cohorts\nCox p = 0.479 (GSE141198)\nCohort-specific marker","ATF4 strongest TF\nrho=+0.439, FDR<0.0001\nMYC pathway rho=+0.613\nISR/UPR axis implicated"),
  col=c("#2E7D32","#7B1FA2","#F57C00","#1565C0"),
  stringsAsFactors=FALSE)
for(i in 1:4) {
  p <- p + annotate("rect", xmin=fn$cx[i]-0.10, xmax=fn$cx[i]+0.10, ymin=0.26, ymax=0.42, fill="white", color="#43A047", linewidth=1.2)
  p <- p + annotate("text", x=fn$cx[i], y=0.395, label=fn$ico[i], size=7, color=fn$col[i], fontface="bold")
  p <- p + annotate("text", x=fn$cx[i], y=0.36, label=fn$t[i], fontface="bold", size=3.5, color=fn$col[i], lineheight=0.9)
  p <- p + annotate("text", x=fn$cx[i], y=0.29, label=fn$d[i], size=2.4, color="grey40", lineheight=0.9)
}

# === ROW 4: CONCLUSION (y: 0.08-0.18) ===
p <- p + annotate("rect", xmin=0.02, xmax=0.98, ymin=0.06, ymax=0.20, fill="#FFEBEE", alpha=0.6, color="#E53935", linewidth=1.3)
p <- p + annotate("text", x=0.50, y=0.16, label="Translation co-expression is conserved across HF and HCC", fontface="bold", size=5, color="#B71C1C")
p <- p + annotate("text", x=0.50, y=0.11, label="Perturbation direction is disease-specific, driven by distinct upstream regulators (ATF4/ISR + MYC)", size=3.8, color="#C62828")

# === ARROWS: Analysis -> Findings ===
for(i in 1:4) {
  p <- p + annotate("segment", x=ana$cx[i], xend=fn$cx[i],
                    y=0.48, yend=0.425, arrow=arrow(length=unit(0.03,"npc")),
                    color="#43A047", linewidth=0.7)
}

p <- p + coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  theme_void() +
  theme(plot.background=element_rect(fill="white", color=NA))

ggsave(file.path(PROJ_DIR, "figures", "Figure1_Study_Design.png"), p, width=20, height=14, dpi=300)
cat("Figure 1 clean saved\n")
