library(ggplot2)

PROJ_DIR <- "D:/R_projects/revision_analysis"

col_data <- "#5C6BC0"
col_ana  <- "#FF9800"
col_find <- "#43A047"
col_conc <- "#E53935"

p <- ggplot()

# DATA row background
p <- p + annotate("rect", xmin=0.5, xmax=11.5, ymin=9.5, ymax=12.8, fill="#E8EAF6", alpha=0.5, color=NA)
p <- p + annotate("text", x=0.8, y=12.5, label="DATA", fontface="bold", size=5.5, color=col_data, hjust=0)

# Three data boxes
p <- p + annotate("rect", xmin=1.2, xmax=4.0, ymin=9.8, ymax=12.2, fill="white", color=col_data, linewidth=1.5)
p <- p + annotate("text", x=2.6, y=11.5, label="GSE57338", fontface="bold", size=4.5, color="#283593")
p <- p + annotate("text", x=2.6, y=10.8, label="HF (n=313)\nDCM + ICM + NF\nAffymetrix 1.1 ST", size=3.0, lineheight=0.9, color="grey40")

p <- p + annotate("rect", xmin=4.8, xmax=7.6, ymin=9.8, ymax=12.2, fill="white", color=col_data, linewidth=1.5)
p <- p + annotate("text", x=6.2, y=11.5, label="TCGA-LIHC", fontface="bold", size=4.5, color="#283593")
p <- p + annotate("text", x=6.2, y=10.8, label="HCC (n=424)\n371 Tumor + 50 Normal\nRNA-seq Illumina", size=3.0, lineheight=0.9, color="grey40")

p <- p + annotate("rect", xmin=8.4, xmax=11.2, ymin=9.8, ymax=12.2, fill="white", color=col_data, linewidth=1.5)
p <- p + annotate("text", x=9.8, y=11.5, label="GSE141198", fontface="bold", size=4.5, color="#283593")
p <- p + annotate("text", x=9.8, y=10.8, label="HCC Validation\nn=148 RNA-seq\n94 OS events", size=3.0, lineheight=0.9, color="grey40")

# Arrows data->analysis
p <- p + annotate("segment", x=4.0, xend=4.6, y=11.0, yend=8.5, arrow=arrow(length=unit(0.1,"cm")), color=col_ana, linewidth=0.8)
p <- p + annotate("segment", x=6.2, xend=6.2, y=9.75, yend=8.9, arrow=arrow(length=unit(0.1,"cm")), color=col_ana, linewidth=0.8)
p <- p + annotate("segment", x=8.0, xend=8.4, y=11.0, yend=8.5, arrow=arrow(length=unit(0.1,"cm")), color=col_ana, linewidth=0.8)

# ANALYSIS row background
p <- p + annotate("rect", xmin=0.5, xmax=11.5, ymin=5.8, ymax=8.8, fill="#FFF3E0", alpha=0.5, color=NA)
p <- p + annotate("text", x=0.8, y=8.55, label="ANALYSIS", fontface="bold", size=5.5, color=col_ana, hjust=0)

# Analysis boxes
an <- data.frame(
  x1=c(1.2,3.6,6.0,8.4), x2=c(3.2,5.6,8.0,10.4),
  cx=c(2.2,4.6,7.0,9.4),
  t=c("WGCNA","ssGSEA","Survival\nValidation","Upstream\nTF Prediction"),
  s=c("Module discovery\nHub gene ID\nGO enrichment","81 pathways\nHallmark+KEGG+Reactome\nCross-disease ES","TGS (6-gene score)\nKM + Cox regression\n3 independent cohorts","19 TF candidates\nTF-TGS correlation\nATF4/MYC/E2F/mTORC1"),
  stringsAsFactors=FALSE
)

for(i in 1:4) {
  p <- p + annotate("rect", xmin=an$x1[i], xmax=an$x2[i], ymin=6.1, ymax=8.6, fill="white", color=col_ana, linewidth=1.3)
  p <- p + annotate("text", x=an$cx[i], y=8.1, label=an$t[i], fontface="bold", size=3.8, color="#E65100", lineheight=0.9)
  p <- p + annotate("text", x=an$cx[i], y=6.8, label=an$s[i], size=2.5, color="grey40", lineheight=0.9)
  p <- p + annotate("segment", x=an$cx[i], xend=an$cx[i], y=6.05, yend=5.2, arrow=arrow(length=unit(0.1,"cm")), color=col_find, linewidth=0.8)
}

# FINDINGS row background
p <- p + annotate("rect", xmin=0.5, xmax=11.5, ymin=2.0, ymax=5.1, fill="#E8F5E9", alpha=0.5, color=NA)
p <- p + annotate("text", x=0.8, y=4.85, label="FINDINGS", fontface="bold", size=5.5, color=col_find, hjust=0)

fn <- data.frame(
  cx=c(2.2,4.6,7.0,9.4),
  ic=c("✓","!","✗","★"),
  t=c("Module\nConservation","Mirror\nEffect","Prognostic\nSpecificity","ATF4/ISR +\nMYC Axis"),
  d=c("Translation module\nreplicated in GSE141198\nblue module p=6.3e-12\n4/7 hub genes co-localized","Cross-disease pathway\nES rho=-0.598\np=0.0003\nHCC up vs HF down","TGS NOT prognostic\nin 3 external cohorts\nCox p=0.479 (GSE141198)\nCohort-specific","ATF4 strongest TF\nrho=+0.439, FDR<0.0001\nMYC pathway rho=+0.613\nISR/UPR axis implicated"),
  stringsAsFactors=FALSE
)

for(i in 1:4) {
  p <- p + annotate("rect", xmin=fn$cx[i]-1.2, xmax=fn$cx[i]+1.2, ymin=2.3, ymax=4.9, fill="white", color=col_find, linewidth=1.3)
  p <- p + annotate("text", x=fn$cx[i], y=4.4, label=fn$ic[i], size=10, color=col_find, fontface="bold")
  p <- p + annotate("text", x=fn$cx[i], y=3.7, label=fn$t[i], fontface="bold", size=3.5, color="#2E7D32", lineheight=0.9)
  p <- p + annotate("text", x=fn$cx[i], y=2.8, label=fn$d[i], size=2.5, color="grey40", lineheight=0.85)
}

# CONCLUSION banner
p <- p + annotate("rect", xmin=0.5, xmax=11.5, ymin=0.3, ymax=1.7, fill="#FFEBEE", alpha=0.7, color=col_conc, linewidth=1.5)
p <- p + annotate("text", x=6.0, y=1.2, label="Translation co-expression is conserved across HF and HCC, but perturbation direction is disease-type-specific", fontface="bold", size=4.2, color="#B71C1C")
p <- p + annotate("text", x=6.0, y=0.6, label="ATF4/ISR stress adaptation + MYC proliferative program jointly regulate HCC translational program", size=3.3, color="#C62828")

# Title
p <- p + annotate("text", x=6.0, y=13.3, label="Cross-Disease Transcriptomic Analysis of Translation Regulation", fontface="bold", size=6.5, color="grey20")

p <- p + coord_cartesian(xlim=c(0,12), ylim=c(0,13.8))
p <- p + theme_void()
p <- p + theme(plot.background=element_rect(fill="#FAFAFA", color=NA), plot.margin=margin(15,15,15,15))

ggsave(file.path(PROJ_DIR, "figures", "Figure1_Study_Design.png"), p, width=16, height=12, dpi=300)
message("Figure 1 saved")
