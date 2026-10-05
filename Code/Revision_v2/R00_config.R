# =============================================================================
# R00_config.R -- FROZEN configuration for the PLOS ONE major revision
# =============================================================================
# Every downstream script (R01..R12) sources this file. Nothing here may be
# overridden per-script: the whole point of this revision is that the two
# disease pipelines are provably identical and every gene set is declared
# a priori instead of being scraped by regex.
#
# Provenance: every set name below was verified to EXIST in msigdbr 26.1.0
# (db_version 2026.1.Hs) by Code/Revision_v2/00_probe4.R. Do not add a name
# without re-running that probe.
# =============================================================================

CONFIG_VERSION <- "v2.0-2026-10-01"

# ---------------------------------------------------------------- paths ------
PROJ       <- "D:/R_projects/revision_analysis"
DATA_RAW   <- "D:/R_projects"                     # GSE57338, TCGA_LIHC_se.rds
VALID_DIR  <- file.path(PROJ, "validation_data")
V2         <- file.path(PROJ, "Code", "Revision_v2")
OUT        <- file.path(PROJ, "v2_output")
for (d in c(OUT, file.path(OUT, "tables"), file.path(OUT, "figures"),
            file.path(OUT, "intermediate"), file.path(OUT, "logs"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ------------------------------------------------------------ WGCNA ----------
# Identical for BOTH diseases. assert_wgcna_parity() enforces this at runtime.
WGCNA <- list(
  networkType        = "signed",
  TOMType            = "signed",
  minModuleSize      = 30,
  deepSplit          = 2,
  mergeCutHeight     = 0.25,
  pamRespectsDendro  = FALSE,
  # Memory: a single block of all |U| = 10,213 genes needs ~830 MB per
  # blockSize^2 double (adjacency + TOM + dissTOM + as.dist), i.e. a 4-5 GB
  # peak. This workstation has 15.3 GB total with only ~3.3 GB available, so
  # genes are processed in blocks of 5,000 (peak ~1.5 GB/block), which
  # blockwiseModules then merges by eigengene correlation at mergeCutHeight.
  # IDENTICAL for both diseases -- parity is asserted, not assumed.
  maxBlockSize       = 5000,
  numericLabels      = FALSE,
  verbose            = 3,
  # soft-threshold selection rule (NEVER hard-code a beta):
  powerVector        = c(1:10, seq(12, 20, by = 2)),
  powerR2Cut         = 0.85   # smallest beta reaching scale-free R^2 >= this.
  # There is deliberately NO fallback beta. v1 hard-coded beta=12 and admitted
  # in a public script that it did so "to ensure Figure 2 legend matches the
  # text". If no power reaches powerR2Cut, pick_beta() takes the argmax and
  # flags reached_cut = FALSE so the shortfall is reported, not papered over.
)

# ------------------------------------------------------- statistics ----------
STATS <- list(
  SEED             = 42,
  MIN_EFFECT       = 0.20,     # |Hedges g| gate for "mirror"
  N_BOOT           = 2000,     # per-pathway bootstrap CI
  N_PERM_NULL      = 10000,    # sample-level permutation null
  N_PERM_PRESERV   = 200,      # modulePreservation permutations
  FDR_METHOD       = "BH",
  ALPHA            = 0.05
)

# ----------------------------------------------------------- cohorts ---------
# disease_samples_only = TRUE  -> network built on diseased tissue only, so
# disease status cannot drive module structure (reviewer R1 #1).
COHORTS <- list(
  GSE57338 = list(
    label = "HF (GSE57338)", disease = "HF", tech = "array",
    platform = "GPL11532 (HuGene 1.1 ST)", n_samples = 313,
    file = file.path(PROJ, "GSE57338_series_matrix.txt.gz"),
    role = "discovery", disease_samples_only = TRUE,
    n_disease = 177, n_control = 136),
  GSE141198 = list(
    label = "HCC (GSE141198)", disease = "HCC", tech = "RNA-seq",
    platform = "Illumina HiSeq", n_samples = 148,
    file = file.path(PROJ, "GSE141198_vst.rds"),
    role = "discovery", disease_samples_only = TRUE,
    n_disease = 148, n_control = 0),
  TCGA_LIHC = list(
    label = "HCC (TCGA-LIHC)", disease = "HCC", tech = "RNA-seq",
    platform = "Illumina", n_samples = 424,
    file = file.path(DATA_RAW, "TCGA_LIHC_se.rds"),
    role = "validation_same_tech_for_GSE141198", disease_samples_only = TRUE,
    n_disease = 371, n_control = 50),
  GSE116250 = list(
    label = "HF (GSE116250)", disease = "HF", tech = "RNA-seq",
    platform = "Illumina HiSeq 2500", n_samples = NA,
    file = file.path(VALID_DIR, "GSE116250_rpkm.txt.gz"),
    role = "validation_same_tech_for_GSE57338", disease_samples_only = FALSE),
  GSE14520 = list(
    label = "HCC (GSE14520)", disease = "HCC", tech = "array",
    platform = "GPL3921 + GPL571", n_samples = NA,
    file = c(file.path(DATA_RAW, "GSE14520-GPL3921_series_matrix.txt.gz"),
             file.path(DATA_RAW, "GSE14520-GPL571_series_matrix.txt.gz")),
    annot = c(file.path(DATA_RAW, "GPL3921.soft.gz"),
              file.path(DATA_RAW, "GPL571.soft.gz")),
    role = "validation_cross_tech", disease_samples_only = FALSE),
  GSE76427 = list(
    label = "HCC (GSE76427)", disease = "HCC", tech = "array",
    platform = "GPL10558 (HT-12 V4)", n_samples = NA,
    file = file.path(DATA_RAW, "GSE76427_series_matrix.txt.gz"),
    annot = file.path(DATA_RAW, "GPL10558.soft.gz"),
    role = "validation_cross_tech", disease_samples_only = FALSE),
  GSE141910 = list(
    label = "HF (GSE141910)", disease = "HF", tech = "RNA-seq",
    platform = "Illumina", n_samples = NA,
    file = file.path(VALID_DIR, "GSE141910_series_matrix.txt.gz"),
    raw_dir = file.path(VALID_DIR, "GSE141910"),
    role = "validation_same_tech_for_GSE57338", disease_samples_only = FALSE)
)

# Co-primary discoveries (the two networks we compare)
DISCOVERY <- c("GSE57338", "GSE141198")

# ------------------------------------------ canonical translation sets -------
# THREE TIERS. Tier A is the primary a priori definition. Tier C is the exact
# set used in manuscript v1, kept ONLY for continuity reporting -- it is
# contaminated (viral / post-translational / proliferation / single-gene miRNA)
# and is NOT a defensible "translation" definition.

# --- Tier A: core translation machinery (a priori, PRIMARY) -----------------
SETS_TIER_A <- c(
  "KEGG_RIBOSOME",
  "KEGG_AMINOACYL_TRNA_BIOSYNTHESIS",
  "REACTOME_TRANSLATION",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "REACTOME_TRNA_AMINOACYLATION",
  "REACTOME_CYTOSOLIC_TRNA_AMINOACYLATION",
  "REACTOME_RRNA_PROCESSING",
  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
  "REACTOME_RIBOSOME_ASSOCIATED_QUALITY_CONTROL",
  "REACTOME_RIBOSOME_QUALITY_CONTROL_RQC_COMPLEX_EXTRACTS_AND_DEGRADES_NASCENT_PEPTIDE",
  "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "REACTOME_MITOCHONDRIAL_TRANSLATION_ELONGATION",
  "REACTOME_MITOCHONDRIAL_TRNA_AMINOACYLATION",
  "REACTOME_RRNA_PROCESSING_IN_THE_MITOCHONDRION"
)

# --- Tier B: A + tRNA/rRNA processing & modification (sensitivity) ----------
SETS_TIER_B_EXTRA <- c(
  "REACTOME_TRNA_PROCESSING",
  "REACTOME_TRNA_PROCESSING_IN_THE_NUCLEUS",
  "REACTOME_TRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL",
  "REACTOME_RRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL",
  "REACTOME_TRNA_MODIFICATION_IN_THE_MITOCHONDRION",
  "REACTOME_TRNA_PROCESSING_IN_THE_MITOCHONDRION"
)

# --- Tier C: EXACTLY the 33 sets named in manuscript v1 Table S2 -------------
# Continuity only. 13 of these are genuine translation machinery; 20 are not.
SETS_TIER_C_V1 <- c(
  "REACTOME_CARBOXYTERMINAL_POST_TRANSLATIONAL_MODIFICATIONS_OF_TUBULIN",
  "HALLMARK_MYC_TARGETS_V1",
  "REACTOME_CYTOSOLIC_TRNA_AMINOACYLATION",
  "REACTOME_TRNA_AMINOACYLATION",
  "REACTOME_RIBOSOME_QUALITY_CONTROL_RQC_COMPLEX_EXTRACTS_AND_DEGRADES_NASCENT_PEPTIDE",
  "REACTOME_RIBOSOME_ASSOCIATED_QUALITY_CONTROL",
  "REACTOME_REGULATION_OF_PD_L1_CD274_POST_TRANSLATIONAL_MODIFICATION",
  "REACTOME_DENGUE_VIRUS_GENOME_TRANSLATION_AND_REPLICATION",
  "REACTOME_RESPIRATORY_SYNCYTIAL_VIRUS_RSV_GENOME_REPLICATION_TRANSCRIPTION_AND_TRANSLATION",
  "HALLMARK_MYC_TARGETS_V2",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
  "REACTOME_RRNA_PROCESSING",
  "REACTOME_SARS_COV_2_MODULATES_HOST_TRANSLATION_MACHINERY",
  "REACTOME_REGULATION_OF_CDH1_POSTTRANSLATIONAL_PROCESSING_AND_TRAFFICKING_TO_PLASMA_MEMBRANE",
  "REACTOME_TRANSLATION_OF_SARS_COV_2_STRUCTURAL_PROTEINS",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
  "REACTOME_TRANSLATION",
  "REACTOME_SARS_COV_1_MODULATES_HOST_TRANSLATION_MACHINERY",
  "HALLMARK_MTORC1_SIGNALING",
  "REACTOME_MITOCHONDRIAL_TRNA_AMINOACYLATION",
  "REACTOME_TRANSLATION_OF_REPLICASE_AND_ASSEMBLY_OF_THE_REPLICATION_TRANSCRIPTION_COMPLEX",
  "REACTOME_REGULATION_OF_CDH1_MRNA_TRANSLATION_BY_MICRORNAS",
  "REACTOME_POST_TRANSLATIONAL_MODIFICATION_SYNTHESIS_OF_GPI_ANCHORED_PROTEINS",
  "REACTOME_TRANSCRIPTIONAL_AND_POST_TRANSLATIONAL_REGULATION_OF_MITF_M_EXPRESSION_AND_ACTIVITY",
  "REACTOME_REGULATION_OF_CDH11_MRNA_TRANSLATION_BY_MICRORNAS",
  "REACTOME_TRANSLATION_OF_SARS_COV_1_STRUCTURAL_PROTEINS",
  "REACTOME_REGULATION_OF_PTEN_MRNA_TRANSLATION",
  "REACTOME_MITOCHONDRIAL_TRANSLATION_ELONGATION",
  "REACTOME_RRNA_PROCESSING_IN_THE_MITOCHONDRION",
  "REACTOME_REGULATION_OF_NPAS4_MRNA_TRANSLATION",
  "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "REACTOME_REGULATION_OF_PD_L1_CD274_TRANSLATION"
)

# --- a priori EXCLUSIONS (documented so the choice is auditable) -------------
# Applied when deriving Tier A/B; Tier C is reported unfiltered for continuity.
SETS_EXCLUDED <- list(
  post_translational = c(
    "REACTOME_CARBOXYTERMINAL_POST_TRANSLATIONAL_MODIFICATIONS_OF_TUBULIN",
    "REACTOME_REGULATION_OF_PD_L1_CD274_POST_TRANSLATIONAL_MODIFICATION",
    "REACTOME_REGULATION_OF_CDH1_POSTTRANSLATIONAL_PROCESSING_AND_TRAFFICKING_TO_PLASMA_MEMBRANE",
    "REACTOME_POST_TRANSLATIONAL_MODIFICATION_SYNTHESIS_OF_GPI_ANCHORED_PROTEINS",
    "REACTOME_TRANSCRIPTIONAL_AND_POST_TRANSLATIONAL_REGULATION_OF_MITF_M_EXPRESSION_AND_ACTIVITY",
    "REACTOME_POST_TRANSLATIONAL_PROTEIN_MODIFICATION"),
  viral = c(
    "REACTOME_DENGUE_VIRUS_GENOME_TRANSLATION_AND_REPLICATION",
    "REACTOME_RESPIRATORY_SYNCYTIAL_VIRUS_RSV_GENOME_REPLICATION_TRANSCRIPTION_AND_TRANSLATION",
    "REACTOME_SARS_COV_1_MODULATES_HOST_TRANSLATION_MACHINERY",
    "REACTOME_SARS_COV_2_MODULATES_HOST_TRANSLATION_MACHINERY",
    "REACTOME_TRANSLATION_OF_SARS_COV_1_STRUCTURAL_PROTEINS",
    "REACTOME_TRANSLATION_OF_SARS_COV_2_STRUCTURAL_PROTEINS",
    "REACTOME_TRANSLATION_OF_REPLICASE_AND_ASSEMBLY_OF_THE_REPLICATION_TRANSCRIPTION_COMPLEX"),
  proliferation = c(
    "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_MTORC1_SIGNALING", "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
  single_gene_mirna = c(
    "REACTOME_REGULATION_OF_CDH1_MRNA_TRANSLATION_BY_MICRORNAS",
    "REACTOME_REGULATION_OF_CDH11_MRNA_TRANSLATION_BY_MICRORNAS",
    "REACTOME_REGULATION_OF_PTEN_MRNA_TRANSLATION",
    "REACTOME_REGULATION_OF_NPAS4_MRNA_TRANSLATION",
    "REACTOME_REGULATION_OF_PD_L1_CD274_TRANSLATION")
)

# --- upstream TF panel (Fig 2 / TATS work) ----------------------------------
TFS_PANEL <- c("ATF4", "MYC", "XBP1", "ATF6", "EIF4EBP1", "EIF2AK3", "ERN1",
               "DDIT3", "EIF4E", "EIF4G1", "MTOR", "RPTOR", "TSC1", "TSC2")

# --- signature gene sets for module labelling (R03) -------------------------
SIGNATURE_SETS <- list(
  Translation = SETS_TIER_A,
  Ribosome_biogenesis = c("KEGG_RIBOSOME", "REACTOME_RRNA_PROCESSING"),
  OXPHOS = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
  EMT = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  Cell_cycle = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"),
  Immune = c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  Myogenesis = c("HALLMARK_MYOGENESIS"),
  UPR = c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
  Fatty_acid = c("HALLMARK_FATTY_ACID_METABOLISM"),
  Bile_acid = c("HALLMARK_BILE_ACID_METABOLISM"),
  Xenobiotic = c("HALLMARK_XENOBIOTIC_METABOLISM")
)

# --- marker genes for cell-composition adjustment (R11) ---------------------
# Hand-curated; no extra packages required.
MARKERS <- list(
  immune_pan     = c("PTPRC","CD3D","CD3E","CD3G","CD2","CD8A","CD4","CD19","MS4A1",
                     "NKG7","GNLY","KLRD1","CD68","CD14","FCGR3A","ITGAM","LYZ",
                     "AIF1","TYROBP","CXCL9","CXCL10","IL2RA","FOXP3","CTLA4"),
  myeloid        = c("CD68","CD14","FCGR3A","ITGAM","LYZ","AIF1","TYROBP","CSF1R","ITGAX"),
  lymphocyte     = c("CD3D","CD3E","CD2","CD8A","CD4","CD19","MS4A1","NKG7","GNLY","KLRD1"),
  fibroblast     = c("COL1A1","COL1A2","COL3A1","COL5A1","FN1","DCN","LUM","FAP",
                     "PDGFRB","ACTA2","POSTN","THY1","TAGLN","ELN"),
  endothelial    = c("PECAM1","VWF","CDH5","KDR","ENG","TEK","CLDN5","EGFL7"),
  hepatocyte     = c("ALB","APOA1","APOB","TTR","TF","SERPINA1","CYP3A4","ASGR1",
                     "HP","FGA","FGB","FGG","APOC3","CYP2E1"),
  cardiomyocyte  = c("MYH7","MYH6","TNNT2","TNNC1","TNNI3","ACTC1","MYL2","MYL7",
                     "RYR2","ATP2A2","NPPA","NPPB","TTN","MYBPC3"),
  proliferation  = c("MKI67","TOP2A","PCNA","CCNB1","CDK1","BUB1","AURKA","PLK1",
                     "CCNA2","RRM2","TYMS")
)

# --- ESTIMATE-style rank score gene lists -----------------------------------
ESTIMATE_LISTS <- list(
  stromal = MARKERS$fibroblast,
  immune  = MARKERS$immune_pan
)

# =============================================================================
# loader + assertions
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(msigdbr)
})

#' Build (and cache) the combined MSigDB lookup table. Pulling C2:CP:REACTOME
#' plus KEGG plus Hallmark plus GO:BP takes ~2 min, so cache it once per
#' msigdbr version and reuse across every downstream script.
.msig_cache <- new.env(parent = emptyenv())
msig_table <- function(species = "Homo sapiens") {
  key <- paste0("tbl_", species)
  if (!is.null(.msig_cache[[key]])) return(.msig_cache[[key]])
  cache_file <- file.path(OUT, "intermediate",
                          sprintf("msigdb_%s_%s.rds",
                                  as.character(packageVersion("msigdbr")),
                                  gsub("[^A-Za-z]", "_", species)))
  if (file.exists(cache_file)) {
    tbl <- readRDS(cache_file)
    .msig_cache[[key]] <- tbl
    log_msg("msigdb table loaded from cache: ", basename(cache_file))
    return(tbl)
  }
  pulls <- list(c("H"), c("C2", "CP:REACTOME"), c("C2", "CP:KEGG_LEGACY"),
                c("C2", "CP:KEGG_MEDICUS"), c("C5", "GO:BP"))
  tbl <- do.call(rbind, lapply(pulls, function(p) {
    a <- list(db_species = "HS", species = species, collection = p[1])
    if (length(p) > 1) a$subcollection <- p[2]
    d <- do.call(msigdbr, a)
    data.frame(gs_name = d$gs_name, gene_symbol = d$gene_symbol,
               stringsAsFactors = FALSE)
  }))
  saveRDS(tbl, cache_file)
  .msig_cache[[key]] <- tbl
  log_msg("msigdb table built and cached: ", basename(cache_file),
          " (", nrow(tbl), " rows)")
  tbl
}

#' Fetch MSigDB sets, resolving each name to its subcollection automatically.
#' Errors loudly if any requested name is absent -- a silently missing set
#' would change every downstream count.
resolve_sets <- function(set_names, species = "Homo sapiens") {
  all <- msig_table(species)
  avail <- unique(all$gs_name)
  missing <- setdiff(set_names, avail)
  if (length(missing)) {
    stop("These gene-set names do NOT exist in msigdbr ",
         as.character(packageVersion("msigdbr")), ":\n  ",
         paste(missing, collapse = "\n  "),
         "\nRe-run Code/Revision_v2/00_probe4.R and update R00_config.R.")
  }
  sub <- all[all$gs_name %in% set_names, ]
  lst <- split(sub$gene_symbol, sub$gs_name)
  lst[set_names]
}

#' Hard parity check: the two discovery networks must use byte-identical WGCNA
#' parameters. Called by R02 before either network is built.
assert_wgcna_parity <- function(par_a, par_b) {
  keys <- c("networkType","TOMType","minModuleSize","deepSplit","mergeCutHeight",
            "pamRespectsDendro","maxBlockSize","numericLabels")
  bad <- keys[!vapply(keys, function(k)
    identical(par_a[[k]], par_b[[k]]), logical(1))]
  if (length(bad))
    stop("WGCNA parameter mismatch between the two diseases: ",
         paste(bad, collapse = ", "))
  invisible(TRUE)
}

#' Pick beta by rule, never by hand.
#'
#' WGCNA convention: fitIndices$SFT.R.sq is an unsigned R^2 and the sign of the
#' slope restores the direction, so the usable criterion is
#' -sign(slope) * SFT.R.sq.
#'
#' Two rule-based outcomes, both explicit -- there is NO hard-coded fallback:
#'   "cut"     : smallest power whose fit reaches `cut`
#'   "max_r2"  : the data never reaches `cut`, so take the power that maximises
#'               the fit (still rule-derived; reported as a limitation)
pick_beta <- function(sft, cut = WGCNA$powerR2Cut) {
  fit    <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  powers <- sft$fitIndices[, 1]
  ok <- which(fit >= cut)
  if (length(ok)) {
    i <- ok[1]
    return(list(power = powers[i], r2 = fit[i], method = "cut",
                fit = fit, powers = powers, reached_cut = TRUE))
  }
  i <- which.max(fit)
  list(power = powers[i], r2 = fit[i], method = "max_r2",
       fit = fit, powers = powers, reached_cut = FALSE)
}

#' For a coordinated cross-disease comparison the two networks should be built
#' under the SAME stringency. `choose_common_beta` takes the per-disease picks
#' and returns the most stringent rule-implied power, so neither network is
#' built more permissively than its own data requires.
choose_common_beta <- function(picks) {
  p <- max(vapply(picks, function(x) x$power, numeric(1)))
  list(power = p,
       rationale = paste(vapply(names(picks), function(nm)
         sprintf("%s chose beta=%s (%s, R2=%.3f)", nm, picks[[nm]]$power,
                 picks[[nm]]$method, picks[[nm]]$r2), character(1)),
         collapse = "; "))
}

#' Standard signed Cohen's d / Hedges g with n-weighted pooled SD.
#' Replaces the equal-weight sqrt((s1^2+s2^2)/2) used in v1.
cohens_d <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(c(d = NA, g = NA, sd_pool = NA))
  sd_pool <- sqrt(((n1 - 1) * stats::var(x) + (n2 - 1) * stats::var(y)) /
                  (n1 + n2 - 2))
  if (!is.finite(sd_pool) || sd_pool == 0) return(c(d = NA, g = NA, sd_pool = sd_pool))
  d  <- (mean(x) - mean(y)) / sd_pool
  df <- n1 + n2 - 2
  J  <- 1 - 3 / (4 * df - 1)          # small-sample correction
  c(d = d, g = J * d, sd_pool = sd_pool)
}

#' Li & Ji (2005) effective number of independent tests from a correlation
#' matrix, via the eigenvalues of its correlation matrix.
n_effective <- function(cor_mat) {
  m <- nrow(cor_mat)
  if (m < 2) return(m)
  ev <- eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 0)
  if (sum(ev) == 0) return(m)
  sum(ev)^2 / sum(ev^2)
}

#' Jaccard distance between gene sets, for redundancy clustering.
jaccard_dist <- function(sets) {
  n <- length(sets)
  d <- matrix(0, n, n, dimnames = list(names(sets), names(sets)))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i >= j) next
    a <- sets[[i]]; b <- sets[[j]]
    u <- length(union(a, b)); k <- length(intersect(a, b))
    d[i, j] <- d[j, i] <- if (u == 0) 1 else 1 - k / u
  }
  stats::as.dist(d)
}

#' Provenance stamp appended to every written table.
stamp <- function(df) {
  attr(df, "generated_by") <- "Code/Revision_v2"
  attr(df, "config_version") <- CONFIG_VERSION
  attr(df, "generated_on") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  df
}

write_table <- function(df, name) {
  p <- file.path(OUT, "tables", name)
  utils::write.csv(df, p, row.names = FALSE)
  cat("[written]", p, "\n")
  invisible(p)
}

log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(OUT, "logs", "run.log"), append = TRUE)
}

set.seed(STATS$SEED)
options(stringsAsFactors = FALSE)

# Shared GEO / expression I/O helpers (parse_series_matrix, parse_soft,
# collapse_by_mean, soft_probe2symbol, ensg2symbol, filter_*). Kept in a
# separate file so R01 and R05 cannot diverge.
source(file.path(V2, "R00_io.R"))

log_msg("R00 config loaded: ", CONFIG_VERSION,
        " | Tier A=", length(SETS_TIER_A),
        " Tier B=", length(c(SETS_TIER_A, SETS_TIER_B_EXTRA)),
        " Tier C(v1)=", length(SETS_TIER_C_V1))
