"""
Generate final submission docx with figure/table placeholders.
Replaces ASCII art blocks with clean [FIGURE X — INSERT HERE] markers.
"""
import re
import subprocess
import os

INPUT_MD   = r"D:\R_projects\revision_analysis\Manuscript_EN_v1_draft.md"
TEMP_MD    = r"D:\R_projects\revision_analysis\Manuscript_temp.md"
OUTPUT_DOCX = r"D:\R_projects\revision_analysis\Manuscript_EN_v1_draft.docx"
PANDOC     = os.path.expandvars(r"%LOCALAPPDATA%\Pandoc\pandoc.exe")

# ── Figure/Table registry ───────────────────────────────────
# Format: (pattern_to_match, placeholder_block)
# pattern_to_match matches the **Figure X.** caption line

FIGURE_PLACEHOLDERS = [
    # (marker_caption_line, number, filename, short_desc)
    ("**Figure 1. Study design overview.**",
     "FIGURE 1", "figures/Figure1_Study_Design.png",
     "Five-layer cross-disease analytical framework for testing the mirror regulation hypothesis."),

    ("**Figure 2. Independent replication and conservation of the translation "
     "co-expression program across heart failure and hepatocellular carcinoma.**",
     "FIGURE 2", "figures/Figure2_CrossDisease_Conservation.png",
     "Independent replication and conservation of the translation co-expression program across HF and HCC. "
     "(A) Gene dendrogram and module assignment from signed WGCNA of GSE57338 (HF, n=313). "
     "The black module (n=99 genes) was identified as the translation-associated module (β=12, signed R²=0.79). "
     "(B) Gene dendrogram and module assignment from signed WGCNA of GSE141198 (HCC, n=148). "
     "The blue module (n=1,315 genes) was identified as the translation-associated module (β=4, signed R²=0.84). "
     "(C) Shared functional enrichment of translation-associated modules across diseases. "
     "GO Biological Process terms enriched in both HF black and HCC blue modules; "
     "translation-related terms highlighted. These analyses demonstrate that translation-associated "
     "co-expression programs identified independently in HF and HCC share conserved functional characteristics."),

    ("**Figure 3. Cross-disease comparison of pathway effect sizes.**",
     "FIGURE 3", "figures/Figure_ssGSEA_cross_disease.png",
     "Cross-disease comparison of pathway effect sizes. Scatter plot of Cohen's d for 81 pathways in HCC (TCGA-LIHC) vs. HF (GSE57338). Translation pathways (n=33) in red; non-translation (n=50) in gray. Spearman ρ = −0.598 (p = 0.0003; permutation p = 0.0079)."),

    ("**Figure 4. TGS prognostic validation and hub gene cross-disease directionality.**",
     "FIGURE 4", "figures/Figure4_TGS_Validation.png",
     "TGS prognostic validation and hub gene cross-disease directionality. (A) Kaplan–Meier survival curve for TGS in GSE141198 (n=148). Log-rank p = 0.479. (B) Cross-disease log2FC of seven hub genes in HCC vs. HF."),

    ("**Figure 5. Upstream regulator analysis of the translation co-expression program.**",
     "FIGURE 5", "figures/Figure_TF_TGS_correlation.png\n│  File (Panel B): figures/Figure_Pathway_TGS_correlation.png",
     "(A) Spearman correlation of 19 candidate TFs with TGS in TCGA-LIHC (n=371). ATF4 is the strongest correlate (ρ=+0.439, FDR<0.0001). (B) Spearman correlation of Hallmark pathway ssGSEA scores with TGS. MYC Targets V2 shows the strongest correlation (ρ=+0.613, p<0.0001)."),

    ("**Figure 6. Cross-disease mirror regulation framework",
     "FIGURE 6", "figures/Figure6_Mechanistic_Model.png",
     "Cross-disease mirror regulation framework for translation in HF versus HCC. Left (HF): energy deficit and mTORC1 suppression → coordinated downregulation (↓). Right (HCC): ISR activation + MYC-driven proliferative signaling → upregulation (↑). Center: structurally conserved translation co-expression module."),

    # Supplementary figures
    ("**Supplementary Figure S1. Cross-disease directional consistency",
     "FIGURE S1", "figures/Figure_S1_Direction_Consistency.png",
     "Cross-disease directional consistency of six directionally concordant hub genes (EEF1A1, FAU, RPL39, RPL3, RPL32, RPL41). Scatter plot of log2FC in HCC vs. HF. Spearman ρ = 0.486 (p = 0.356). RPS28 (red) excluded from TGS due to opposite directionality; seven-gene ρ = 0.126 (p = 0.788)."),

    ("**Supplementary Figure S2. RPL39 single-gene expression",
     "FIGURE S2", "figures/Figure_S2_GSE141198_Survival.png",
     "RPL39 single-gene expression in GSE141198. Kaplan–Meier curve for RPL39 (median dichotomization) in GSE141198 (n=148)."),

    ("**Supplementary Figure S3. Independence of hub genes",
     "FIGURE S3", "figures/Figure_S3_TF_Fisher_enrichment.png",
     "Independence of hub genes from canonical TF target gene sets. Fisher's exact test results for enrichment of seven hub genes in Hallmark TF target gene sets. All p = 1.00."),

    ("**Supplementary Figure S4. Soft threshold selection",
     "FIGURE S4", "figures/Figure_S4A_SoftThreshold_GSE141198.png\n│  File (Panel B): figures/Figure_S4B_SoftThreshold_GSE57338.png",
     "Soft threshold selection for WGCNA. Scale-free topology model fit (R²) and mean connectivity vs. soft threshold power (β). (A) GSE141198: β=4 (R²=0.84). (B) GSE57338: β=12 (R²=0.79)."),
]

def make_figure_placeholder(num, fname, desc):
    """Generate a clean figure placeholder block."""
    border = "─" * 60
    return f"""
<div style="border: 2px solid #333; padding: 12px; margin: 18px 0; background: #f5f5f5;">

**⬇ [{num} — INSERT IMAGE HERE] ⬇**

**File:** `{fname}`

**Caption:** {desc}

</div>
"""

def make_table_placeholder(num, fname, desc):
    """Generate a clean table placeholder block."""
    return f"""
<div style="border: 2px dashed #666; padding: 12px; margin: 18px 0; background: #fafafa;">

**⬇ [{num} — INSERT TABLE HERE] ⬇**

**Source:** `{fname}`

**Content:** {desc}

</div>
"""

# ── Main processing ─────────────────────────────────────────

with open(INPUT_MD, 'r', encoding='utf-8') as f:
    text = f.read()

lines = text.split('\n')

# ── Step 1: Remove ASCII art blocks and insert figure placeholders ──
new_lines = []
i = 0
while i < len(lines):
    line = lines[i]

    # Check if this line starts an ASCII art block
    stripped = line.strip()
    if stripped.startswith('╔'):
        # Skip all lines until we exit the ASCII block
        while i < len(lines):
            s = lines[i].strip()
            if s == '' or s == '```':
                i += 1
                break
            i += 1
        # Also skip trailing blank lines / ```
        while i < len(lines) and lines[i].strip() in ('', '```'):
            i += 1
        continue

    # Skip bare ``` lines
    if stripped == '```':
        i += 1
        continue

    # Check if this line is a figure caption that needs a placeholder
    inserted = False
    for marker, num, fname, desc in FIGURE_PLACEHOLDERS:
        if stripped.startswith(marker):
            # Insert placeholder BEFORE the caption
            placeholder = make_figure_placeholder(num, fname, desc)
            new_lines.append(placeholder)
            inserted = True
            break

    new_lines.append(line)
    i += 1

# ── Step 2: Remove the "Figure Legends" and "Table Legends" sections ──
# (since placeholders now contain this info)
clean_lines = []
skip_until_end = False
for line in new_lines:
    if line.strip().startswith('## Figure Legends'):
        skip_until_end = True
        continue
    if skip_until_end and line.strip().startswith('## Table 2'):
        skip_until_end = False
        # Add a table placeholder before Table 2
        clean_lines.append('')
        clean_lines.append(make_table_placeholder(
            "TABLE 2", "embedded below",
            "Hub Gene Characteristics and Cross-Disease Module Assignment. Seven hub genes with log2FC, FDR, direction concordance, and module assignment in HF (GSE57338) and HCC (TCGA-LIHC/GSE141198). RPS28 excluded from TGS due to opposite directionality."
        ))
        clean_lines.append('')
    if skip_until_end:
        continue
    clean_lines.append(line)

# ── Step 3: Insert Table 1 and Table 3 placeholders ──
# Table 1: After Methods 4.1 (after the datasets paragraph, before 4.2)
# Table 3: After Results 2.4, before Discussion

final_lines = []
for line in clean_lines:
    final_lines.append(line)
    # Table 1: after "Processing and TGS calculation were performed identically across all validation cohorts."
    if 'Processing and TGS calculation were performed identically across all validation cohorts.' in line:
        final_lines.append('')
        final_lines.append(make_table_placeholder(
            "TABLE 1", "Table1_Cohort_Characteristics.txt",
            "Cohort Characteristics. Summary of three primary cohorts (GSE57338, TCGA-LIHC, GSE141198): sample size, platform, disease/control counts, and key clinical annotations."
        ))
        final_lines.append('')
    # Table 3: after the TF-TGS summary sentence
    if 'Complete TF–TGS correlation statistics for all nineteen candidate TFs are summarized in Table 3.' in line:
        final_lines.append('')
        final_lines.append(make_table_placeholder(
            "TABLE 3", "TF_upstream_result.rds (render as formatted table)",
            "TF–TGS Correlation Summary. Spearman ρ, 95% CI, FDR-adjusted p-values for all 19 candidate transcription factors correlated with TGS in TCGA-LIHC tumors (n=371)."
        ))
        final_lines.append('')

# ── Step 4: Append Supplementary Table placeholders after Table 2 ──
supp_table_placeholders = [
    ("TABLE S1", "Table_S1_Module_Comparison.txt",
     "Module Composition and GO Enrichment Comparison. Side-by-side comparison of translation module characteristics in GSE57338 (HF) and GSE141198 (HCC): module size, gene composition (including snoRNA/scaRNA dominance in HF black module), GO enrichment, hub gene co-localization, and cross-disease synthesis. Note: GSE141198 turquoise module (n=1,858) was flagged by initial GO screen but blue (n=1,665) was selected as the primary translation module based on canonical cytoplasmic translation GO terms."),
    ("TABLE S2", "Table_S2_ssGSEA_Effect_Sizes.csv",
     "Complete ssGSEA Pathway Effect Sizes. Cohen's d, group means, and pooled SD for all 81 pathways (48 Hallmark + 1 KEGG ribosome + 32 Reactome translation) in TCGA-LIHC and GSE57338."),
    ("TABLE S3", "Table_S3_Direction_Consistency.csv",
     "Hub Gene Cross-Disease Directionality. Log2FC, FDR, and direction concordance annotation for seven hub genes in HF vs. HCC. Used for Supplementary Figure S1 and TGS gene selection."),
]

for num, fname, desc in supp_table_placeholders:
    final_lines.append('')
    final_lines.append(make_table_placeholder(num, fname, desc))
    final_lines.append('')

# ── Write temp markdown ──
temp_text = '\n'.join(final_lines)

with open(TEMP_MD, 'w', encoding='utf-8') as f:
    f.write(temp_text)

print(f"Temp markdown written: {TEMP_MD}")
print(f"Lines: {len(final_lines)}")

# ── Convert to docx with pandoc ──
cmd = [
    PANDOC,
    TEMP_MD,
    '-o', OUTPUT_DOCX,
    '--from', 'markdown',
    '--to', 'docx',
    '-V', 'mainfont=Times New Roman',
    '-V', 'fontsize=11',
    '--metadata', 'title=Conserved Translational Co-expression Programs Exhibit Mirror Regulation Across Heart Failure and Hepatocellular Carcinoma',
]

result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode == 0:
    print(f"[OK] DOCX saved: {OUTPUT_DOCX}")
else:
    print(f"Pandoc error: {result.stderr}")

# Clean up temp
# os.remove(TEMP_MD)
print("Done.")
