# =============================================================================
# Figure renumbering: Delete Fig 1, renumber Fig 2->1, 3->2, promote S5->3, S6->4
# Combine Fig 4A+4B+Fig5 -> new 3-panel Fig 5
# Renumber supplementary: S4->S3, S7->S4
# =============================================================================

import re

MANUSCRIPT = r"D:\R_projects\revision_analysis\Manuscript_EN_v1_draft.md"

with open(MANUSCRIPT, 'r', encoding='utf-8') as f:
    text = f.read()

print(f"Original length: {len(text)} chars")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1: Delete old Fig 1 (Study Design workflow)
# ═══════════════════════════════════════════════════════════════════════════════

# Remove "(Fig. 1)" parenthetical reference in Introduction
text = text.replace('cohorts (Fig. 1).', 'cohorts.')

# Remove image link line
text = text.replace('![Figure 1](figures/Figure1_Study_Design.png)\n\n', '')

# Remove caption (the block starting with **Fig. 1. ...)
text = re.sub(
    r'\*\*Fig\. 1\. Study design overview\.\*\* .*?perturbation direction is disease-specific\.\n\n',
    '',
    text,
    flags=re.DOTALL
)

print("STEP 1: Deleted Fig 1 (Study Design)")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2: Supplementary renumbering (S4->S3, S7->S4)
# Do these FIRST before main figure renumbering to avoid conflicts
# ═══════════════════════════════════════════════════════════════════════════════

# S4 -> S3
text = text.replace('Supplementary Fig. S4', '⟨⟨⟨S4_TEXT⟩⟩⟩')
text = text.replace('Supplementary Figure S4A', '⟨⟨⟨S4A_IMG⟩⟩⟩')
text = text.replace('Supplementary Figure S4B', '⟨⟨⟨S4B_IMG⟩⟩⟩')
text = text.replace('**Fig. S4.**', '⟨⟨⟨S4_CAP⟩⟩⟩')
text = text.replace('- **Fig. S4.**', '⟨⟨⟨S4_LIST⟩⟩⟩')

# S7 -> S4
text = text.replace('Supplementary Fig. S7', '⟨⟨⟨S7_TEXT⟩⟩⟩')
text = text.replace('- **Fig. S7.**', '⟨⟨⟨S7_LIST⟩⟩⟩')

# Resolve temp markers
text = text.replace('⟨⟨⟨S4_TEXT⟩⟩⟩', 'Supplementary Fig. S3')
text = text.replace('⟨⟨⟨S4A_IMG⟩⟩⟩', 'Supplementary Figure S3A')
text = text.replace('⟨⟨⟨S4B_IMG⟩⟩⟩', 'Supplementary Figure S3B')
text = text.replace('⟨⟨⟨S4_CAP⟩⟩⟩', '**Fig. S3.**')
text = text.replace('⟨⟨⟨S4_LIST⟩⟩⟩', '- **Fig. S3.**')
text = text.replace('⟨⟨⟨S7_TEXT⟩⟩⟩', 'Supplementary Fig. S4')
text = text.replace('⟨⟨⟨S7_LIST⟩⟩⟩', '- **Fig. S4.**')

print("STEP 2: Supplementary renumbered (S4->S3, S7->S4)")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: Promote S5 -> Fig 3 (Cohort dependency)
# ═══════════════════════════════════════════════════════════════════════════════

text = text.replace('(Supplementary Fig. S5)', '⟨⟨⟨S5_TEXT⟩⟩⟩')
text = text.replace('![Supplementary Figure S5]', '⟨⟨⟨S5_IMG⟩⟩⟩')
text = text.replace('**Fig. S5.**', '⟨⟨⟨S5_CAP⟩⟩⟩')
text = text.replace('- **Fig. S5.**', '⟨⟨⟨S5_LIST⟩⟩⟩')  # remove from supp list later

text = text.replace('⟨⟨⟨S5_TEXT⟩⟩⟩', '(Fig. 3)')
text = text.replace('⟨⟨⟨S5_IMG⟩⟩⟩', '![Figure 3]')
text = text.replace('⟨⟨⟨S5_CAP⟩⟩⟩', '**Fig. 3.**')
# S5_LIST will be handled in supplementary list cleanup

print("STEP 3: Promoted S5 -> Fig 3")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: Promote S6 -> Fig 4 (TATS KM + hub GS)
# ═══════════════════════════════════════════════════════════════════════════════

text = text.replace('(Fig. S6A)', '⟨⟨⟨S6_TEXT⟩⟩⟩')
text = text.replace('![Supplementary Figure S6]', '⟨⟨⟨S6_IMG⟩⟩⟩')
text = text.replace('**Fig. S6.**', '⟨⟨⟨S6_CAP⟩⟩⟩')
text = text.replace('- **Fig. S6.**', '⟨⟨⟨S6_LIST⟩⟩⟩')

text = text.replace('⟨⟨⟨S6_TEXT⟩⟩⟩', '(Fig. 4A)')
text = text.replace('⟨⟨⟨S6_IMG⟩⟩⟩', '![Figure 4]')
text = text.replace('⟨⟨⟨S6_CAP⟩⟩⟩', '**Fig. 4.**')
# S6_LIST will be handled in supplementary list cleanup

print("STEP 4: Promoted S6 -> Fig 4")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: Renumber Fig 3 -> Fig 2 (ssGSEA mirror scatter)
# ═══════════════════════════════════════════════════════════════════════════════

text = text.replace('(Fig. 3)', '⟨⟨⟨F3_TEXT⟩⟩⟩')
text = text.replace('![Figure 3]', '⟨⟨⟨F3_IMG⟩⟩⟩')
text = text.replace('**Fig. 3.**', '⟨⟨⟨F3_CAP⟩⟩⟩')

text = text.replace('⟨⟨⟨F3_TEXT⟩⟩⟩', '(Fig. 2)')
text = text.replace('⟨⟨⟨F3_IMG⟩⟩⟩', '![Figure 2]')
text = text.replace('⟨⟨⟨F3_CAP⟩⟩⟩', '**Fig. 2.**')

print("STEP 5: Renumbered Fig 3 -> Fig 2")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6: Renumber Fig 2 -> Fig 1 (Translation Module ID)
# ═══════════════════════════════════════════════════════════════════════════════

text = text.replace('(Fig. 2)', '⟨⟨⟨F2_TEXT⟩⟩⟩')
text = text.replace('![Figure 2]', '⟨⟨⟨F2_IMG⟩⟩⟩')
text = text.replace('**Fig. 2.**', '⟨⟨⟨F2_CAP⟩⟩⟩')

text = text.replace('⟨⟨⟨F2_TEXT⟩⟩⟩', '(Fig. 1)')
text = text.replace('⟨⟨⟨F2_IMG⟩⟩⟩', '![Figure 1]')
text = text.replace('⟨⟨⟨F2_CAP⟩⟩⟩', '**Fig. 1.**')

print("STEP 6: Renumbered Fig 2 -> Fig 1")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7: Renumber Fig 4A -> Fig 5A, Fig 4B -> Fig 5B
# ═══════════════════════════════════════════════════════════════════════════════

text = text.replace('(Fig. 4A)', '⟨⟨⟨F4AT⟩⟩⟩')
text = text.replace('(Fig. 4B)', '⟨⟨⟨F4BT⟩⟩⟩')
text = text.replace('![Figure 4A]', '⟨⟨⟨F4AI⟩⟩⟩')
text = text.replace('![Figure 4B]', '⟨⟨⟨F4BI⟩⟩⟩')

text = text.replace('⟨⟨⟨F4AT⟩⟩⟩', '(Fig. 5A)')
text = text.replace('⟨⟨⟨F4BT⟩⟩⟩', '(Fig. 5B)')
text = text.replace('⟨⟨⟨F4AI⟩⟩⟩', '![Figure 5A]')
text = text.replace('⟨⟨⟨F4BI⟩⟩⟩', '![Figure 5B]')

print("STEP 7: Renumbered Fig 4A->5A, Fig 4B->5B")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 8: Renumber standalone Fig 5 -> combined Fig 5 (Discussion reference)
# The (Fig. 5) in Discussion now references the combined figure
# ═══════════════════════════════════════════════════════════════════════════════

# The Discussion reference "(Fig. 5)" — this now refers to the combined Fig 5
# which includes panels A (TF-TATS), B (Pathway-TATS), C (mechanistic model)
# Keep as "(Fig. 5)" since it references the whole figure
# (No change needed for the text reference)

# But the standalone Fig 5 image and caption need to be merged into the combined layout
# Remove standalone Fig 5 image link:
text = text.replace('![Figure 5](figures/Figure5_Mechanistic_Model.png)\n\n', '')

print("STEP 8: Removed standalone Fig 5 image (will be combined)")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 9: Merge old Fig 4 caption + old Fig 5 caption -> new combined Fig 5 caption
# Old Fig 4 caption starts with "**Fig. 4. Transcription factor..."
# Old Fig 5 caption starts with "**Fig. 5. Disease-context-dependent remodeling..."
# ═══════════════════════════════════════════════════════════════════════════════

# First, find and fix the old Fig 4 caption header
text = text.replace(
    '**Fig. 4. Transcription factor and pathway activity associations with TATS in HCC.**',
    '**Fig. 5. Transcription factor associations, pathway activity, and conceptual model of translation-related transcriptional remodeling in HCC.**'
)

# Now merge the rest: the old Fig 4 caption text + old Fig 5 caption text
# Old Fig 4 caption continues: "(A) Spearman correlation of 19 candidate..."
# Old Fig 5 caption was: "**Fig. 5. Disease-context-dependent remodeling...** Schematic summary..."

# Remove the old Fig 5 caption header and merge its content into Fig 4's caption
# Find: "**Fig. 5. Disease-context-dependent remodeling...** Schematic summary..."
old_fig5_caption_pattern = r'\*\*Fig\. 5\. Disease-context-dependent remodeling of translation-related transcriptional programs\.\*\* '
text = re.sub(old_fig5_caption_pattern, '', text)

print("STEP 9: Merged Fig 4 + Fig 5 captions")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 10: Clean up supplementary figures list
# Remove promoted S5 and S6 entries, update numbering
# ═══════════════════════════════════════════════════════════════════════════════

# Remove S5 entry line (now promoted to Fig 3)
text = text.replace('⟨⟨⟨S5_LIST⟩⟩⟩ Cohort dependency and disease-context analyses of cross-disease mirror perturbation. Four-panel scatter plot of Cohen\'s d effect sizes (x-axis: HCC, TCGA-LIHC) for 33 translation pathways (red) vs. 48 other pathways (gray). (A) HF discovery (GSE57338): ρ = −0.598, 24/33 mirror. (B) Cardiac same-organ control (GSE141910, HCM vs. NF): ρ = −0.036 (p = 0.84), no mirror perturbation. (C) Liver same-organ control (GSE89377, cirrhosis vs. normal): ρ = +0.402 (p = 0.021). (D) Independent HF cohort (GSE116250): all-pathway ρ = −0.283 (p = 0.011); translation subset shows cohort-dependent heterogeneity (ρ = +0.249, 14/33 mirror).\n', '')

# Remove S6 entry line (now promoted to Fig 4)
text = text.replace('⟨⟨⟨S6_LIST⟩⟩⟩ TATS exploratory clinical association and hub gene significance. (A) Kaplan–Meier overall survival curves for TATS (median split) in the GSE141198 HCC cohort (n = 148, 94 events). Log-rank p = 0.303; Cox HR = 1.68, 95% CI [0.67, 4.21], p = 0.265. (B) Gene Significance (Pearson correlation with HF status) for the top 10 green module hub genes in GSE57338.\n', '')

print("STEP 10: Cleaned up supplementary list")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 11: Final cleanup — check for any remaining temp markers
# ═══════════════════════════════════════════════════════════════════════════════

remaining = re.findall(r'⟨⟨⟨[^⟩]+⟩⟩⟩', text)
if remaining:
    print(f"WARNING: {len(remaining)} unresolved temp markers: {remaining}")
else:
    print("STEP 11: All temp markers resolved. OK")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 12: Fix any double spaces introduced by removals
# ═══════════════════════════════════════════════════════════════════════════════

text = re.sub(r'  +', ' ', text)
# Fix ".\n\n\n" -> ".\n\n" (excess blank lines)
text = re.sub(r'\n{4,}', '\n\n\n', text)

# ═══════════════════════════════════════════════════════════════════════════════
# Write back
# ═══════════════════════════════════════════════════════════════════════════════

with open(MANUSCRIPT, 'w', encoding='utf-8', newline='\n') as f:
    f.write(text)

print(f"Final length: {len(text)} chars")
print("═══ Renumbering complete ═══")
