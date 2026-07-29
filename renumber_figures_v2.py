# =============================================================================
# Figure renumbering v2: Single-pass old->temp, then temp->new
# Fixes double-replacement bug from v1
# =============================================================================

import re

MANUSCRIPT = r"D:\R_projects\revision_analysis\Manuscript_EN_v1_draft.md"

with open(MANUSCRIPT, 'r', encoding='utf-8') as f:
    text = f.read()

print(f"Original length: {len(text)} chars")

# =============================================================================
# PHASE 1: All old patterns -> unique temp markers
# Ordered: longer strings first within each category to avoid partial matches
# =============================================================================

old_to_temp = [
    # --- Delete Fig 1 ---
    ('consistent across independent cohorts (Fig. 1).',
     'consistent across independent cohorts.'),
    ('![Figure 1](figures/Figure1_Study_Design.png)\n\n'
     '**Fig. 1. Study design overview.** The analytical framework comprises five layers: '
     '(A) WGCNA-based identification of translation-related co-expression modules in HF '
     'myocardium (GSE57338, n = 313) with systematic annotation-based module scoring; '
     '(B) independent WGCNA in HCC (GSE141198, n = 148); (C) parallel ssGSEA pathway '
     'activity scoring across both diseases with quantitative cross-disease effect-size '
     'comparison; (D) exploratory evaluation of the translation-associated transcriptional '
     'score (TATS) in three HCC cohorts; (E) systematic screening for upstream '
     'transcriptional regulators. The framework tests whether translation-related programs '
     'exhibit shared network organization or disease-context-dependent remodeling, and '
     'whether perturbation direction is disease-specific.\n\n',
     '\n'),

    # --- S4 -> S3 (supplementary renumbering) ---
    ('Supplementary Fig. S4',   '%%%S4_TEXT%%%'),
    ('Supplementary Figure S4A', '%%%S4A%%%'),
    ('Supplementary Figure S4B', '%%%S4B%%%'),
    ('Fig. S4.',                 '%%%S4_CAP%%%'),

    # --- S7 -> S4 ---
    ('Supplementary Fig. S7',    '%%%S7%%%'),

    # --- S5 -> Fig 3 (PROMOTE to main) ---
    ('(Supplementary Fig. S5)',  '%%%S5_TEXT%%%'),
    ('Supplementary Figure S5',  '%%%S5_IMG%%%'),
    ('Fig. S5.',                 '%%%S5_CAP%%%'),

    # --- S6 -> Fig 4 (PROMOTE to main) ---
    ('(Fig. S6A)',               '%%%S6_TEXT%%%'),
    ('Supplementary Figure S6',  '%%%S6_IMG%%%'),
    ('Fig. S6.',                 '%%%S6_CAP%%%'),

    # --- Fig 4A/B -> Fig 5A/B (before Fig 3/2 to avoid substring issues) ---
    ('(Fig. 4A)',      '%%%F4A_TEXT%%%'),
    ('(Fig. 4B)',      '%%%F4B_TEXT%%%'),
    ('Figure 4A',      '%%%F4A_IMG%%%'),
    ('Figure 4B',      '%%%F4B_IMG%%%'),

    # --- Fig 3 -> Fig 2 ---
    ('(Fig. 3)',       '%%%F3_TEXT%%%'),
    ('Figure 3',       '%%%F3_IMG%%%'),
    ('Fig. 3.',        '%%%F3_CAP%%%'),

    # --- Fig 2 -> Fig 1 ---
    ('(Fig. 2)',       '%%%F2_TEXT%%%'),
    ('Figure 2',       '%%%F2_IMG%%%'),
    ('Fig. 2.',        '%%%F2_CAP%%%'),
]

missing = []
for old, temp in old_to_temp:
    if old in text:
        text = text.replace(old, temp)
    else:
        missing.append(old[:70])

if missing:
    print(f"WARNING: {len(missing)} patterns not found:")
    for m in missing:
        print(f"  - {m}")
else:
    print("PHASE 1: All patterns matched.")

# Verify no old figure numbers remain (check for common patterns)
for check in ['(Fig. 1)', '(Fig. 2)', '(Fig. 3)', '(Fig. 4A)', '(Fig. 4B)',
              'Figure 1', 'Figure 2', 'Figure 3', 'Fig. S4.', 'Fig. S5.', 'Fig. S6.',
              'Supplementary Fig. S4', 'Supplementary Fig. S5', 'Supplementary Fig. S7']:
    if check in text:
        print(f"  RESIDUAL: '{check}' still present!")

# =============================================================================
# PHASE 2: All temp markers -> final values
# =============================================================================

temp_to_new = {
    '%%%S4_TEXT%%%':  'Supplementary Fig. S3',
    '%%%S4A%%%':      'Supplementary Figure S3A',
    '%%%S4B%%%':      'Supplementary Figure S3B',
    '%%%S4_CAP%%%':   'Fig. S3.',
    '%%%S7%%%':       'Supplementary Fig. S4',

    '%%%S5_TEXT%%%':  '(Fig. 3)',
    '%%%S5_IMG%%%':   'Figure 3',
    '%%%S5_CAP%%%':   'Fig. 3.',

    '%%%S6_TEXT%%%':  '(Fig. 4A)',
    '%%%S6_IMG%%%':   'Figure 4',
    '%%%S6_CAP%%%':   'Fig. 4.',

    '%%%F3_TEXT%%%':  '(Fig. 2)',
    '%%%F3_IMG%%%':   'Figure 2',
    '%%%F3_CAP%%%':   'Fig. 2.',

    '%%%F2_TEXT%%%':  '(Fig. 1)',
    '%%%F2_IMG%%%':   'Figure 1',
    '%%%F2_CAP%%%':   'Fig. 1.',

    '%%%F4A_TEXT%%%': '(Fig. 5A)',
    '%%%F4B_TEXT%%%': '(Fig. 5B)',
    '%%%F4A_IMG%%%':  'Figure 5A',
    '%%%F4B_IMG%%%':  'Figure 5B',
}

for temp, new in temp_to_new.items():
    text = text.replace(temp, new)

# Verify no temp markers remain
residual = re.findall(r'%%%[^%]+%%%', text)
if residual:
    print(f"ERROR: {len(residual)} unresolved temp markers: {residual}")
else:
    print("PHASE 2: All temp markers resolved.")

# =============================================================================
# PHASE 3: Merge captions (Fig 4 + Fig 5 -> combined Fig 5)
# =============================================================================

# Replace old Fig 4 caption header
text = text.replace(
    '**Fig. 4. Transcription factor and pathway activity associations with TATS in HCC.**',
    '**Fig. 5. Transcription factor associations, pathway activity, and conceptual model '
    'of translation-related transcriptional remodeling in HCC.**'
)

# Remove standalone Fig 5 image link (now part of combined Fig 5)
text = text.replace(
    '![Figure 5](figures/Figure5_Mechanistic_Model.png)\n\n',
    ''
)

# Remove old Fig 5 caption header (content already merged)
text = text.replace(
    '**Fig. 5. Disease-context-dependent remodeling of translation-related '
    'transcriptional programs.** ',
    ''
)

# Add combined Fig 5 image reference before the merged caption
# The TF and Pathway images are separate; add a combined reference
# Actually, keep the individual 5A/5B images and add the 5C (model) after caption
# Or restructure: add 5C image after the caption
# For now, just clean up

print("PHASE 3: Captions merged.")

# =============================================================================
# PHASE 4: Clean up Supplementary Figures list
# =============================================================================

# Remove promoted S5 entry
text = re.sub(
    r'- \*\*Fig\. S5\.\*\* Cohort dependency and disease-context analyses.*?'
    r'supporting disease-context specificity\.\n',
    '',
    text
)
# Remove promoted S6 entry
text = re.sub(
    r'- \*\*Fig\. S6\.\*\* TATS exploratory clinical association and hub gene significance\.'
    r'[^\n]*\n[^\n]*\n',
    '',
    text
)
# Renumber S4->S3, S7->S4 in supplementary list
text = text.replace('- **Fig. S4.**', '- **Fig. S3.**')
text = text.replace('- **Fig. S7.**', '- **Fig. S4.**')

print("PHASE 4: Supplementary list updated.")

# =============================================================================
# PHASE 5: Final cleanup
# =============================================================================

# Fix double spaces
text = re.sub(r'  +', ' ', text)
# Fix excess blank lines
text = re.sub(r'\n{4,}', '\n\n\n', text)

# =============================================================================
# Write
# =============================================================================

with open(MANUSCRIPT, 'w', encoding='utf-8', newline='\n') as f:
    f.write(text)

print(f"Final length: {len(text)} chars")
print("Done.")
