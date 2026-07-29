"""
Preprocess manuscript markdown: remove ASCII figure placeholders, clean up formatting.
"""
import re

INPUT_MD  = r"D:\R_projects\revision_analysis\Manuscript_EN_v1_draft.md"
OUTPUT_MD = r"D:\R_projects\revision_analysis\Manuscript_EN_clean.md"

with open(INPUT_MD, 'r', encoding='utf-8') as f:
    text = f.read()

# Remove ASCII figure placeholder blocks (the ╔═══ ... ╚═══ blocks)
# Pattern: starts with optional whitespace + ╔, ends with ╝ or after ```
# These blocks contain lines starting with ╔ ║ ╚
lines = text.split('\n')
cleaned = []
in_block = False

for line in lines:
    stripped = line.strip()

    # Detect start of ASCII art block
    if stripped.startswith('╔'):
        in_block = True
        continue

    if in_block:
        # End of block — blank line or ```
        if stripped == '' or stripped == '```':
            in_block = False
        continue

    # Skip bare ``` lines (code fence markers)
    if stripped == '```':
        continue

    # Remove ``` that are part of figure placeholders (check for ║ lines nearby — already handled above)
    cleaned.append(line)

# Write cleaned markdown
with open(OUTPUT_MD, 'w', encoding='utf-8') as f:
    f.write('\n'.join(cleaned))

print(f"Cleaned markdown written to: {OUTPUT_MD}")
print(f"Original lines: {len(lines)}, Cleaned lines: {len(cleaned)}")
