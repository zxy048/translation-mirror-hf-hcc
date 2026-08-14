# -*- coding: utf-8 -*-
import re

with open("Manuscript_HC_v1.md", encoding="utf-8") as f:
    text = f.read()

# Abstract word count (exclude the "**Background:**" style markers? HC counts body of abstract)
def extract_section(title):
    m = re.search(rf'^## {title}\n(.*?)(?=\n## |\Z)', text, re.S | re.M)
    return m.group(1).strip() if m else ""

abstract = extract_section("Abstract")
# strip markdown emphasis markers for word counting
abstract_clean = re.sub(r'\*\*|\*', '', abstract)
print(f"Abstract words (incl. labels): {len(abstract.split())}")
print(f"Abstract words (labels excluded): {len(abstract_clean.split())}")

# Count references
m = re.search(r'^## References\n(.*?)(?=\n## |\Z)', text, re.S | re.M)
refs = m.group(1).strip()
ref_lines = [l for l in refs.split('\n') if l.strip().startswith('[')]
print(f"References: {len(ref_lines)} (HC limit 50)")

# Verify all in-text citation numbers <= 26 and match a reference
body = extract_section("1. Introduction") + extract_section("2. Methods") + extract_section("3. Results") + extract_section("4. Discussion")
cited = set()
for n in re.findall(r'\[(\d+)[\],]', body):
    cited.add(int(n))
for rng in re.findall(r'\[(\d+)-(\d+)\]', body):
    for n in range(int(rng[0]), int(rng[1])+1):
        cited.add(n)
print(f"In-text citation numbers used: {sorted(cited)}")
print(f"Any citation > 26? {[n for n in cited if n > 26]}")

# Check figure legend count
leg = re.findall(r'\*\*Fig\. \d+\.\*\*', text)
print(f"Figure legends found: {len(leg)} -> {leg}")
