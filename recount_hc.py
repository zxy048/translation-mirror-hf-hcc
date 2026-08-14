# -*- coding: utf-8 -*-
import re

with open("Manuscript_HC_v1.md", encoding="utf-8") as f:
    text = f.read()

def extract_section(title):
    m = re.search(rf'^## {title}\n(.*?)(?=\n## |\Z)', text, re.S | re.M)
    return m.group(1).strip() if m else ""

def wc(s):
    return len(s.split())

def section_line_count(title):
    m = re.search(rf'^## {title}\n(.*?)(?=\n## |\Z)', text, re.S | re.M)
    if not m:
        return 0, 0
    block = m.group(1)
    start = text[:m.start()].count('\n') + 2
    end = start + block.count('\n')
    return start, end

for t in ["1. Introduction", "2. Methods", "3. Results", "4. Discussion"]:
    s, e = section_line_count(t)
    print(f"{t:20s} words={wc(extract_section(t)):5d}  lines {s}-{e}")

intro = extract_section("1. Introduction")
methods = extract_section("2. Methods")
results = extract_section("3. Results")
disc = extract_section("4. Discussion")
total = wc(intro) + wc(methods) + wc(results) + wc(disc)
print(f"\nBODY total: {total} (limit 5000, need to cut {total-5000})")
