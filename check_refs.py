# -*- coding: utf-8 -*-
import re

with open("Manuscript_HC_v1.md", encoding="utf-8") as f:
    text = f.read()

def extract_section(title):
    m = re.search(rf'^## {title}\n(.*?)(?=\n## |\Z)', text, re.S | re.M)
    return m.group(1).strip() if m else ""

body = extract_section("1. Introduction") + extract_section("2. Methods") + extract_section("3. Results") + extract_section("4. Discussion")

# Extract ALL numbers from ALL bracket groups, including [n,m] and [n-m]
cited = set()
for m in re.finditer(r'\[([0-9,\-\s]+)\]', body):
    inner = m.group(1)
    for part in inner.split(','):
        part = part.strip()
        if '-' in part:
            a, b = part.split('-')
            cited.update(range(int(a), int(b)+1))
        elif part:
            cited.add(int(part))

ref_count = len(re.findall(r'^\[\d+\] ', text, re.M))
uncited = [n for n in range(1, ref_count+1) if n not in cited]
print(f"Refs listed: {ref_count}")
print(f"Cited: {sorted(cited)}")
print(f"UNCITED refs: {uncited}")

# Show the context of uncited refs
for n in uncited:
    m = re.search(rf'^\[{n}\] (.*)$', text, re.M)
    if m:
        print(f"  [{n}] {m.group(1)[:90]}")
