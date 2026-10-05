import sys, re
import pandas as pd, fitz
sys.path.insert(0, "D:/R_projects/revision_analysis/Code/Revision_v2")
from R19_build_SI_pdf import wrap_cell, TAB

df = pd.read_csv(f"{TAB}/Table_S24_composition_shift.csv")
print("columns:", list(df.columns))
doc = fitz.open("D:/R_projects/revision_analysis/v2_output/Supporting_Information_v2.pdf")
flat = re.sub(r"\s+", " ", "\n".join(doc[i].get_text() for i in range(doc.page_count)))
last = df.columns[-1]
print("last column:", last)
for i, r in df.iterrows():
    v = wrap_cell(r[last])
    if v and v not in flat:
        print(f"  MISSING row {i}: {v!r}  | raw={r[last]!r}")
