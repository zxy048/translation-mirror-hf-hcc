import fitz, os

PDF = "D:/R_projects/revision_analysis/v2_output/Supporting_Information_v2.pdf"
OUT = "D:/R_projects/revision_analysis/v2_output/_si_preview"
os.makedirs(OUT, exist_ok=True)

doc = fitz.open(PDF)
print("pages:", doc.page_count)

# find the first page of each item we care about
wanted = {
    "S1_Fig": "S1 Fig.",
    "S1_Table": "S1 Table",
    "S19_Table": "S19 Table",
    "S7_Table": "S7 Table",
    "S4_Fig_placeholder": "S4 Fig.",
}
found = {}
for i in range(doc.page_count):
    t = doc[i].get_text()
    for key, needle in wanted.items():
        if key not in found and needle in t:
            found[key] = i

print("located:", found)
for key, idx in found.items():
    for extra in (0, 1):
        p = idx + extra
        if p >= doc.page_count:
            continue
        pix = doc[p].get_pixmap(dpi=105)
        name = f"{OUT}/{key}_p{p}.png"
        pix.save(name)
        print("  wrote", name)
