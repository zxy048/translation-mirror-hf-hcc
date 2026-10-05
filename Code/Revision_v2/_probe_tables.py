import glob, os
import pandas as pd

files = sorted(glob.glob("D:/R_projects/revision_analysis/v2_output/tables/Table_S*.csv"))
print(f"{'file':52s} {'rows':>6s} {'cols':>5s} {'maxcell':>8s}  longest header")
worst_cell = (0, "")
worst_cols = (0, "")
for f in files:
    df = pd.read_csv(f)
    mx = 0
    for c in df.columns:
        for v in df[c].head(400):
            v = "" if v is None else str(v)
            if len(v) > mx:
                mx = len(v)
    lh = max((len(str(c)) for c in df.columns), default=0)
    print(f"{os.path.basename(f):52s} {len(df):6d} {len(df.columns):5d} {mx:8d}  {lh}")
    if len(df.columns) > worst_cols[0]:
        worst_cols = (len(df.columns), os.path.basename(f))
    if mx > worst_cell[0]:
        worst_cell = (mx, os.path.basename(f))

print()
print("widest table :", worst_cols)
print("longest cell :", worst_cell)
