"""
Resize all manuscript figures to JTM specifications:
  - Full width: 170 mm (2008 px @ 300 dpi)
  - Max height: 225 mm (2657 px @ 300 dpi)
  - Resolution: 300 dpi
  - Format: PNG (lossless)
Uses Lanczos resampling for high-quality downscaling.
"""
from PIL import Image
import os

JTM_FULL_WIDTH_PX = 2008   # 170 mm @ 300 dpi
JTM_MAX_HEIGHT_PX = 2657   # 225 mm @ 300 dpi
OUT_DIR = r"D:\R_projects\revision_analysis\figures"
DPI = (300, 300)

# Map: (output filename, source filename)
FIGURES = [
    ("Fig1_JTM.png", "Figure2_Translation_Module_Identification.png"),
    ("Fig2_JTM.png", "Figure_ssGSEA_cross_disease.png"),
    ("Fig3_JTM.png", "Figure_S5_Validation_SameOrgan_Controls.png"),
    ("Fig4_JTM.png", "Figure4_TATS_Validation.png"),
    ("Fig5_JTM.png", "Figure5_Combined.png"),
]

for out_name, src_name in FIGURES:
    src_path = os.path.join(OUT_DIR, src_name)
    out_path = os.path.join(OUT_DIR, out_name)

    img = Image.open(src_path)
    orig_w, orig_h = img.size

    # Calculate target size (maintain aspect ratio)
    scale = JTM_FULL_WIDTH_PX / orig_w
    new_w = JTM_FULL_WIDTH_PX
    new_h = int(orig_h * scale)

    # Warn if height exceeds max
    if new_h > JTM_MAX_HEIGHT_PX:
        print(f"  WARNING: {out_name} height {new_h}px ({(new_h/300*25.4):.0f}mm) exceeds max {JTM_MAX_HEIGHT_PX}px (225mm)")
        new_h = JTM_MAX_HEIGHT_PX
        new_w = int(orig_w * (JTM_MAX_HEIGHT_PX / orig_h))
        print(f"    Clamped to width {new_w}px")

    # Resize with Lanczos
    img_resized = img.resize((new_w, new_h), Image.LANCZOS)
    img_resized.save(out_path, dpi=DPI)

    size_kb = os.path.getsize(out_path) / 1024
    w_mm = new_w / 300 * 25.4
    h_mm = new_h / 300 * 25.4
    print(f"  {out_name}: {orig_w}x{orig_h} -> {new_w}x{new_h}px ({w_mm:.0f}x{h_mm:.0f}mm) | {size_kb:.0f} KB")

print("\nDone. JTM-sized figures written to:")
for name, _ in FIGURES:
    print(f"  {os.path.join(OUT_DIR, name)}")
