# =============================================================================
# NCVR-specific docx formatting:
#   Body: 11pt, 1.5 line spacing
#   H1: 12pt bold, H2: 11pt bold, H3: 11pt bold
#   Figure captions: 9pt
#   A4, margins: top/bottom 2.5cm, left/right 2cm
# =============================================================================

from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
import sys

FONT_NAME = "Times New Roman"

BODY_SIZE = Pt(11)
H1_SIZE = Pt(12)
H2_SIZE = Pt(11)
H3_SIZE = Pt(11)
CAPTION_SIZE = Pt(9)


def format_ncvr(path):
    doc = Document(path)

    # --- Document-level styles ---
    style_map = {
        "Normal":    (BODY_SIZE, False, 1.5),
        "Heading 1": (H1_SIZE, True, 1.5),
        "Heading 2": (H2_SIZE, True, 1.5),
        "Heading 3": (H3_SIZE, True, 1.5),
    }

    for style_name, (size, bold, line_spacing) in style_map.items():
        style = doc.styles[style_name]
        font = style.font
        font.name = FONT_NAME
        font.size = size
        font.bold = bold
        pf = style.paragraph_format
        pf.space_before = Pt(0)
        pf.space_after = Pt(6)
        pf.line_spacing = line_spacing

    # --- Paragraphs ---
    for para in doc.paragraphs:
        para.paragraph_format.line_spacing = 1.5

        # Detect figure captions (Fig. N | ... or Fig. N. ...)
        is_caption = False
        text = para.text.strip()
        if text.startswith('Fig.') and ('|' in text or text.startswith('Fig. ')):
            is_caption = True

        for run in para.runs:
            run.font.name = FONT_NAME

            if para.style.name.startswith('Heading 1'):
                run.font.size = H1_SIZE
                run.font.bold = True
            elif para.style.name.startswith('Heading 2'):
                run.font.size = H2_SIZE
                run.font.bold = True
            elif para.style.name.startswith('Heading 3'):
                run.font.size = H3_SIZE
                run.font.bold = True
            elif is_caption:
                run.font.size = CAPTION_SIZE
                run.font.bold = True  # figure captions bold per NCVR
            else:
                run.font.size = BODY_SIZE

        if is_caption:
            para.alignment = WD_ALIGN_PARAGRAPH.LEFT

    # --- Page layout: A4, NCVR margins ---
    for section in doc.sections:
        section.page_width = Cm(21.0)
        section.page_height = Cm(29.7)
        section.top_margin = Cm(2.5)
        section.bottom_margin = Cm(2.5)
        section.left_margin = Cm(2.0)
        section.right_margin = Cm(2.0)

    doc.save(path)
    print(f"  NCVR formatted: {path}")
    print(f"    Body={BODY_SIZE.pt:.0f}pt, H1={H1_SIZE.pt:.0f}pt, H2={H2_SIZE.pt:.0f}pt, "
          f"H3={H3_SIZE.pt:.0f}pt, Caption={CAPTION_SIZE.pt:.0f}pt, 1.5 spacing")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        format_ncvr(sys.argv[1])
    else:
        format_ncvr(r"D:\R_projects\revision_analysis\Manuscript_NCVR_v1.docx")
