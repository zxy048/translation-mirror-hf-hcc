# =============================================================================
# JTM (Journal of Translational Medicine) docx formatting:
#   Body: 11pt, double (2.0) line spacing
#   H1: 12pt bold, H2: 11pt bold, H3: 11pt bold
#   Figure captions: 9pt
#   A4, margins: 2.54cm all sides (1 inch)
#   Line numbering: continuous
#   No page breaks
#   All text black
# =============================================================================

from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from lxml import etree
import sys

FONT_NAME = "Times New Roman"

BODY_SIZE = Pt(11)
H1_SIZE = Pt(12)
H2_SIZE = Pt(11)
H3_SIZE = Pt(11)
CAPTION_SIZE = Pt(9)


def add_line_numbering(section):
    """Add continuous line numbering to a section."""
    sect_pr = section._sectPr
    if sect_pr is None:
        return
    # Create line numbering element
    nsmap = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    lnNumType = etree.SubElement(sect_pr, qn('w:lnNumType'))
    lnNumType.set(qn('w:countBy'), '1')
    lnNumType.set(qn('w:start'), '1')
    lnNumType.set(qn('w:restart'), 'continuous')


def format_jtm(path):
    doc = Document(path)

    # --- Document-level styles ---
    style_map = {
        "Normal":    (BODY_SIZE, False, 2.0),
        "Heading 1": (H1_SIZE, True, 2.0),
        "Heading 2": (H2_SIZE, True, 2.0),
        "Heading 3": (H3_SIZE, True, 2.0),
    }

    for style_name, (size, bold, line_spacing) in style_map.items():
        style = doc.styles[style_name]
        font = style.font
        font.name = FONT_NAME
        font.size = size
        font.bold = bold
        font.color.rgb = RGBColor(0, 0, 0)
        pf = style.paragraph_format
        pf.space_before = Pt(0)
        pf.space_after = Pt(6)
        pf.line_spacing = line_spacing

    # --- Paragraphs ---
    for para in doc.paragraphs:
        para.paragraph_format.line_spacing = 2.0

        # Detect figure captions (Fig. N | ...)
        is_caption = False
        text = para.text.strip()
        if text.startswith('Fig.') and '|' in text:
            is_caption = True

        for run in para.runs:
            run.font.name = FONT_NAME
            run.font.color.rgb = RGBColor(0, 0, 0)

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
                run.font.bold = True
            else:
                run.font.size = BODY_SIZE
                run.font.bold = False

        if is_caption:
            para.alignment = WD_ALIGN_PARAGRAPH.LEFT

    # --- Page layout: A4, 2.54cm margins, line numbering ---
    for section in doc.sections:
        section.page_width = Cm(21.0)
        section.page_height = Cm(29.7)
        section.top_margin = Cm(2.54)
        section.bottom_margin = Cm(2.54)
        section.left_margin = Cm(2.54)
        section.right_margin = Cm(2.54)
        add_line_numbering(section)

    doc.save(path)
    print(f"  JTM formatted: {path}")
    print(f"    Body={BODY_SIZE.pt:.0f}pt, H1={H1_SIZE.pt:.0f}pt, H2={H2_SIZE.pt:.0f}pt, "
          f"H3={H3_SIZE.pt:.0f}pt, Caption={CAPTION_SIZE.pt:.0f}pt, double spacing, "
          f"line numbers, A4 2.54cm margins")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        format_jtm(sys.argv[1])
    else:
        format_jtm(r"D:\R_projects\revision_analysis\Manuscript_JTM_v1.docx")
