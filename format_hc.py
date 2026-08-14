# =============================================================================
# Hepatology Communications (LWW/AASLD) docx formatting:
#   Body: Times New Roman 12pt, double (2.0) line spacing
#   H1: 14pt bold, H2/H3: 12pt bold
#   Letter page size, 1-inch (2.54 cm) margins all around
#   Line numbers: continuous
#   All text black
# =============================================================================

from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from lxml import etree
import sys

FONT_NAME = "Times New Roman"

BODY_SIZE = Pt(12)
H1_SIZE = Pt(14)
H2_SIZE = Pt(12)
H3_SIZE = Pt(12)
CAPTION_SIZE = Pt(10)


def add_line_numbering(section):
    sect_pr = section._sectPr
    if sect_pr is None:
        return
    lnNumType = etree.SubElement(sect_pr, qn('w:lnNumType'))
    lnNumType.set(qn('w:countBy'), '1')
    lnNumType.set(qn('w:start'), '1')
    lnNumType.set(qn('w:restart'), 'continuous')


def format_hc(path):
    doc = Document(path)

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

    for para in doc.paragraphs:
        para.paragraph_format.line_spacing = 2.0

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

    # Page layout: Letter (8.5 x 11 in), 1-inch margins
    for section in doc.sections:
        section.page_width = Cm(21.59)
        section.page_height = Cm(27.94)
        section.top_margin = Cm(2.54)
        section.bottom_margin = Cm(2.54)
        section.left_margin = Cm(2.54)
        section.right_margin = Cm(2.54)
        add_line_numbering(section)

    doc.save(path)
    print(f"  HC formatted: {path}")
    print(f"    Body={BODY_SIZE.pt:.0f}pt, H1={H1_SIZE.pt:.0f}pt, double spacing, "
          f"1-inch margins, Letter, line numbers")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        format_hc(sys.argv[1])
    else:
        format_hc(r"D:\R_projects\revision_analysis\Manuscript_HC_v1.docx")
