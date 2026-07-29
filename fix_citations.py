"""Post-process NCVR docx: convert [N]-style citations to superscript numbers."""
import re
from docx import Document
from docx.shared import Pt

DOCX_PATH = r"D:\R_projects\revision_analysis\Manuscript_NCVR_v1.docx"

doc = Document(DOCX_PATH)

# Pattern: [1], [1,2], [3-5], [1,3,5], [1-3,5]
CITE_PATTERN = re.compile(r'\[([\d,\-–\s]+)\]')

for para in doc.paragraphs:
    # Process each paragraph's runs
    full_text = para.text

    # Quick check: does this paragraph contain citation patterns?
    if not CITE_PATTERN.search(full_text):
        continue

    # Collect all runs into a single text representation
    # We need to preserve non-citation formatting (bold, italic)
    matches = list(CITE_PATTERN.finditer(full_text))
    if not matches:
        continue

    # Work backwards through matches to avoid offset issues
    for match in reversed(matches):
        start, end = match.start(), match.end()
        cite_text = match.group(1)  # e.g., "1,2" or "3-5" or "5"

        # Find which run(s) contain this match
        char_pos = 0
        for run_idx, run in enumerate(para.runs):
            run_text = run.text
            run_end = char_pos + len(run_text)

            if char_pos <= start < run_end:
                # This run contains the match start
                local_start = start - char_pos
                local_end = min(end, run_end) - char_pos

                # Split: before_cite + cite + after_cite
                before = run_text[:local_start]
                after = run_text[local_end:]

                # The citation part - create as superscript, keeping brackets: [1], [1,2]
                run.text = before

                # Create superscript run for the citation including brackets
                sup_run = para._element.makeelement(
                    '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}r',
                    {}
                )
                # Copy run properties from original run
                rPr = run._r.find('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}rPr')
                if rPr is not None:
                    sup_run.append(rPr)  # copy existing formatting

                sup_run_elem = sup_run.makeelement(
                    '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}rPr',
                    {}
                )
                vertAlign = sup_run_elem.makeelement(
                    '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}vertAlign',
                    {'{http://schemas.openxmlformats.org/wordprocessingml/2006/main}val': 'superscript'}
                )
                sup_run_elem.append(vertAlign)
                sup_run.insert(0, sup_run_elem)

                sup_text = sup_run.makeelement(
                    '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t',
                    {'{http://schemas.xmlsoap.org/soap/envelope/}space': 'preserve'}
                )
                sup_text.text = '[' + cite_text + ']'
                sup_run.append(sup_text)

                # Insert superscript run after current run
                run._r.addnext(sup_run)

                # Handle any remaining text from this run
                if after:
                    after_run = para._element.makeelement(
                        '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}r',
                        {}
                    )
                    after_text = after_run.makeelement(
                        '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t',
                        {'{http://schemas.xmlsoap.org/soap/envelope/}space': 'preserve'}
                    )
                    after_text.text = after
                    after_run.append(after_text)
                    sup_run.addnext(after_run)

                break

            char_pos = run_end

doc.save(DOCX_PATH)
print(f"Citations converted to superscript in: {DOCX_PATH}")
