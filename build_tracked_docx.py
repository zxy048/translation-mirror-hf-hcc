# -*- coding: utf-8 -*-
"""
Build a Word track-changes docx for the v1 -> v2 revision, without Word.

Why not Word COM: on this machine Application.CompareDocuments and
Document.Compare both hang -- two attempts, 180 s each, each leaving an orphaned
WINWORD process behind, and Application.CompareDocuments additionally rejects a
bare path string for OriginalDocument ("type mismatch"). Constructing the
revision marks in OOXML is reproducible, inspectable, and does not contend with
the Word instance the user is working in.

Algorithm
  1. Read word/document.xml from both docx packages.
  2. Align paragraphs: difflib at block level, then an order-preserving DP inside
     each replaced block, so that only genuinely corresponding paragraphs get a
     character-level diff. Unpaired paragraphs become clean whole-paragraph
     insert/delete rather than a mangled character-level smear, and so does a
     paired paragraph whose diff fragments past FRAGMENT_CHARS_PER_CHANGE --
     a character-level diff is the right tool for an edited sentence and the
     wrong one for a rewritten one.
  3. Rebuild each changed paragraph's runs, wrapping inserted text in <w:ins> and
     deleted text in <w:del>/<w:delText>. Per-run rPr is carried across so that
     superscripts and subscripts survive the rebuild -- a rebuild that dropped
     them would silently print "Zsummary" and "I2" in the marked-up manuscript.
  4. Write the result by copying v2's package and replacing only document.xml.
  5. Verify by simulating accept-all (must reproduce v2) and reject-all (must
     reproduce v1). A revision file that fails either check is worse than none,
     because it looks authoritative.

Run: python build_tracked_docx.py
"""
import copy
import difflib
import shutil
import sys
import zipfile

from lxml import etree

ROOT = r"D:\R_projects\revision_analysis"
V1 = ROOT + r"\Manuscript_PLOSONE_v1.docx"
V2 = ROOT + r"\Manuscript_PLOSONE_v2.docx"
OUT = ROOT + r"\Manuscript_PLOSONE_v2_tracked.docx"

W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
W = "{%s}" % W_NS
XML_NS = "http://www.w3.org/XML/1998/namespace"
XML_SPACE = "{%s}space" % XML_NS

AUTHOR = "Revision v2"
DATE = "2026-10-03T00:00:00Z"

# A paragraph pair is only diffed word-by-word above this similarity; below it,
# calling the two paragraphs "the same paragraph, edited" produces a wall of
# interleaved insert/delete that hides the actual rewrite.
PAIR_THRESHOLD = 0.55

# A character-level diff of a substantially rewritten paragraph comes out as a
# run of short alternating fragments. In a diff tool that is readable; in Word it
# renders as interleaved rubble -- the revised title, on the first attempt, came
# out as "[-D][+Shared architecture, d]isease-context-dependent [+directi]o[+n: a
# ha]r[-ga][+mo]ni[-z]...", which is not a representation of a title change that
# anyone can read.
#
# The test is change *density* across the changed span, not the raw number of
# hunks and not the paragraph length. Paragraph length is the tempting metric
# and it is the wrong one: the worst smear in the file is a 726-character
# paragraph whose 18 fragments all sit inside one rewritten sentence, so
# measuring against 726 hides it behind 400 characters of untouched text.
#
# Measured over the v1 -> v2 diff, characters per fragment inside the changed
# span separate cleanly, with nothing in between: the fragments that render as
# rubble come in at 7.2, 9.9, 10.3, 12.7, 12.9, 14.3, 26.0 and 28.7, while the
# paragraphs that read fine come in at 65.0 and 147.5. The threshold sits in
# that gap. Above it the paragraph was rewritten rather than edited, and the
# honest representation is the whole old text struck through followed by the
# whole new text -- which is also what a reader would conclude unaided.
FRAGMENT_CHARS_PER_CHANGE = 40


def read_document(path):
    with zipfile.ZipFile(path) as z:
        return etree.fromstring(z.read("word/document.xml"))


def register_namespaces(root):
    for prefix, uri in (root.nsmap or {}).items():
        if prefix:
            etree.register_namespace(prefix, uri)


def para_text(p):
    """Concatenated visible text of a paragraph."""
    out = []
    for node in p.iter():
        if node.tag == W + "t" or node.tag == W + "delText":
            out.append(node.text or "")
    return "".join(out)


def body_children(root):
    return list(root.find(W + "body"))


def char_style_map(p):
    """[(char, rPr_element_or_None)] in document order, across nested runs."""
    out = []
    for r in p.iter(W + "r"):
        rPr = r.find(W + "rPr")
        for t in r.iter():
            if t.tag in (W + "t", W + "delText"):
                for ch in (t.text or ""):
                    out.append((ch, rPr))
    return out


def make_run(text, rPr, deleted):
    r = etree.Element(W + "r")
    if rPr is not None:
        r.append(copy.deepcopy(rPr))
    node = etree.SubElement(r, W + ("delText" if deleted else "t"))
    node.set(XML_SPACE, "preserve")
    node.text = text
    return r


class RevIds:
    def __init__(self):
        self.n = 1

    def next(self):
        self.n += 1
        return self.n


def wrap(run, kind, ids):
    el = etree.Element(W + kind)  # w:ins or w:del
    el.set(W + "id", str(ids.next()))
    el.set(W + "author", AUTHOR)
    el.set(W + "date", DATE)
    el.append(run)
    return el


def rebuild_paragraph(p_v2, text_v1, map_v1, ids):
    """Return a new <w:p> carrying v2's content with v1's deletions marked.

    pPr is taken verbatim from the v2 paragraph, so heading levels, list
    numbering and the PLOS body style all survive.
    """
    pPr = p_v2.find(W + "pPr")
    map_v2 = char_style_map(p_v2)
    text_v2 = "".join(c for c, _ in map_v2)

    # (kind, text, rPr) triples, kind in {"", "ins", "del"}. Each opcode's
    # character range is sliced again at every rPr change: taking one rPr per
    # opcode would silently flatten any superscript or subscript that happens to
    # sit inside an otherwise unchanged stretch.
    pieces = []

    def slices(cmap, lo, hi, kind):
        if hi <= lo:
            return
        start, cur = lo, cmap[lo][1]
        for k in range(lo + 1, hi):
            if cmap[k][1] is not cur:
                pieces.append((kind, "".join(c for c, _ in cmap[start:k]), cur))
                start, cur = k, cmap[k][1]
        pieces.append((kind, "".join(c for c, _ in cmap[start:hi]), cur))

    sm = difflib.SequenceMatcher(None, text_v1, text_v2, autojunk=False)
    opcodes = sm.get_opcodes()
    changed = [o for o in opcodes if o[0] != "equal"]
    if len(changed) > 1:
        # Span in both coordinate systems and take the wider, so that a
        # paragraph which is mostly insertion or mostly deletion is measured
        # the same way as one that is mostly replacement.
        span = max(changed[-1][2] - changed[0][1],
                   changed[-1][4] - changed[0][3])
        if span < len(changed) * FRAGMENT_CHARS_PER_CHANGE:
            # Too fragmented to read. Fall back to one delete and one insert;
            # the rest of this function then renders them as struck-through old
            # text followed by underlined new text, which is the reading a
            # rewritten paragraph should get.
            opcodes = [("delete", 0, len(text_v1), 0, 0),
                       ("insert", 0, 0, 0, len(text_v2))]

    for tag, i1, i2, j1, j2 in opcodes:
        if tag == "equal":
            slices(map_v2, j1, j2, "")
        elif tag == "delete":
            slices(map_v1, i1, i2, "del")
        elif tag == "insert":
            slices(map_v2, j1, j2, "ins")
        else:  # replace: deletions first, so the paragraph reads del-then-ins
            slices(map_v1, i1, i2, "del")
            slices(map_v2, j1, j2, "ins")

    new_p = etree.Element(W + "p")
    if pPr is not None:
        new_p.append(copy.deepcopy(pPr))

    # merge adjacent pieces that share a kind and rPr serialisation, so the
    # marked-up file does not explode into one run per character
    merged = []
    for kind, text, rPr in pieces:
        key = (kind, etree.tostring(rPr) if rPr is not None else None)
        if merged and merged[-1][0] == kind and merged[-1][3] == key[1] and text:
            merged[-1][1] += text
        else:
            merged.append([kind, text, rPr, key[1]])

    for kind, text, rPr, _ in merged:
        if not text:
            continue
        if kind == "":
            new_p.append(make_run(text, rPr, deleted=False))
        elif kind == "ins":
            new_p.append(wrap(make_run(text, rPr, deleted=False), "ins", ids))
        else:
            new_p.append(wrap(make_run(text, rPr, deleted=True), "del", ids))
    return new_p


def whole_paragraph(p_src, kind, ids, text=None):
    """A paragraph whose entire content is inserted or deleted."""
    src_pPr = p_src.find(W + "pPr")
    pPr = copy.deepcopy(src_pPr) if src_pPr is not None else etree.Element(W + "pPr")
    new_p = etree.Element(W + "p")
    new_p.append(pPr)

    # mark the paragraph mark itself as inserted/deleted, which is what tells
    # Word the whole paragraph appears or disappears rather than merely having
    # all of its text changed. w:ins/w:del must come first in a paragraph's rPr
    # (CT_ParaRPr sequences them ahead of rStyle/rFonts/b/i), so it is inserted
    # at position 0 rather than appended.
    rPr = pPr.find(W + "rPr")
    if rPr is None:
        rPr = etree.Element(W + "rPr")
        pPr.insert(0, rPr)
    marker = etree.Element(W + kind)
    marker.set(W + "id", str(ids.next()))
    marker.set(W + "author", AUTHOR)
    marker.set(W + "date", DATE)
    rPr.insert(0, marker)

    src_map = char_style_map(p_src)
    src_text = "".join(c for c, _ in src_map) if text is None else text
    if src_text:
        rPr0 = src_map[0][1] if src_map else None
        new_p.append(wrap(make_run(src_text, rPr0, deleted=(kind == "del")), kind, ids))
    return new_p


def align_block(v1_block, v2_block):
    """Order-preserving alignment of two paragraph lists.

    Returns a list of ("same"|"edit"|"del"|"ins", v1_para_or_None,
    v2_para_or_None). "same" means identical text; "edit" means paired and worth
    diffing.
    """
    n, m = len(v1_block), len(v2_block)
    sim = [[0.0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            a, b = v1_block[i], v2_block[j]
            if not a.strip() or not b.strip():
                continue
            sim[i][j] = difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()

    NEG = -1e9
    dp = [[NEG] * (m + 1) for _ in range(n + 1)]
    back = [[None] * (m + 1) for _ in range(n + 1)]
    dp[0][0] = 0.0
    for i in range(n + 1):
        for j in range(m + 1):
            if dp[i][j] == NEG:
                continue
            if i < n and j < m and sim[i][j] >= PAIR_THRESHOLD:
                sc = dp[i][j] + sim[i][j] - 0.5
                if sc > dp[i + 1][j + 1]:
                    dp[i + 1][j + 1] = sc
                    back[i + 1][j + 1] = ("pair", i, j)
            if i < n:
                sc = dp[i][j] - 0.35
                if sc > dp[i + 1][j]:
                    dp[i + 1][j] = sc
                    back[i + 1][j] = ("del", i, j)
            if j < m:
                sc = dp[i][j] - 0.35
                if sc > dp[i][j + 1]:
                    dp[i][j + 1] = sc
                    back[i][j + 1] = ("ins", i, j)

    ops = []
    i, j = n, m
    while (i, j) != (0, 0):
        kind, pi, pj = back[i][j]
        if kind == "pair":
            ops.append(("pair", pi, pj))
            i, j = pi, pj
        elif kind == "del":
            ops.append(("del", pi, None))
            i, j = pi, pj
        else:
            ops.append(("ins", None, pj))
            i, j = pi, pj
    ops.reverse()
    return ops


def main():
    r1, r2 = read_document(V1), read_document(V2)
    register_namespaces(r2)

    c1, c2 = body_children(r1), body_children(r2)
    p1 = [el for el in c1 if el.tag == W + "p"]
    p2 = [el for el in c2 if el.tag == W + "p"]
    t1 = [para_text(p) for p in p1]
    t2 = [para_text(p) for p in p2]
    print(f"v1: {len(p1)} paragraphs, v2: {len(p2)} paragraphs")

    ids = RevIds()
    new_body = etree.Element(W + "body")
    # index into p1/p2 of the paragraphs consumed so far
    k1 = k2 = 0
    n_ins = n_del = n_edit = n_same = 0

    # Non-paragraph children (the trailing sectPr, and any table) are kept
    # verbatim: sectPr carries the page setup format_plos installed, and the
    # paragraph diff has nothing to say about either.
    non_paras = [copy.deepcopy(el) for el in c2 if el.tag != W + "p"]

    sm = difflib.SequenceMatcher(None, [x.strip() for x in t1],
                                 [x.strip() for x in t2], autojunk=False)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            for k in range(i2 - i1):
                new_body.append(copy.deepcopy(p2[j1 + k]))
                n_same += 1
        elif tag == "insert":
            for k in range(j1, j2):
                new_body.append(whole_paragraph(p2[k], "ins", ids))
                n_ins += 1
        elif tag == "delete":
            for k in range(i1, i2):
                new_body.append(whole_paragraph(p1[k], "del", ids))
                n_del += 1
        else:  # replace
            v1_block = [t1[k] for k in range(i1, i2)]
            v2_block = [t2[k] for k in range(j1, j2)]
            for kind, a, b in align_block(v1_block, v2_block):
                if kind == "pair":
                    ra, rb = i1 + a, j1 + b
                    if t1[ra].strip() == t2[rb].strip():
                        new_body.append(copy.deepcopy(p2[rb]))
                        n_same += 1
                    else:
                        new_body.append(
                            rebuild_paragraph(p2[rb], t1[ra], char_style_map(p1[ra]), ids))
                        n_edit += 1
                elif kind == "ins":
                    new_body.append(whole_paragraph(p2[j1 + b], "ins", ids))
                    n_ins += 1
                else:
                    new_body.append(whole_paragraph(p1[i1 + a], "del", ids))
                    n_del += 1

    for el in non_paras:
        new_body.append(el)

    r2.remove(r2.find(W + "body"))
    r2.append(new_body)

    # count the marks we actually produced
    n_ins_marks = len(r2.findall(".//" + W + "ins"))
    n_del_marks = len(r2.findall(".//" + W + "del"))
    print(f"paragraphs: same={n_same} edited={n_edit} inserted={n_ins} deleted={n_del}")
    print(f"revision marks: w:ins={n_ins_marks} w:del={n_del_marks}")
    if n_ins_marks == 0 or n_del_marks == 0:
        sys.exit("ERROR: one side of the revision is empty -- refusing to write")

    # ---- verify by simulating Word's accept-all / reject-all ----
    accept = copy.deepcopy(r2)
    for el in accept.findall(".//" + W + "del"):
        el.getparent().remove(el)
    for el in accept.findall(".//" + W + "ins"):
        parent = el.getparent()
        idx = list(parent).index(el)
        for child in list(el):
            parent.insert(idx, child)
            idx += 1
        parent.remove(el)
    accept_txt = [para_text(p) for p in accept.find(W + "body")
                  if p.tag == W + "p"]

    reject = copy.deepcopy(r2)
    for el in reject.findall(".//" + W + "ins"):
        el.getparent().remove(el)
    for el in reject.findall(".//" + W + "del"):
        parent = el.getparent()
        idx = list(parent).index(el)
        for child in list(el):
            for dt in child.iter(W + "delText"):
                dt.tag = W + "t"
            parent.insert(idx, child)
            idx += 1
        parent.remove(el)
    reject_txt = [para_text(p) for p in reject.find(W + "body")
                  if p.tag == W + "p"]

    def norm(xs):
        return [x.strip() for x in xs if x.strip()]

    a_ok = norm(accept_txt) == norm(t2)
    r_ok = norm(reject_txt) == norm(t1)
    print(f"accept-all reproduces v2: {a_ok}")
    print(f"reject-all reproduces v1: {r_ok}")
    if not a_ok:
        for x, y in zip(norm(accept_txt), norm(t2)):
            if x != y:
                print(f"  first accept mismatch:\n    got {x[:110]!r}\n    want {y[:110]!r}")
                break
    if not r_ok:
        for x, y in zip(norm(reject_txt), norm(t1)):
            if x != y:
                print(f"  first reject mismatch:\n    got {x[:110]!r}\n    want {y[:110]!r}")
                break
    if not (a_ok and r_ok):
        sys.exit("ERROR: round-trip check failed -- not writing a file that "
                 "misrepresents the revision")

    shutil.copy2(V2, OUT)
    new_xml = etree.tostring(r2, xml_declaration=True, encoding="UTF-8",
                             standalone=True)
    # rewrite the package with the new document.xml, preserving every other part
    tmp = OUT + ".tmp"
    with zipfile.ZipFile(V2) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/document.xml":
                data = new_xml
            zout.writestr(item, data)
    shutil.move(tmp, OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
