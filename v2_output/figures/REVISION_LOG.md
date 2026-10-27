# Figure revision log

Every change to a rendered main figure, newest first. The purpose is the same as
`archive/README.md`: a render that differs from the last reviewed one must be
traceable to a decision, not to drift. Each entry records the defect, its cause,
the code change, and measurements taken off the PNGs before and after, so the
next audit can tell a deliberate change from a re-render artefact.

Measurements are in pixels on the 300 dpi PNG (2250 × 2625 px = 7.50 × 8.75 in,
the PLOS ceiling); 1 px = 0.24 pt.

---

## 2026-10-27 — Fig 1, panel d: legend moved inside the panel

**Defect (author-reported, confirmed).** Panel d's colour legend — "HCC network
tested in HF" / "HF network tested in HCC" — was laid out as a single row above
the panel: 903 px wide (3.01 in), x 1024–1927, y 1532–1554. Panel d's own plot
area is ~600 px wide and its column starts at x ≈ 1106, so the row began 75 px
(0.25 in, 18 pt) **left of panel b's right edge** (x 1099, the grey HF/HCC
strip), in the gutter between the two columns, at the same height as the HCC
half of panel b's bar chart. It read as part of panel b rather than as panel d's
key, and at 600 dpi the clearance to b was under 20 pt with nothing separating
them.

**Cause.** `th + theme(legend.position = "top")`. The legend was a row of the
gtable above the panel; ggplot sizes it independently of the panel, so a legend
wider than the panel is drawn at full width with the panel centred underneath.
patchwork aligns panels, not legends, so nothing constrained it to its column.

**Fix.** `Code/Revision_v2/R17_figures.R`, panel d only:

- `legend.position = "bottom"`, `legend.box = "vertical"`,
  `legend.title = element_blank()`;
- `guides(colour = guide_legend(ncol = 1))` — `legend.box` does **not** stack the
  entries of a single colour guide (there is only one legend box), so
  one-key-per-line has to be set on the guide itself;
- `legend.key.height = 9 pt`, `legend.margin = 2/0/0/0`,
  `legend.spacing.y = 1 pt`, so the two-line legend does not eat the panel.

**After.** Two rows totalling 427 px (1.42 in), x 1262–1689, at y 2001–2024 and
2059–2082 — centred under panel d's plot. Left edge 163 px (0.54 in, 39 pt) to
the right of panel b's strip.

**Knock-on reflow — measured, not assumed.** Removing the over-wide legend
changed the widths patchwork allocates, so more moved than the legend itself:

| element | before | after |
|---|---|---|
| panel b block | x 316–1099, y 859–1959 | x 316–1062, y 859–1920 |
| panel c block | x 1135–1734, y 840–1281 | x 1135–1734, y 840–1285 |
| panel d block | x 1080–1947, y 1330–1939 | x 1106–1947, y 1334–1934 |

Panel b's plotting area is 37 px (0.12 in) narrower and 39 px (0.13 in) shorter.
Nothing was lost with it: the bar chart, its 43 module labels, the grey
direction strips and the axis title all still fit, and the title block moved by
4 px. No plotted value changed — the figure is read from the same RDS files and
`R17_figures.R` carries no hand-entered numbers.

**Panel-clearance sweep, all three main figures.** Widest fully-white run inside
each inter-panel gutter (i.e. the actual clearance between neighbouring panels):

| figure | gutter | px | in | pt |
|---|---|---|---|---|
| Fig 1 | a \| row (b/c/d) | 49 | 0.163 | 11.8 |
| Fig 1 | b \| c–d (column gutter) | 43 | 0.143 | 10.3 |
| Fig 1 | c \| d | 48 | 0.160 | 11.5 |
| Fig 1 | d legend \| e | 72 | 0.240 | 17.3 |
| Fig 1 | row (b/c/d) \| e | 66 | 0.220 | 15.8 |
| Fig 2 | a \| b | 74 | 0.247 | 17.8 |
| Fig 2 | row 1 \| row 2 | 54 | 0.180 | 13.0 |
| Fig 2 | c \| d | 73 | 0.243 | 17.5 |
| Fig 2 | d \| colour bar | 61 | 0.203 | 14.6 |
| Fig 2 | colour bar \| right edge | 66 | 0.220 | 15.8 |
| Fig 3 | a \| b | 50 | 0.167 | 12.0 |
| Fig 3 | row 1 \| row 2 | 54 | 0.180 | 13.0 |
| Fig 3 | c \| d | 43 | 0.143 | 10.3 |

Minimum across the three figures is 10.3 pt (3.6 mm); nothing required
adjustment.

The two elements specifically flagged for this sweep:

- **Fig 2c legend** ("Observed" / "Permutation null"): x 289–812, y 1583–1625,
  wholly inside panel c's column block (x 131–917). No overflow.
- **Fig 3b legend box** ("concordant" / "discordant"): border x 1918–2129,
  y 929–1070, inside panel b's frame (right edge x 2204, bottom edge y 1172) —
  74 px (17.8 pt) to the right, 102 px (24.5 pt) below. The nearest other
  element is the last green set label, whose bottom is y 897: 32 px (7.7 pt)
  above the box border, and no collision. No data point is hidden under the box
  — colour segmentation finds only the box's own two keys inside its rectangle.

**Verification.**

- `build_figure_tiffs.py`: Fig1–3 at 4500 × 5250 px, 600 dpi, 1.95 / 1.91 /
  1.62 MB.
- `verify_all.py` figure checks: PASS on all three TIFFs (7.50 × 8.75 in,
  600 dpi, under the 10 MB ceiling).
- Text layer of all three figure PDFs scanned for v1 titles ("Translation
  Module", "TATS", "Figure 1/2/3", "Figure S", "WGCNA"): none. Every internal
  title is the v2 one.
- Fig 2 and Fig 3 were re-rendered from unchanged code. Frame-to-frame pixel
  difference peaks at 21/255 grey on antialiased edges with no element moving:
  re-render noise, not a change.

**Not changed.** Manuscript text and the Fig 1 caption (L418–449); the response
letter's Minor 1 / Minor 2 statements; `01`/`02`/`03`/`04` and `S1_File.pdf`.

**Superseded render.** `v2_output/audit/pre_fig1d_legend_20261027/` holds the
PNG and PDF as they were before this change. That directory is not published
(`/v2_output/audit/` is gitignored), so the measurements above are the durable
record; the previous bytes are recoverable from the commit before this one.
