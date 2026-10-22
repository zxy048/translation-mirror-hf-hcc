<!-- pre-publish-check: policy-documentation -->
<!-- This file describes the exclusion policy. It quotes the phrase "decision
     letter" in prose and therefore trips the content scan; the sentinel above
     exempts it. It contains no peer-review material. -->

# Archive: the superseded analysis round

Everything indexed here produced the **first submission**. It is retained for
provenance and must not be reused or cited. The authoritative pipeline is
`Code/Revision_v2/`; see the repository [README](../Code/README.md) for what the
first analysis got wrong and why it was rebuilt.

## Why the files are not physically in this directory

They bind to their own location by absolute path — all 22 superseded scripts
under `Code/` contain references of the form
`D:/R_projects/revision_analysis/...`, including the `source()` calls that chain
the original run together. Moving them here would break exactly the provenance
this archive exists to preserve, and would leave a set of scripts that look
runnable and are not.

An index that points at code still in place is worth more than a relocation that
quietly breaks it. If the repository is ever reorganised for another reason, the
paths below are what must be rewritten in the same commit.

## Index

### Superseded pipeline — `Code/`

| Path | What it was |
|---|---|
| `Code/WGCNA/` | co-expression network construction, both diseases |
| `Code/ssGSEA/` | pathway scoring |
| `Code/TF_analysis/` | transcription-factor correlation, hub-gene direction, external validation |
| `Code/Figure_generation/` | main and supplementary figures, manuscript docx assembly |

These built the two networks **incomparably** — different gene filters, different
soft-threshold rules, different detection functions. The comparison they support
is not interpretable. See the README section "What the first analysis got wrong".

### Superseded pipeline — repository root

The numbered scripts `00_MASTER_execution_guide.R` through
`19b_fix_mycv2_ssgsea_bootstrap.R` (42 files) are the original run order, chained
by absolute-path `source()` calls from `00_MASTER_execution_guide.R`.

The file most worth reading is **`09_HF_WGCNA_final.R`**, which sets the HF soft
threshold to β = 12 by hand. Its own comment records that the value was chosen
"to ensure Figure 2 legend matches the text", and it computes R² = 0.788 where
the first submission reported 0.865. This is the clearest single example of why
the rebuild was necessary, and it is why the disclosure in the main README exists
rather than a silent correction.

### Submission packages

`PLOS_ONE_Submission/`, `PLOS_ONE_Submission_v2/`, `Frontiers_Submission/`,
`HC_Submission/`, `SciRep_Submission/`, and `BBA_Submission/` are working
directories for manuscript and cover-letter packages, **not** analysis
artifacts. Most are excluded by `.gitignore` and are not published; only
`HC_Submission/` is visible because it was tracked historically.

> **Peer-review correspondence is not published.** The editor's decision letter,
> the reviewers' comments, and the response letters are confidential and are
> excluded via `.gitignore`. `pre_publish_check.py` enforces the same rule by
> reading what git would actually stage, because a path list only covers names
> someone remembered to add.

### The pre-submission check

`qa_plos.py` is the checker that let all of this through. It reads exactly two
files — a manuscript `.md` and its `.docx` — and asserts twelve properties of
them, every one a property of the manuscript text alone. It is kept here as the
worked example of the failure: **a check confined to the document cannot see a
defect that lives in the data.** `verify_all.py` replaces it. The root-level
copy is labelled historical in its header; it targets the v1 files by absolute
path and is not the current QA entry point.
