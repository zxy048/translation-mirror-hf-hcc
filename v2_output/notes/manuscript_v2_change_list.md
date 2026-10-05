# Manuscript v1 → v2 change list

Source: `Manuscript_PLOSONE_v1.md` (249 lines, the version submitted to PLOS ONE).
Built by reading v1 in full against the R00–R13 outputs and `Table_S1`–`Table_S31`.

**Nothing in this list has been applied yet.** It is the inventory the rewrite works from.

Reading key: **[N]** = number changes · **[C]** = conclusion changes · **[W]** = wording only
· **[?]** = needs a decision before it can be written.

---

## 1. Title, abstract, author block (lines 1–17)

| line | v1 | action |
|---|---|---|
| 1 | Title: "Disease-context-dependent organization of translation-related transcriptional programs…" | **[C]** Survives, but only in the narrowed sense: the shared object is the module **architecture** (R04 shows the HF translation module IS preserved in HCC, Z = 11.03, rank 2/23); what differs is the **direction of association with disease**. The title cannot imply a difference in module structure. Keep as-is only together with a Results ordering that states preservation before the mirror. |
| 17 | Abstract: "Translation-related modules **differed** between diseases (HCC blue, 1,315 genes…; HF green, 227 genes…), with limited gene-level overlap (62/227; Fisher's OR = 1.7, P = 0.039)" | **[C][N]** All four numbers are unreproducible under the harmonised pipeline. Fisher's test is removed entirely (R2 M4/R3 M10). Replace with: identical rules on a common universe + preservation result + the corrected module sizes from `Table_S4`. |
| 17 | "24 of 33 … mirror perturbation (ρ = −0.598, P = 0.0003)" | **[N]** → Tier A 10/15 mirror sets, sample-level permutation **P = 1e-4**; cluster level 3, P = 0.0031. Tier C 21/33, P = 1e-4. See `R09_mirror.rds`. |
| 17 | "absent in same-organ controls (HCM ρ = −0.036; cirrhosis ρ = +0.402)" | **[N]** v2 has no cirrhosis cohort in R06 (GSE89377 not loaded). Either load it or drop the cirrhosis control — **decision needed**. HCM arm: ρ = −0.100 under v2's set definition. |
| 17 | "ATF4 … (ρ = +0.500); MYC target … (ρ = +0.753)" | **[N]** Recompute from `Table_S30`. Also **[C]** the words "dominant regulatory dimensions" must go (R4 #2). |
| 19–25 | Introduction | **[W]** The framing "does the same program operate in both, with context determining only direction" is now *supported* rather than merely posed — R04 answers the structural half. One paragraph should say so. |
| 3–13 | Author block | unchanged. Corresponding author ORCID is present — satisfies JR #3. |

## 2. Methods

| line | v1 | action |
|---|---|---|
| 35 | TCGA "3 recurrent tumor samples were excluded"; blind = TRUE, nsub = 1000 | **[W]** Keep the numbers (verified: 424 = 371 + 3 + 50; analysed 421). Reword to state the *mechanism* (sample-type filter) rather than a design decision, plus the verified independence reason, plus disclose that 53 patients contribute both tumour and adjacent normal. See response letter R3 minor 9. |
| 37 | "GSE141198 (HCC **validation set**)" | **[W]** → HCC discovery / network-construction cohort (R2 minor 4). |
| 39 | GSE14520 "n = 221 HCC tumors" (GPL3921 only); GSE76427 "n = 115 tumors" | **[N]** v2 loads GSE14520 split by platform: GPL3921 n = 445, GPL571 n = 41; GSE76427 n = 167. State per-platform and per-group counts from `Table_S11`. |
| 43 | "top 6,000 most variable genes retained… **GH** β = 17 vs **HCC** β = 4… parameters were **independently optimized for each dataset**, rather than being artificially matched" | **[N][C]** Entire paragraph replaced. Both networks now use **disease samples only** under identical frozen settings; power selected by one rule (≥ 0.85 scale-free R²); variance pre-filtering removed. Current `Table_S5` carries the fitted values for both. This sentence is the one R1 #1 and R3 #3 were written about. |
| 45 | "manual stepwise procedure… for platform compatibility" vs "blockwiseModules one-step" | **[N][C]** Both diseases now use the same function. Paragraph removed/rewritten. |
| 47 | GO keyword translation-enrichment score; "blue module identified… using the same keyword-based approach" | **[N][C]** Replaced by ssGSEA against a priori canonical sets; the full module × set score matrix is published (`Table_S8`), module designation from `Table_S9`. |
| 49 | Hub genes: \|MM\| > 0.80 & \|GS\| > 0.20; "top ten hub genes … S1 Table" | **[N]** Effect-size convention changed; GS threshold: v1's \|GS\| > 0.20 → **FDR < 0.05** (R3 M4). `Table_S31` shows this is *more* stringent than \|GS\| > 0.50 on this data (GSE57338 purple: 9 genes at \|GS\|>0.50 vs 107 at FDR<0.05). Recompute the hub table. |
| 51 | Module internal robustness: intra-modular TOM Z = 62.8, jackknife CV = 1.6%, cross-disease coherence Z = 2.9, P = 0.008 | **?** **This analysis is not in the R00–R13 pipeline.** v1 reports it; v2 has no script or table for it. Either re-run it on the new modules (TOM is already computed in R02, so this is cheap and reproduces v1's own claim rather than adding new biology) or drop it and say so. **Decision needed.** |
| 55 | "msigdbr… all 50 Hallmark + KEGG ribosome + Reactome **keyword-filtered** (32 sets) → 83 input → 81 pathways" | **[N][C]** Keyword filtering is exactly what R2 M6 / R3 M5 objected to. Replaced by a priori tiers (A/B/C). State the tier definitions and the exact counts from `Table_S2`, `Table_S13`. |
| 57 | Cohen's *d*; "10,000-iteration permutation test… where the vector of **HF pathway effect sizes** was randomly shuffled" | **[N][C]** Three changes: Hedges g with n-weighted pooled SD + bootstrap CI (`R07`); permutation of **disease labels at the sample level** (`R09`); mirror definition now requires opposite signs **AND** \|g\| ≥ 0.20 in both **AND** CI excluding 0 in both. |
| 61 | TATS = mean z-score of 227 green module genes; KM/Cox | **[N][C]** TATS replaced as primary by **TCS** (a priori canonical sets, disease-independent); module-eigengene scores reported in both directions as secondary. **Every cohort's exact contributing gene count** must appear as a table column (`Table_S27`), not as a single number in the text. |
| 65 | "Nineteen candidate TFs" | unchanged in spirit; **[C]** results now described as partial correlations adjusting for purity/proliferation/immune (`Table_S30`), and the words "regulatory"/"drives"/"mediates" removed (R4 #2). |
| 71 | Direction consistency: log2FC vectors; "the **HCC-side log2FC vector** was randomly shuffled" | **[N][C]** Same label-permutation defect as line 57. Replaced by the sample-level null in `R09`. |
| 75 | GSE116250 "**partial replication**"; n = 64 = 14 NF + 37 DCM + 13 ICM | **[N][C]** "Partial replication" removed everywhere (R4 #1). Report: Tier A 2/14, P = 0.174 — non-significant, i.e. not replicated. |
| 79 | "hypertrophic cardiomyopathy (HCM) vs. non-failing myocardium (GSE141910; n = 28 HCM, n = 166 NF)" | **[N][C]** Full composition must be disclosed: 366-sample cardiomyopathy collection = **166 DCM + 28 HCM + 6 PPCM vs 166 NF**. v1's HCM-only use must be justified or replaced by the pooled 200. R2 minor 5 / R3 minor 8. |
| 81–83 | LLM usage: "GPT-4 … for English language editing" | **[W]** Keep. Confirm it still describes what was done. |
| 91 | Statistical analysis: "two-sided P < 0.05", "10,000 iterations, set.seed(42)", "BH FDR", package versions | **[N]** Add: Hedges g + J correction; n-weighted pooled SD; bootstrap CI method; sample-level permutation; Jaccard clustering + n_eff (Nyholt 2004; Li & Ji 2005); I²/Q for heterogeneity; effective-test correction; the frozen-config hash. Update the package version list to the versions actually used by R00–R13 (`R00_config.R` records them). |

## 3. Results

| line | v1 | action |
|---|---|---|
| 95 | Heading "**Disease-specific** network remodeling identifies…" | **[C]** Heading must change. R04 shows the translation module is preserved; nothing here is "disease-specific" in the structural sense. |
| 97 | GSE57338 n = 313, top 6,000 genes, β = 17, "eleven co-expression modules" | **[N][C]** Network now built on the **177 failing samples only**. New module count and sizes from `Table_S4` (HF: 23 named modules). |
| 99 | green module 227 genes; TOM score 695.1; GO ribosome biogenesis adj P = 2.3e-13; eigengene r = −0.521 | **[N]** All superseded: new module, new size, new GO table. Keep the *structure* of the argument (a translation-enriched module exists and is negatively correlated with HF). |
| 101 | hub genes; "none of the seven previously reported canonical hub genes met these criteria" | **[N]** Recompute under FDR < 0.05. The negative statement may or may not survive — must be re-derived, not carried over. |
| 103 | HCC blue module 1,315 genes; "62 of 227… Fisher OR = 1.7, P = 0.039"; "**limited gene-level overlap suggests… disease-context-dependent network configurations**" | **[C][N]** The Fisher test is removed and the inference it supported is **reversed by R04**. Rewrite around: identical rules, common universe, preservation in both directions (HF purple→HCC Z = 11.03 rank 2/23 **strong**; HCC magenta→HF Z = 7.09 **weak-moderate**), asymmetric, and ~half of modules show no evidence either way. |
| 105 | "their gene-module organization … was **disease-specific**, consistent with disease-context-dependent network organization rather than conservation of a single co-expression module" | **[C]** Must be rewritten. The defensible claim: translation-module architecture is **conserved**; the **disease association** reverses. The TOM robustness sentence must move with the §51 decision. |
| 107 | Fig 1 legend: "62 of 227 genes found; Fisher OR = 1.7, P = 0.039" | **[N]** Redraw panel d as the preservation plot (Zsummary/medianRank both directions). Legend rewritten. Also fixes R2 minor 2 (the uploaded Fig 1 was internally titled "Figure 2"). |
| 111 | "81 pathways (48 Hallmark + 1 KEGG ribosome + 32 Reactome…)" | **[N]** New tier counts. |
| 113 | ρ = −0.290 (all 81); ρ = −0.598 (33 pathways, P = 0.0003, CI [−0.756, −0.340]); perm P = 0.0079 | **[N]** Recompute from `R07`/`R09`. The *shape* of the claim (negative cross-disease correlation, stronger on translation sets) should be re-derived from the new numbers before being asserted. |
| 115 | Fig 2 legend: "Twenty-four of 33 … (**binomial test P < 0.0001**)" | **[N]** Binomial removed. Exact binomial for 24/33 is 0.0135, and on v2 counts Tier A is P = 1.0000. Replace with permutation P and state why the binomial null (p = 0.5) is wrong here. |
| 117 | "**all 33 translation-related pathways showed positive effect sizes in HCC**" | **[C]** **False** — `Table_S2` in the v1 submission has six negative values. Sentence removed (R2 M8). |
| 117 | "24 … mirror (binomial test, P < 0.0001; S4 Table)" | **[N][C]** Binomial removed; permutation test substituted. |
| 117 | "Proliferation-associated pathways (E2F, G2M, MYC V1) largest positive in HCC (d = +1.83 to +2.35)"; HF top: Bile Acid +0.79, IFN-α +0.67 | **[N]** Recompute. |
| 119 | "represents **the strongest quantitative signal in this study** (ρ = −0.598, perm P = 0.0079)" | **[N][C]** Re-derive. Under the corrected null the strongest statement is Tier A P = 1e-4 — still strong, but the number and the framing change. |
| 123–125 | GSE116250: all-pathway ρ = −0.283; translation subset ρ = **+0.249**, P = 0.162; "**This discrepancy may reflect (i) etiology and severity… (ii) platform… (iii) statistical power**" | **[C]** The three speculative explanations are replaced by measurements: median I² = 87.3% (HF) / 88.5% (HCC), median p_Q < 0.001; aetiology tested and **excluded** (GSE141910 DCM arm alone 1/13; DCM-vs-DCM across cohorts only 4/15 same direction); platform tested via same-technology pairings (R3 M1) and **not** a clean explanation, since the primary comparison is itself cross-technology and is the strongest row. |
| 125 | "mirror perturbation is a disease-state-associated transition whose magnitude varies with cohort characteristics rather than a binary feature" | **[C]** Now understated rather than overstated: state plainly that the mirror is a property of the GSE57338-anchored comparisons and is **not recovered in either independent HF RNA-seq cohort**. |
| 127 | cardiac control ρ = −0.036 (HCM, 14/33); liver control ρ = +0.402 (cirrhosis, 9/33, P = 0.021) | **[N]** HCM arm re-derived (ρ = −0.100, 2/13 gated). Cirrhosis: **GSE89377 is not in the v2 pipeline** — load it or drop the control. **Decision needed.** |
| 129 | Fig 3 legend | **[N][C]** Redraw; remove P-values that no longer hold; drop the cirrhosis panel if the cohort is dropped. |
| 131–135 | TATS survival: log-rank P = 0.303; Cox HR 1.68 [0.67, 4.21] | **[N][C]** Recompute under TCS. **[W]** R2 minor 8: shorten this section substantially and move the detail to Supporting Information. Keep the null explicit. |
| 137 | Fig 4 legend | **[N]** Recompute; rescale fonts (R2 minor 1). |
| 141–147 | TF–TATS: MYC V2 ρ = +0.753; ATF4 ρ = +0.500; E2F1 +0.352; NFE2L2 −0.280; MYCL −0.388; "11 of 19 significant" | **[N][C]** Recompute (`Table_S30`), now as partial correlations adjusting for purity/proliferation/immune. **[W]** Remove "consistent with its role as a central transcription factor in the integrated stress response" and any causal framing (R3 M6, R4 #2). Note "de-overlap sensitivity analysis" should be recomputed, not carried over. |
| 145 | Fig 5 legend (a/b/c) — contains "ATF4/ISR (stress-responsive) and MYC (proliferative) represent complementary **transcriptional dimensions**", "**partial replication**", "Fisher's exact test OR = 1.7" | **[C]** All three must go. Panel c (conceptual model) must be redrawn around preservation-then-mirror. |
| 149 | "two **regulatory dimensions**" | **[C]** Removed (R4 #2). |

## 4. Discussion and Conclusions

| line | v1 | action |
|---|---|---|
| 153–155 | "…their gene composition, hub gene identity, and regulatory architecture were **disease-specific**, suggesting that distinct regulatory networks may independently converge…" | **[C]** The central paragraph of the Discussion and it is now wrong in its premise. Rewrite around: preserved architecture + reversed disease association. |
| 157 | same-organ controls argument | **[N][C]** Depends on the cirrhosis decision. |
| 159 | "cohort-dependent heterogeneity… severity, etiology composition" | **[N]** Replace speculation with the measured I², the aetiology exclusion, and the platform analysis. |
| 161 | ATF4/MYC paragraph with ISR framing | **[C]** Remove "pathway-level engagement of ISR/UPR components", "distinct regulatory dimensions", "neither subordinate to the other". Keep as association only. |
| 163 | TATS not prognostic | **[W]** Fine in substance; rename the score and shorten. |
| 165 | Limitations — **"Formal module preservation analysis was not performed owing to cross-platform interpretability challenges; gene-level overlap quantification (Fisher OR = 1.7, P = 0.039) provides an empirical assessment"** | **[C]** This sentence is the direct subject of R1 #2 / R2 M4 and must be replaced with the R04 result. Add: composition adjustment attenuates but does not remove the signal (55/103 pairs survive, median attenuation 34.1%); GSE57338 retains only 10/24 immune markers so its "immune-adjusted" model is myeloid-only; modulePreservation subsamples modules >1,000 genes. |
| 167 | Conclusions: "**ATF4/ISR engagement** suggests that the failing heart activates conserved stress-adaptive programs…", "**disease-stage-specific vulnerability**", "**potential therapeutic entry points**" | **[C]** All three phrases are named in R2 minor 9 and are removed by name. Rewrite the closing paragraph around what an observational cross-disease bulk-transcriptome comparison can and cannot support. |
| 171 | Data availability: github.com/zxy048/translation-mirror-hf-hcc | **[W]** Repo must be restructured first (R3 minor 11): v2 pipeline authoritative, superseded scripts moved to `archive/` with a README stating the hard-coded β = 12 problem. The link stays; its content changes. |

## 5. Supporting Information captions (lines 243–249)

| issue | action |
|---|---|
| Captions exist for **S1 Fig, S3 Fig, S4 Fig only** | S2 Fig caption is missing and **S2 Fig was never embedded** (R2 minor 3). Must be restored and captioned. |
| **No table captions at all** | S1–S7 Tables are cited in the text but never captioned. The regenerated SI has 31 tables (`Table_S1`–`Table_S31`); every one needs a caption and the in-text citations need renumbering. |
| S3 Fig legend: "GSE141198 (HCC, β = 4)"; "GSE57338 (HF, β = 17)" | **[N]** Both powers change; regenerate. |
| S4 Fig legend: "green module (227 genes)… Z = 62.8… cross-disease coherence Z = 2.9, P = 0.008" | **[N] or [?]** Tied to the §51 decision. |
| Numbering continuity | Verify by script at build time (R2 minor 2). |

## 6. References

- Add DOI **10.1002/imt2.70157** (iMeta, "Human Biomarker Navigator") if the recommendation is accepted — one sentence in the Discussion (R1 #6). Verified to be a cross-system biomarker-atlas review, only tangentially relevant. **Decision needed.**
- Cite the methods the revision adds: Hedges (1981) / the J correction; Nyholt (2004); Li & Ji (2005); the I²/Q heterogeneity method; Langfelder et al. on `modulePreservation` (the WGCNA preservation paper is a separate reference from refs 15/16).
- Reference numbering must be rechecked after insertion.

---

## Decisions required before the rewrite can start

1. **TOM robustness analysis (v1 line 51)** — re-run on the new modules, or drop it and disclose?
   Recommendation: re-run. TOM is already computed in R02; this reproduces v1's own claim under the new pipeline rather than adding new biology.
2. **GSE89377 cirrhosis control** — load it into the v2 pipeline, or drop the control?
   Recommendation: load it. It is already on local disk and dropping a same-organ control weakens the disease-context argument at exactly the point R2 M3 and R4 #1 press.
3. **iMeta citation** — accept and cite, or decline?
   Recommendation: cite once in the Discussion. Costs nothing; JR #4 does not require it.
4. **GSE141910 primary grouping** — pooled 200 (current R05) or the DCM arm as primary?
   Recommendation: keep the pooled 200 as primary (every subset fails equally, so no conclusion moves), disclose the composition, and report the stratified arms as sensitivity (already done, `Table_S22b/c/d`).

## Ordering requirement that is not negotiable

The Results must state the **preservation** result before the **mirror** result. If the mirror
comes first, the reader infers a structural difference between the two diseases' translation
modules that R04 shows does not exist. This ordering is the manuscript-level expression of
the R04 finding and it is the one thing that, if got wrong, re-creates the defect v1 was
criticised for.
