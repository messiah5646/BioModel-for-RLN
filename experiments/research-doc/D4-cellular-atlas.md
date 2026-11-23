# D4: Cellular atlas

Q20 / Round 110. Mirror of the tab "D4: Cellular atlas" in the claude.ai research doc (artifact Mv7f8uVm8ELsioUZqncEBH), written 2026-11-23. Literature synthesis: no pre-registered test, nothing scored.

First pass, 2026-11-23. Built from about 20 primary, trial or review papers reached through PubMed (abstracts), plus the project's own single-cell tests (rounds 89–104). Evidence labels on every claim; mouse findings marked.

**The answer in one line:** in human SLE and LN, the cells that change *with activity* are the B-lineage effectors (activated naive → DN2/ABC → plasmablast) and their T helpers (Tph), while the interferon-stimulated state is spread across almost every cell type and behaves as a patient trait. Inside the kidney the effector shifts from monocytes/macrophages (early, active) to lymphoid cells (late, scarred). Long-lived plasma cells sit outside all standard therapy and are the clearest cellular reason for refractoriness.

## 1. Atlas table

| Cell state | Where | What changes in SLE / LN | Evidence | Tracks | Drugs that hit it |
|---|---|---|---|---|---|
| Activated naive (aNAV) → DN2 / ABC B cells (CD11c⁺ T-bet⁺ CXCR5⁻) | Blood, kidney | Expanded; TLR7-hyper-responsive, IL-21-driven pre-plasma cells; dominate in African-American patients with active disease and nephritis | Supported | Activity, LN, anti-Sm/RNP | CD19/CD20 depletion, CAR-T; TLR7/8 inhibitors (Hypothetical) |
| Plasmablasts (CD27hi CD38hi) | Blood | Expanded at activity; partly mucosal (IgA, CCR10, β7 in ~40–48%) | Established (expansion), Supported (mucosal share) | Activity (best blood activity signal) | CD19 (not CD20), CD38, BCMA |
| Long-lived plasma cells (CD19⁻ BCMA⁺) | Bone marrow, possibly inflamed tissue | Secrete stable autoantibodies; unresponsive to standard immunosuppression and to CD20 depletion | Supported | Refractory serology | CD38 (daratumumab), BCMA; CD19 CAR-T only partly |
| Tph (PD-1hi CXCR5⁻ CD4) and Tfh | Blood, kidney | Tph highly expanded, correlate with activity and CD11c⁺ B cells; help B cells via IL-21 and MAF; SAP⁺ Tph higher in LN | Supported | Activity, LN | Indirect (CAR-T, IL-21 blockade untested in LN) |
| Tregs | Blood | Relative IL-2 deficiency, Tregs functionally short; low-dose IL-2 expands them | Supported | Activity (inverse) | Low-dose IL-2 |
| Naive CD4 T cells | Blood | Reduced; the fall correlates with monocyte ISG level | Supported | IFN level | — |
| GZMH⁺ CD8 (blood); GZMK⁺ CD8 (kidney) | Blood; kidney interstitium | Blood: expanded, repertoire-restricted cytotoxic clones. Kidney: dominant T cell, sits with atypical B cells, much of it bystander, cytokine-driven (IL-15, IFN-β) | Supported (blood), Emerging (kidney) | Chronicity (kidney) | None specific |
| Classical monocytes (ISG-high) | Blood | Main carrier of the blood IFN signature | Established | Patient IFN set point, not flare | Anifrolumab (signal), CAR-T (falls) |
| Infiltrating monocytes → disease macrophages (TREM2, SPP1, APOE, CD9) | Kidney | Dominate early, active proliferative LN; contract as chronicity rises | Supported | Activity, histological class | None specific |
| pDCs | Blood (fewer), tissue (skin, kidney) | Main type I IFN producers after nucleic-acid immune complexes; sparse in blood | Supported | IFN | Anti-BDCA2 (litifilimab: cutaneous lupus) |
| Low-density granulocytes (LDGs) | Blood (in PBMC fraction) | Pro-inflammatory, make type I IFN, kill endothelium, spontaneous NETosis with oxidised mtDNA | Supported | Vascular damage, IFN | None specific |
| Kidney resident cells (podocytes, tubular, parietal epithelium) | Kidney | IFN response in tubules; podocyte injury via C5a/C5aR1; parietal cells signalled by disease macrophages → sclerosis (see D2) | Supported/Emerging | Response, chronicity | SGLT2i, complement (see D2) |

## 2. B lineage: the activity compartment

- **Supported.** SLE B-cell expansions are mostly double-negative (IgD⁻ CD27⁻) cells. The DN2 subset (CXCR5⁻ CD11c⁺, T-bet network) is a pre-plasma-cell state. It lacks the TLR brake TRAF5, is hyper-responsive to TLR7, and becomes antibody-secreting cells in an IL-21-dependent way. Clonal sharing links activated naive cells, DN2 cells and plasma cells in vivo. DN2 cells predominated in African-American patients with active disease, nephritis and anti-Sm/anti-RNA antibodies ([Jenks 2018, *Immunity*](https://doi.org/10.1016/j.immuni.2018.08.015)).
- **Supported (review).** Human SLE flares feature naive-derived, extrafollicular effector B cells rather than only germinal-centre output ([Immunol Rev 2019](https://doi.org/10.1111/imr.12741)). This fits D1's loop: TLR7 sits at the hinge between innate sensing and B-cell output.
- **Supported.** Blood plasmablasts in SLE are partly mucosal: IgA, CCR10 and β7 integrin on 48%, 40% and 38% of them ([Mei 2017, *Arthritis Rheumatol*](https://doi.org/10.1002/art.40181)). So not every plasmablast is a kidney-relevant autoantibody factory, which blunts plasmablast counts as an LN-specific marker.
- **Supported.** Long-lived plasma cells resist standard immunosuppression. Daratumumab (anti-CD38) depleted them in two life-threatening refractory cases, with falls in type I IFN activity ([Ostendorf 2020, *NEJM*](https://doi.org/10.1056/NEJMoa2023325)). In refractory SLE, blood CD19⁺ B cells and marrow CD19⁻ BCMA⁺ plasma cells are both dominant autoantibody sources; co-infused CD19 + BCMA CAR-T gave LLDAS and DORIS remission in 12 of 15 patients by week 12 ([*Nat Med* 2025](https://doi.org/10.1038/s41591-025-03937-8), phase 1, uncontrolled). A BCMA-CD19 compound CAR cleared autoantibodies "including those derived from long-lived plasma cells" within 3 months in 11 patients ([*Ann Rheum Dis* 2024](https://doi.org/10.1136/ard-2024-225785), uncontrolled).
  - Reading for refractory LN: a refractory patient whose pathogenic antibody comes from CD19⁻ plasma cells is predicted to fail CD20 and partly fail CD19 depletion. This is a testable split (D7 route 2, incomplete target engagement). **Hypothetical** until a trial reports response by baseline plasma-cell source.

## 3. T cells: helpers, a missing brake, and bystanders

- **Supported.** PD-1hi CXCR5⁻ T peripheral helper (Tph) cells are highly expanded in SLE blood (52 patients, several cohorts). Their frequency, not that of Tfh, correlated with disease activity and with CD11c⁺ B cells; they drive B-cell differentiation through IL-21 and MAF, and PD-1hi CD4 cells in LN kidneys correlated with kidney B-cell numbers ([Bocharnikov 2019, *JCI Insight*](https://doi.org/10.1172/jci.insight.130062)). SAP⁺ Tph were higher in LN than non-renal SLE after adjustment for complement and anti-dsDNA (n 30, [*Front Immunol* 2024](https://doi.org/10.3389/fimmu.2024.1327437); Emerging, small).
- **Supported.** IL-2 deficiency limits Tregs. In a 152-patient phase IIb trial, low-dose IL-2 gave dose-dependent SRI-4 response at week 12 (69.7% at 1 M IU vs 23.5% placebo) and expanded Tregs ([*Nat Commun* 2026](https://doi.org/10.1038/s41467-026-72245-7)). In an earlier 60-patient RCT the week-12 primary endpoint was missed (55% vs 30%, p 0.052); LN complete remission was 7/13 vs 2/12 ([*Ann Rheum Dis* 2020](https://doi.org/10.1136/annrheumdis-2019-215396); LN subgroup tiny).
- **Supported.** In 1.2 M PBMCs (162 cases, 99 controls), SLE showed fewer naive CD4 T cells, in proportion to monocyte ISG level, and expansion of repertoire-restricted GZMH⁺ CD8 cytotoxic cells ([Perez 2022, *Science*](https://doi.org/10.1126/science.abf1970)).
- **Emerging.** In the kidney, GZMK⁺ CD8 T cells dominate and many expanded clones are not kidney-specific (D2, section 3d). The T-cell infiltrate is partly bystander.

## 4. Myeloid and granulocytes: the innate amplifiers

- **Established.** The blood IFN signature lives mostly in monocytes but is shared by small ISG-high sub-populations of almost every cell type, including T, NK, DC, B and especially plasma cells; expansion of ISG-high and monogenic-lupus-gene-rich subpopulations marked the most active children ([Nehar-Belaid 2020, *Nat Immunol*](https://doi.org/10.1038/s41590-020-0743-0); 33 children, 11 controls, adult replication).
- **Supported.** In LN kidneys, 21 leukocyte subsets, stepwise monocyte-to-macrophage differentiation in situ, an IFN response in most cells, and close kidney–urine agreement ([Arazi 2019, *Nat Immunol*](https://doi.org/10.1038/s41590-019-0398-x)). Disease macrophages expand with activity and contract as chronicity rises (D2, 3c).
- **Supported (review + primary).** pDCs are the specialist type I IFN producers, triggered through endosomal TLR7/9 by nucleic-acid immune complexes ([*Annu Rev Immunol* 2024](https://doi.org/10.1146/annurev-immunol-090122-041105)). Blocking BDCA2 on pDCs reduced cutaneous lupus activity in phase 2 ([LILAC, *NEJM* 2022](https://doi.org/10.1056/NEJMoa2118024)). Whether pDCs matter in the *kidney* is **unknown**: the open kidney immune atlases contain few or no pDCs.
- **Supported.** Low-density granulocytes in SLE PBMC fractions are activated, secrete type I IFN, TNF and IFN-γ, kill endothelial cells and block endothelial repair ([Denny 2010, *J Immunol*](https://doi.org/10.4049/jimmunol.0902199)). Their spontaneous NETs carry oxidised mitochondrial DNA that drives IFN through STING ([Lood 2016, *Nat Med*](https://doi.org/10.1038/nm.4027); mouse in vivo part). Review: LDGs are a distinct pro-inflammatory neutrophil subset tied to vascular damage ([*Immunol Rev* 2023](https://doi.org/10.1111/imr.13161)).
  - Caution: whether NET clearance failure drives *nephritis* is contested (D2, 3b).

## 5. What the project's own data add

- **Supported (own tests).** Blood IFN in monocytes behaves as a patient set point: it did not rise before flares (round 89, n 6), falls by a shared multiplier after B-cell depletion (round 96), and is *lower* in biopsy-proven LN than in active non-renal SLE (round 99). This matches the atlas: IFN is broadly distributed and slow; B-lineage output is the fast compartment.
- **Emerging (own tests).** The one open LN kidney immune sample re-mapped in round 103 was monocyte-dominated and IFN-high, and the public kidney reference has no neutrophils, so neutrophil/LDG roles in the human kidney cannot be scored from open data yet (Q46).
- **Emerging (own tests).** Cell-mix features from a frozen Perez model carried over only weakly to an independent flare set, and transferred CD8 and plasmablast fractions looked implausible (round 104, Q48). A ribosomal-fraction QC metric carried disease signal (round 102, Q44). Composition estimates need label-free re-annotation before use.

## 6. Mapping to LUMEN

| LUMEN state | Cells it stands for | Speed | Moved by |
|---|---|---|---|
| IFN set point (patient parameter) | ISG-high monocytes, pDCs, LDGs | Slow, trait-like | Anifrolumab; CAR-T (indirectly, via fewer immune complexes) |
| B-effector / activity state | aNAV → DN2 → plasmablast; Tph help | Fast (weeks) | CD20/CD19 depletion, MMF, CYC, belimumab (partial) |
| Plasma-cell reservoir (new, proposed) | CD19⁻ BCMA⁺ long-lived plasma cells | Very slow | CD38, BCMA; not CD20 |
| Kidney activity | Infiltrating monocytes, disease macrophages, neutrophils | Fast | Steroids, MMF, CNI, complement |
| Kidney chronicity / damage | Lymphoid aggregates, GZMK⁺ CD8, parietal-cell sclerosis | Absorbing | Nothing reverses it |

- **Hypothetical.** Splitting a slow plasma-cell reservoir from the fast B-effector state is the cellular route by which LUMEN could represent refractoriness to CD20 depletion. It adds one parameter; it should be added only if a dataset can identify it (serology half-life after CD20 vs CD19 vs BCMA depletion). Follow-up queued as Q56.

## 7. Unknowns

- How much LN kidney antibody comes from local (interstitial germinal-centre) vs marrow plasma cells: *unknown*.
- Whether pDCs or LDGs are present and active in human LN glomeruli at scale: *unknown* in open data.
- Whether DN2 expansion predicts LN response (not just activity): no prospective human test found.

## Sources

PubMed, accessed 2026-11-23 (abstracts). Links inline above. Additional: [Werth 2022 LILAC](https://doi.org/10.1056/NEJMoa2118024); [CAR-T review, *Arthritis Rheumatol* 2024](https://doi.org/10.1002/art.42784); [anti-CD19 rationale, *Arthritis Res Ther*](https://doi.org/10.1186/ar3909); [low-dose IL-2 review, *Curr Opin Rheumatol* 2023](https://doi.org/10.1097/BOR.0000000000000924).
