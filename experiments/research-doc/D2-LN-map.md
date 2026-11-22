# D2: Lupus nephritis map

Status: first pass, 2026-11-22. Built from 20 primary, cohort or consensus papers found through PubMed (abstracts only). Each claim carries an evidence label. Mouse findings are marked as mouse findings.

## 1. Who gets LN, and what it costs

- **Established.** About 38% of people with SLE develop LN. In the SLICC inception cohort (n=1,827), 700 developed LN; 81% of them already had it at enrolment, within 15 months of diagnosis (*Rheumatology* 2016).
  - LN was more common in younger patients, men, and people of African, Asian and Hispanic ancestry.
  - 10-year risk of kidney failure (ESKD) was 10.1% with LN vs 4.3% overall; LN tripled the risk of death (HR 2.98).
- **Established.** Kidney-failure risk stopped improving in the mid-1990s. In developed countries the 5-year ESKD risk fell from 16% (1970s) to 11% (mid-1990s) and then plateaued; class IV carried a 15-year risk of 44% in the 2000s (meta-analysis, 187 studies, 18,309 patients, *Arthritis Rheumatol* 2016).
  - Consequence for LUMEN: the "standard of care" baseline in the trial harness has not moved in 30 years. A new drug has to beat a stable target.

## 2. The coordinate system: biopsy class plus two indices

- **Established.** LN is classified by biopsy (ISN/RPS). The 2018 revision kept classes I–VI, dropped the IV-S / IV-G split, and replaced the old "active / chronic" labels with an NIH **activity index** (0–24) and **chronicity index** (0–12) to be scored in every class (*Kidney Int* 2018).
  - Classes in brief: I–II mesangial; III focal and IV diffuse proliferative (endocapillary); V membranous (subepithelial); VI advanced sclerosis (≥90% sclerosed glomeruli).
- **Supported.** The two indices behave as two separate axes. In 194 biopsies from 2 cohorts, a factor analysis of all lesions recovered two main factors that match activity and chronicity. Endocapillary hypercellularity was in >70% of biopsies and travelled with neutrophils and cellular crescents (*KI Reports* 2026).
  - This is the same activity / damage split the six-state model uses. It is the strongest external support so far that the split is real, not an artefact of the model.
- **Supported.** Rare patterns sit outside the immune-complex picture: lupus podocytopathy (minimal-change or FSGS pattern, no capillary-wall deposits, diffuse foot-process effacement) is about 1% of LN biopsies (*Adv Chronic Kidney Dis* 2019). Thrombotic microangiopathy is another (section 3b).

## 3. The causal chain, compartment by compartment

### 3a. Antibody and antigen reach the kidney

- **Established.** LN starts with immune complexes in the glomerulus. Anti-dsDNA and anti-nucleosome IgG bound to chromatin fragments are the classic deposits.
- **Supported (mouse).** In NZB/W mice, chromatin–IgG complexes first appear in the mesangium together with anti-dsDNA and mild or silent nephritis. Only later does renal DNase1 fall; large chromatin fragments then build up in the glomerular basement membrane, and that marks severe disease (*PLoS One* 2009).
  - Reading: antibody starts the disease; failure to clear chromatin locally may turn mild mesangial disease into capillary-wall disease. Whether human renal DNase1 falls the same way is *unknown*.
- **Supported (mouse, with human association).** Anti-C1q antibodies are strongly associated with nephritis in SLE. In mice they deposit in glomeruli but cause disease only when C1q-fixing immune complexes are already there, and the damage needs both complement and Fcγ receptors (*J Clin Invest* 2004).
  - This fits D3: high anti-C1q predicts *better* response. Anti-C1q marks antibody-and-complement-driven disease, the kind standard therapy treats well.
- **Supported.** Fcγ receptor genetics change the risk. Low *FCGR3B* copy number roughly doubles LN risk (OR about 2.0; meta-analysis of 4 comparisons, 689 LN patients) (*Int J Rheum Dis* 2015). FcγRIIIb sits on neutrophils and helps clear immune complexes.

### 3b. Complement and the first responders

- **Established.** Low serum C3 and C4 track active LN, and glomerular C3 deposits are typical of proliferative disease.
- **Emerging (mouse plus human biopsy).** C5a acting on C5aR1 drives podocyte mitochondrial fission and injury. C5aR1 is raised in podocytes in human LN biopsies; blocking C5aR in lupus mice reduced podocyte damage and proteinuria (*Mol Ther* 2024).
- **Supported (small human cohort).** LN with thrombotic microangiopathy shows stronger alternative-pathway and terminal complement activation in urine (C3a, C5a, Ba, C5b-9) than LN alone, and these fall after treatment (*Clin Rheumatol* 2021, n=20 per group).
- **Supported.** Neutrophils are captured in glomerular capillaries through their Fcγ receptors within minutes of immune-complex deposition (intravital microscopy, mouse). The older view that poor NET clearance drives LN is challenged by genetic mouse data (review, *Curr Opin Rheumatol* 2019). Contested.

### 3c. Myeloid cells: the main intrarenal effector early on

- **Supported.** Kidney single-cell maps (AMP Phase 1) found 21 leukocyte subsets, a stepwise monocyte-to-macrophage differentiation inside the kidney, an interferon response in most cells, and close agreement between kidney and urine immune cells (*Nat Immunol* 2019).
- **Supported.** In the larger AMP data set (156 LN patients, 30 controls, scRNA-seq plus spatial), the immune response changes with damage (*Ann Rheum Dis* 2026):
  - early and active: dominated by monocytes and disease-specific macrophages (DMacs); in proliferative LN, DMac expansion scales with activity;
  - as the chronicity index rises: DMacs contract, and the infiltrate shifts to CD4 T cells, GZMK⁺ CD8 T cells, B cells and dendritic cells, **with a falling interferon response**;
  - DMacs come from both infiltrating monocytes and resident macrophages, and appear to signal to parietal epithelial cells, promoting glomerulosclerosis.
- **Supported (cross-species).** Disease-associated monocytes and macrophages expressing *TREM2, SPP1, APOE, CD9* are conserved across four mouse models and 155 LN patients and link to histological class and progression (*J Exp Med* 2025). This is a reason to trust mouse myeloid experiments more than mouse lymphoid ones.

### 3d. Adaptive immunity inside the kidney

- **Supported.** More than half of 68 LN biopsies had organised T:B aggregates or true germinal centres in the interstitium, with clonal B-cell expansion and ongoing somatic hypermutation. These went with immune complexes on tubular basement membranes (*J Immunol* 2011).
  - So the interstitium can make its own antibody. This is the tissue-resident adaptive axis in D3.
- **Emerging.** Kidney B cells in active LN are mostly atypical (age-associated) B cells and antibody-secreting cells, clonally expanded, IgG1 and IgG3, with fewer mutations than blood ASCs. This points to an extrafollicular response. GZMK⁺ CD8 T cells are the dominant kidney T cells, sit with atypical B cells in tertiary lymphoid structures, and may drive them through IFN-γ and IL-21 (*Ann Rheum Dis* 2025).
- **Emerging.** Kidney-infiltrating CD8 T cells (studied through urine) are recruited from circulating CX3CR1⁺ cells. Their expanded clones did **not** recognise the patient's own tubular cells; about 20% were EBV- or CMV-specific. They were activated instead by local IL-15 and IFN-β (*Sci Transl Med* 2026).
  - Reading: much of the T-cell infiltrate may be bystander, cytokine-driven inflammation, not kidney-specific autoimmunity. That would weaken the case for T-cell-specific therapies and strengthen the case for cytokine blockade.

### 3e. The kidney's own cells respond

- **Supported.** Tubular cells carry an interferon-response signature that distinguishes LN from controls; high tubular IFN and fibrosis signatures each went with failure to respond (*Nat Immunol* 2019, Der et al.).
- **Emerging (mouse plus 9 patients).** Podocytes in LN express SGLT2; empagliflozin reduced podocyte injury in MRL/lpr mice, and proteinuria fell in 9 treated patients (*Ann Rheum Dis* 2023). Uncontrolled.

### 3f. Damage and the genetic floor

- **Established.** Chronicity (glomerulosclerosis, fibrous crescents, tubular atrophy, interstitial fibrosis) predicts kidney failure (D3).
- **Supported.** *APOL1* high-risk genotypes (two G1/G2 alleles) raise the risk of LN-related ESKD in people of African ancestry. A genome-wide study of 901 LN-ESKD cases found no strong second gene beyond *APOL1*, suggesting environmental triggers such as interferon matter more (*Kidney Int* 2018; review *Am J Kidney Dis* 2018). *APOL1* expression is itself interferon-inducible (*Lupus* 2023).
  - Hypothetical link: in high-risk *APOL1* carriers, the lupus interferon set point (round 89) may act on the damage axis directly, not only through immune cells.

## 4. The map in one table

| Stage | Where | What happens | Readout | Evidence |
|---|---|---|---|---|
| 1. Antibody arrives | Mesangium | Chromatin–IgG complexes deposit | Anti-dsDNA, low C3/C4, class I–II | Established |
| 2. Capillary-wall disease | GBM, subendothelium | Large chromatin fragments; anti-C1q amplifies | Class III/IV, activity index | Supported (mouse mechanism) |
| 3. Innate effectors | Glomerulus | Complement (C5a), FcγR-captured neutrophils, monocytes become DMacs | Endocapillary hypercellularity, urine CD163/CD206 | Supported |
| 4. Tissue-resident adaptive | Interstitium | TLS / germinal centres, atypical B cells, ASCs, GZMK⁺ CD8 | Interstitial infiltrate, TBM deposits | Supported / Emerging |
| 5. Parenchymal response | Tubules, podocytes | IFN response, C5aR1 signalling, stress | Tubular IFN score, proteinuria | Supported / Emerging |
| 6. Scar | Whole kidney | DMac–parietal epithelial signalling, sclerosis, IF/TA; myeloid → lymphoid shift; IFN falls | Chronicity index, eGFR | Supported |
| Modifier | Germline | *APOL1*, *FCGR3B* copy number | Ancestry, genotype | Supported |

## 5. What this means for LUMEN

- **Supported.** Biopsy activity and chronicity are two separable axes in human data (section 2). The six-state model's activity / damage split matches the pathology.
- **Emerging, testable.** Inside the kidney, the interferon response *falls* as chronicity rises (AMP 2026). So kidney IFN is not a simple activity gauge: high in early active disease, lower in scarred disease. Prediction for GSE200306 (serial biopsies, round 91): baseline glomerular IFN score correlates negatively with chronicity index. Not yet pre-registered.
- **Hypothetical.** Stage 4 (tissue-resident adaptive immunity) is the one compartment that systemic B-cell depletion may not fully reach. If true, it is where CAR-T and anti-CD20 should differ. No post-CAR-T repeat biopsy series exists to test it (D3 gap).
- **Hypothetical.** If much of the kidney CD8 infiltrate is cytokine-driven bystander inflammation (Sci Transl Med 2026), the model should drive CD8 activity from IFN and IL-15 levels, not from an autoantigen-specific clone.

## 6. Gaps

- Human renal DNase1 and chromatin clearance across classes: *unknown*.
- No study measures all six stages in the same patients over time.
- Membranous (class V) disease has a different, podocyte-facing mechanism and is under-covered here; it needs its own pass.
- Most single-cell work is cross-sectional; the stage order in section 4 is inferred from chronicity, not observed in time.

## Sources (PubMed)

- [Frequency and outcome of LN, SLICC inception cohort, *Rheumatology* 2016](https://doi.org/10.1093/rheumatology/kev311)
- [ESRD risk in LN 1971–2015, meta-analysis, *Arthritis Rheumatol* 2016](https://doi.org/10.1002/art.39594)
- [ISN/RPS classification revision, *Kidney Int* 2018](https://doi.org/10.1016/j.kint.2017.11.023)
- [NIH activity and chronicity indices in 2 cohorts, *KI Reports* 2026](https://doi.org/10.1016/j.ekir.2025.11.011)
- [Lupus podocytopathy overview, *Adv Chronic Kidney Dis* 2019](https://doi.org/10.1053/j.ackd.2019.08.011)
- [Anti-dsDNA and renal DNase1 loss, NZB/W mice, *PLoS One* 2009](https://doi.org/10.1371/journal.pone.0008474)
- [Anti-C1q pathogenic only with glomerular C1q complexes, *J Clin Invest* 2004](https://doi.org/10.1172/JCI21075)
- [FCGR3B copy number loss and LN, meta-analysis, *Int J Rheum Dis* 2015](https://doi.org/10.1111/1756-185X.12342)
- [C5a–C5aR1 and podocyte mitochondrial fission, *Mol Ther* 2024](https://doi.org/10.1016/j.ymthe.2024.03.003)
- [Alternative complement in LN with TMA, *Clin Rheumatol* 2021](https://doi.org/10.1007/s10067-020-05499-1)
- [Neutrophils in LN, review, *Curr Opin Rheumatol* 2019](https://doi.org/10.1097/BOR.0000000000000577)
- [Immune cell landscape in LN kidneys, AMP, *Nat Immunol* 2019](https://doi.org/10.1038/s41590-019-0398-x)
- [Deep myeloid profiling of LN kidneys, AMP, *Ann Rheum Dis* 2026](https://doi.org/10.1016/j.ard.2026.03.004)
- [Human–mouse myeloid atlas, *J Exp Med* 2025](https://doi.org/10.1084/jem.20241873)
- [In situ B-cell responses and tubulointerstitial inflammation, *J Immunol* 2011](https://doi.org/10.4049/jimmunol.1001983)
- [Extrafollicular B cells and GZMK⁺ CD8 T cells in LN kidney, *Ann Rheum Dis* 2025](https://doi.org/10.1136/ard-2024-225876)
- [Cytokine-activated kidney CD8 T cells, *Sci Transl Med* 2026](https://doi.org/10.1126/scitranslmed.adz2015)
- [Tubular cell IFN and fibrosis signatures, *Nat Immunol* 2019](https://doi.org/10.1038/s41590-019-0386-1)
- [SGLT2 inhibition and podocytes, *Ann Rheum Dis* 2023](https://doi.org/10.1136/ard-2023-224242)
- [APOL1 GWAS in LN-ESRD, *Kidney Int* 2018](https://doi.org/10.1016/j.kint.2018.03.017)
- [APOL1-associated nephropathy review, *Am J Kidney Dis* 2018](https://doi.org/10.1053/j.ajkd.2018.06.020)
- [APOL1 and renal outcomes in SLE, *Lupus* 2023](https://doi.org/10.1177/09612033231172660)
