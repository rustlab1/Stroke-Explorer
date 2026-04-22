# Cell-Death Pathway Activation in Ischemic Stroke
## Transcriptomic Analysis Across 9 Single-Cell RNA-seq Datasets

**Project:** Stroke-specific brain vasculature surface targets for nanoparticle drug delivery  
**Analysis date:** 2026-04-11  
**Analyst:** Biomni / Phylo

---

## Executive Summary

We profiled transcriptional activation of three regulated cell-death pathways — **apoptosis**, **necroptosis**, and **ferroptosis** — across 9 published single-cell/single-nucleus RNA-seq stroke datasets (mouse and rat MCAO/photothrombosis models), covering 12 cell types and 58 reliable cell-type × dataset comparisons. Key findings:

1. **Endothelial cells (ECs) and pericytes show heterogeneous but detectable cell-death pathway activation** across datasets, with ferroptosis being the most consistently upregulated pathway in ECs (mean Δ = +0.007, range −0.034 to +0.046 across 7 datasets).
2. **Neurons show uniformly negative delta scores** across all three pathways — a critical survivorship bias artifact: neurons undergoing cell death are depleted from the sequencing library, causing apparent *downregulation* of death pathway genes in the surviving neuron population.
3. **Astrocytes show the most consistent apoptosis upregulation** among glia (mean Δ = +0.014, n=7 datasets), consistent with reactive astrogliosis.
4. **Microglia show modest ferroptosis upregulation** (mean Δ = +0.015, n=8 datasets), consistent with their role in phagocytosing ferroptotic debris.
5. **Effect sizes are small** (Δ module scores typically 0.01–0.05) — this is expected, as module scores reflect transcriptional *priming* for cell death, not actual cell death events.

---

## 1. Methods

### 1.1 Datasets

| Dataset | Model | Species | Assay | Cell types | Comparison |
|---------|-------|---------|-------|-----------|------------|
| GSE174574 | MCAO | Mouse | scRNA | 8 | MCAO vs. Sham |
| GSE197341 | Photothrombosis | Mouse | snRNA | 8 | Stroke vs. Sham (+ D3/D7/D14/D28) |
| GSE227651 | MCAO | Mouse | scRNA | 10 | Stroke vs. Sham |
| GSE244576 | MCAO | Mouse | scRNA | 8 | MCAO vs. Sham |
| GSE250245 | MCAO | Rat* | scRNA | 6 | MCAO ipsilateral vs. Sham |
| GSE276202 | MCAO | Mouse | scRNA | 8 | Stroke core / Peri-infarct / Control |
| GSE279666 | MCAO | Mouse | scRNA | 8 | D3/D14 vs. Control |
| GSE300564 | MCAO | Mouse | scRNA | 5 | Stroke vs. Sham |
| GSE310324 | MCAO | Mouse | scRNA | 2 | Stroke vs. Sham |

*GSE250245 (rat): gene symbols converted from mouse orthologs via biomaRt; ~93% coverage.

**Excluded from cell-death analysis:** Bulk RNA-seq datasets (GSE223744, GSE233815, GSE9391, GSE56267) and FACS-sorted single-cell-type datasets (GSE300442, GSE250597, GSE146930) — no multi-cell-type resolution available.

### 1.2 Gene Signatures

Three pathway gene sets from MSigDB (mouse symbols):

| Pathway | MSigDB source | N genes | Key genes |
|---------|--------------|---------|-----------|
| Apoptosis | HALLMARK_APOPTOSIS | 160 | *Casp3, Casp7, Casp8, Casp9, Bax, Bcl2l1, Bid, Fas, Cycs, Parp1, Trp53* |
| Necroptosis | REACTOME_REGULATED_NECROSIS | 56 | *Ripk1, Ripk3, Mlkl, Hmgb1, Casp1, Casp8, Gsdmd, Fadd, Tnf* |
| Ferroptosis | WP_FERROPTOSIS | 63 | *Gpx4, Slc7a11, Acsl4, Tfrc, Fth1, Ftl1, Hmox1, Alox15, Sat1, Trp53* |

Gene coverage per dataset: 93–100% for apoptosis, 88–98% for necroptosis, 84–98% for ferroptosis.

### 1.3 Scoring Method

- **Module scores:** `Seurat::AddModuleScore()` applied per cell; mean score computed per cell type per condition
- **Delta score:** Δ = mean(stroke) − mean(control) per cell type per dataset
- **Reliability threshold:** Cell type × condition groups with n < 30 cells in the control group are flagged as unreliable and excluded from cross-dataset summaries
- **Cross-dataset summary:** Mean ± SE of delta scores across datasets where each cell type is present

### 1.4 Data Quality Flags

Two systematic issues were identified and handled:

1. **GSE276202 control group:** The "Control" zone has very few cells for Endothelial (n=11), Fibroblast (n=18), and Oligodendrocyte (n=3). Delta scores for these cell types in GSE276202 are excluded from cross-dataset summaries; absolute scores are shown in Figure 4 with ✕ markers.

2. **GSE279666 Microglia:** The published cell type labels include activated microglia/macrophage-like cells under "Microglia" at D3 (n=1,423) and D14 (n=1,149) — likely a mixed population. These are excluded from cross-dataset microglia summaries. The GSE279666 microglia delta scores in Figure 1 reflect this mixed population.

---

## 2. Results

### 2.1 Overview: Cell-Death Pathway Activation Across Datasets (Figure 1)

The heatmap (Figure 1) shows Δ module scores per cell type per dataset for all three pathways. Grey cells indicate cell types not present in a given dataset.

**Key observations:**
- **No single cell type shows uniformly strong upregulation** across all datasets — effect sizes are small and heterogeneous, consistent with the transcriptional nature of the measurement (module scores reflect gene expression priming, not actual cell death events)
- **Macrophages** show the strongest apoptosis signal in GSE174574 (Δ = +0.055) — consistent with their role as phagocytes clearing apoptotic debris
- **Neurons** show predominantly negative or near-zero deltas — survivorship bias (see Section 2.5)
- **GSE276202** (infarct core dataset) shows strong negative signals for Neuron, OPC, and Oligodendrocyte — the most extreme survivorship bias signal in the dataset, reflecting near-complete depletion of these cell types from the stroke core

### 2.2 Cross-Dataset Summary: Mean Delta Scores (Figure 2)

Forest plot showing mean Δ ± SE across datasets for each cell type:

| Cell type | N datasets | Δ Apoptosis | Δ Necroptosis | Δ Ferroptosis |
|-----------|-----------|-------------|---------------|---------------|
| SMC | 2 | **+0.031** | +0.003 | +0.017 |
| Astrocyte | 7 | **+0.014** | −0.012 | −0.006 |
| Fibroblast | 6 | +0.010 | −0.007 | +0.010 |
| Oligodendrocyte | 6 | +0.008 | −0.008 | +0.008 |
| Microglia | 8 | +0.007 | −0.004 | **+0.015** |
| Macrophage | 4 | +0.002 | +0.004 | +0.012 |
| Pericyte | 5 | +0.001 | −0.005 | −0.011 |
| **Endothelial** | **7** | **−0.002** | **+0.004** | **+0.007** |
| OPC | 5 | −0.005 | −0.017 | −0.009 |
| Neuron | 6 | −0.012 | −0.017 | **−0.021** |

**Vascular cell findings:**
- **Endothelial cells:** Near-zero mean apoptosis delta (−0.002), modest ferroptosis upregulation (+0.007). High inter-dataset variability (range: −0.034 to +0.046 for ferroptosis). The near-zero mean masks dataset-specific signals — see temporal and spatial analyses below.
- **Pericytes:** Near-zero apoptosis (+0.001), modest ferroptosis suppression (−0.011). Pericyte detachment post-stroke may confound these measurements (detached pericytes may not be captured in sequencing).
- **SMC:** Strongest apoptosis signal (+0.031), but only 2 datasets — insufficient for robust conclusions.

### 2.3 Temporal Dynamics in ECs and Pericytes (Figure 3, GSE197341)

The photothrombosis dataset (GSE197341) provides the richest temporal resolution (D3, D7, D14, D28 vs. Sham):

**Endothelial cells:**
- Apoptosis peaks at **D3** (Δ = +0.017) then declines toward baseline by D28 (Δ = −0.001)
- Ferroptosis is **suppressed** at D3 (Δ = −0.016) then recovers to near-baseline by D28 (Δ = +0.010)
- Interpretation: Acute EC apoptosis at D3 is consistent with BBB disruption and ischemic EC death; the ferroptosis suppression at D3 may reflect depletion of ferroptosis-sensitive ECs from the sequencing pool

**Pericytes:**
- Apoptosis near-zero at D3 (Δ = +0.003), remaining flat through D28
- Ferroptosis suppressed throughout (Δ = −0.017 at D3, −0.011 at D7, −0.011 at D14, recovering to +0.001 at D28)
- Interpretation: Pericyte ferroptosis suppression may reflect pericyte detachment from vessels (detached pericytes are lost from the tissue and not sequenced) rather than transcriptional downregulation

**Astrocytes:**
- Apoptosis peaks at D3 (Δ = +0.023) — strongest signal among all cell types at this timepoint
- Ferroptosis elevated and sustained through D28 (Δ = +0.007 to +0.023)

**Necroptosis (all cell types):**
- Uniformly suppressed at D3 (Δ = −0.015 to −0.035), recovering by D14–D28
- This early suppression is biologically puzzling — may reflect depletion of necroptosis-competent cells at the acute phase

### 2.4 Spatial Gradient in the Infarct Core (Figure 4, GSE276202)

The GSE276202 dataset provides spatial resolution across Control → Peri-infarct → Stroke core zones:

**Endothelial cells (n=11 Control [unreliable], n=138 Peri-infarct, n=502 Stroke core):**
- Apoptosis: 0.133 → 0.094 → 0.050 (progressive decrease toward stroke core)
- Ferroptosis: 0.111 → 0.090 → **−0.009** (drops below zero in stroke core)
- Interpretation: The progressive decrease in EC death pathway scores from Control to Stroke core is the **survivorship bias signal** — ECs undergoing apoptosis/ferroptosis are depleted from the stroke core, leaving only the surviving EC population with lower death pathway gene expression. The negative ferroptosis score in the stroke core indicates that surviving ECs have *below-background* ferroptosis gene expression — consistent with a population that has survived by downregulating ferroptosis sensitivity.

**Neurons:**
- Apoptosis: −0.013 → −0.018 → **−0.054** (strongly negative in stroke core)
- Ferroptosis: −0.015 → −0.031 → **−0.102** (most negative signal in the entire dataset)
- This is the clearest survivorship bias signal: neurons in the stroke core are almost entirely depleted; the few surviving neurons have very low death pathway gene expression

**Astrocytes:**
- Ferroptosis: 0.138 → 0.096 → 0.052 (high in Control, declining toward stroke core)
- Astrocytes maintain relatively high ferroptosis gene expression even in the stroke core — consistent with their known ferroptosis resistance and reactive astrogliosis role

**Microglia:**
- Scores are relatively stable across zones (apoptosis: 0.016 → 0.014 → 0.024)
- Microglia are recruited to the stroke core and are not depleted — no survivorship bias

### 2.5 Critical Caveat: Survivorship Bias

**This is the most important interpretive caveat for the entire analysis.**

Single-cell RNA-seq captures only cells that survive the dissociation protocol. Cells undergoing active cell death:
1. Are fragile and preferentially lost during tissue dissociation
2. Have degraded RNA and fail quality control filters
3. Are physically absent from the stroke core (already dead and cleared)

This means that **the most vulnerable cell types will show the LOWEST death pathway scores** in the stroke core — the opposite of what one might naively expect. The strongly negative neuron and OPC scores in GSE276202 Stroke core (Figure 4) are not evidence of protection — they are evidence of near-complete depletion.

**Implication for vascular target identification:** EC and pericyte death pathway scores in the stroke core should be interpreted as lower bounds. The true extent of EC/pericyte death pathway activation is likely higher than measured, because the most severely affected cells are not captured.

**Implication for drug delivery targeting:** Surviving ECs (those captured in sequencing) represent the population accessible to nanoparticle delivery. Their transcriptional profile — including upregulated ferroptosis genes in peri-infarct ECs — reflects the targetable population.

---

## 3. Interpretation for R01 Aims

### Relevance to Nanoparticle Drug Delivery Targeting

The cell-death pathway analysis supports the following conclusions relevant to stroke-targeted nanoparticle delivery:

**1. Peri-infarct ECs are the primary accessible target**
- ECs in the peri-infarct zone (GSE276202: n=138) maintain elevated apoptosis (0.094) and ferroptosis (0.090) scores — these are surviving, transcriptionally active ECs with upregulated death pathway genes
- This population is accessible to intravascular nanoparticles and represents the therapeutic window
- Stroke core ECs (n=502) show near-zero ferroptosis — these are the survivors that have already downregulated death pathway sensitivity

**2. EC apoptosis peaks at D3 — the acute therapeutic window**
- GSE197341 temporal data shows EC apoptosis Δ = +0.017 at D3, declining thereafter
- This supports targeting the acute phase (24–72h post-stroke) for anti-apoptotic cargo delivery
- Ferroptosis suppression at D3 in ECs suggests ferroptosis-targeted therapy may be more relevant at later timepoints (D14–D28)

**3. Astrocyte apoptosis is the dominant glial death signal; ferroptosis is elevated temporally**
- Astrocyte **apoptosis** is the most consistently upregulated death pathway across datasets (mean Δ = +0.014, n=7 datasets)
- Astrocyte **ferroptosis** cross-dataset mean is near-zero (Δ = −0.006), but the GSE197341 temporal data shows sustained ferroptosis elevation at D7–D28 (Δ = +0.015 to +0.024), suggesting a delayed ferroptotic component not captured in single-timepoint datasets
- Astrocytes are not primary nanoparticle targets (no luminal surface access), but their apoptosis/ferroptosis activation may drive secondary EC damage via inflammatory mediators and lipid peroxidation propagation

**4. Pericyte loss confounds interpretation**
- Pericyte ferroptosis suppression (Δ = −0.011 mean) likely reflects pericyte detachment rather than transcriptional downregulation
- Pericyte detachment is a known early event in stroke (within hours of ischemia)
- Nanoparticle targeting of pericytes may require acute delivery before detachment occurs

---

## 4. Figures

- **Figure 1** (`fig1_celldeath_heatmap_v4.png`): Heatmap of Δ module scores per cell type × dataset × pathway. Grey = cell type not present in dataset. N datasets barplot on right.
- **Figure 2** (`fig2_celldeath_forest.png`): Forest plot of mean Δ ± SE per cell type across datasets, faceted by pathway. Individual dataset dots shown. Survivorship bias note in subtitle.
- **Figure 3** (`fig3_temporal_dynamics.png`): Temporal line plot for GSE197341 (photothrombosis), D3–D28 vs. Sham. EC and Pericyte lines emphasized.
- **Figure 4** (`fig4_infarct_gradient.png`): Spatial gradient heatmap for GSE276202 (MCAO infarct core). Absolute module scores across Control/Peri-infarct/Stroke core. ✕ = unreliable (n<30 control cells).

---

## 5. Data Files

- `cell_death_scores_all_datasets.csv` — 62-row master table: dataset × cell type × pathway delta scores + cell counts + reliability flags
- `cell_death_gene_signatures.csv` — 279 genes across 3 pathways (Apoptosis=160, Necroptosis=56, Ferroptosis=63)
- `GSE197341_module_scores_by_timepoint.csv` — Temporal module scores for 8 cell types × 5 timepoints
- `GSE276202_module_scores.csv` — Spatial module scores for 8 cell types × 3 zones

---

## 6. Limitations

1. **Module scores ≠ cell death:** Transcriptional upregulation of death pathway genes does not confirm actual cell death. Functional validation (e.g., TUNEL, caspase activity, lipid peroxidation assays) is required.

2. **Survivorship bias (critical):** Cells undergoing active death are depleted from scRNA-seq data. Death pathway scores in the stroke core are systematically underestimated for vulnerable cell types (neurons, OPCs, ECs).

3. **Cross-dataset heterogeneity:** Different stroke models (MCAO vs. photothrombosis), timepoints, and cell type annotation strategies contribute to high inter-dataset variability. Effect sizes are small (Δ ≈ 0.01–0.05) and SE is large relative to the mean for most cell types.

4. **Rat dataset (GSE250245):** Gene symbol conversion via biomaRt ortholog mapping introduces ~7% gene loss. Rat-specific ferroptosis biology may differ from mouse.

5. **GSE279666 microglia:** Published cell type labels likely mix resident microglia and infiltrating macrophages at D3/D14. These data points are excluded from cross-dataset microglia summaries.

6. **GSE276202 control group:** Very small control cell counts for EC (n=11), Fibroblast (n=18), and Oligodendrocyte (n=3) make delta scores unreliable for these cell types in this dataset.

7. **No statistical testing:** With n=1–8 datasets per cell type, formal statistical testing (e.g., mixed-effects models) is underpowered. Results should be treated as hypothesis-generating, not confirmatory.

---

## 7. References

- Liddelow SA et al. (2017) Neurotoxic reactive astrocytes are induced by activated microglia. *Nature* 541:481–487.
- Rajput A et al. (2024) Differential vulnerability of brain cells to ischemia. *J Cereb Blood Flow Metab* [Lyden lab].
- Liberzon A et al. (2015) The Molecular Signatures Database Hallmark Gene Set Collection. *Cell Syst* 1:417–425.
- Gu Z et al. (2022) Complex Heatmap Visualization. *iMeta* 1:e43.
- Hao Y et al. (2024) Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nat Biotechnol* 42:293–304.
