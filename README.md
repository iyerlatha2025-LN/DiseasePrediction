
# Disease Prediction — Clinical AI Portfolio

**Author:** Latha Iyer | SAI Strategic Agentic Intelligence | Louisville, KY
**ORCID:** 0009-0000-8755-8805 | **Zenodo:** 10.5281/zenodo.14837564

---

## Repository Contents

### 1. ALS Biomarker Integration Framework (Python)
`ALS_Framework_V1_0_Latha.ipynb`

Multi-source clinical biomarker analysis for ALS progression prediction.

- 369 reconstructed patient profiles across two published cohorts
  (Lu et al. 2015: 219 patients · Verde et al. 2019: 150 patients)
- Mitochondrial pathway modelling — NAD+, ATP, ROS, VAP biomarkers
- Machine learning models: Random Forest, Gradient Boosting, Logistic Regression
- Cross-pathway mitochondrial target AUC: 0.88
- SHAP analysis identifying key molecular and clinical markers
- Conformal prediction safety wrapper (ACCB) — HALT / REVIEW / PROCEED routing
- Comprehensive sensitivity analysis (6 types)
- Data transparency: profiles computationally reconstructed from published
  cohort parameters — not raw clinical records

---

### 2. CVD Genomic Risk Prediction (R)
`cvd_genomic_risk_prediction_enhancedvf.R`

Genomic risk prediction integrating SNP data and clinical features
for cardiovascular disease classification.

- 46,218 patients with demographic and genomic SNP features
- SNPs analysed: rs10757278, rs8055236, rs4665058, rs1333049
- K-Medoids clustering identifying genetically distinct patient subgroups
- Key finding: Cluster 5 carrying 35% CVD rate vs 1-3% in others
- Gradient Boosting Model AUC: 0.913
- Implemented in R using tidyverse, ggplot2, and caret

---

## Frameworks Built on This Research

**ACCB — Adaptive Confidence Circuit Breaker**
Conformal prediction safety wrapper providing mathematically guaranteed
confidence bounds on clinical AI outputs. Validated on ALS data.
Published: zenodo.org/records/14837564

**CCIM — Context Chain Integrity Monitor**
Preventive AI safety framework verifying reasoning at every stage
of the inference chain — ensuring wrong outputs cannot form rather
than detecting them after the fact.

---

## IEEE Women's Chapter — June 2026
Invited speaker: AI Safety in High-Stakes Clinical Systems
