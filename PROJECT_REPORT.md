
# MHC-I Peptide-MHC Binding Prediction Using SVM and BLOSUM62 Encoding

**Author:** T.D.Varadhavallabha  
**Institution:** Sri Venkateswara College of Engineering  
**Date:** July 2026  
**Repository:** https://github.com/varadha0505-max/mhc-binding-predictor

---

## Abstract

Developed a supervised machine learning model to predict peptide binding 
to the MHC Class I molecule HLA-A*02:01, a key step in neoantigen-based 
personalised cancer vaccine design. Using 8,346 experimentally validated 
peptide-MHC binding measurements from the Immune Epitope Database (IEDB), 
we trained a Support Vector Machine (SVM) with an RBF kernel on 
BLOSUM62-encoded 9-mer peptide features. The model achieved ROC-AUC 0.949 
and 88.6% accuracy on a held-out test set, exceeding our target of AUC > 0.85 
and approaching the performance of published tools such as NetMHCpan. 
Permutation-based feature importance analysis revealed that position P2 
carries the highest predictive weight, and that Leucine at P2 is the single 
most important amino acid feature — consistent with the known HLA-A*02:01 
anchor residue preference (P2=Leu/Met) established by crystallographic 
and biochemical studies. This biological finding emerged from sequence data 
alone, without any prior structural supervision.

---

## 1. Introduction

### 1.1 Biological Background

Cancer immunotherapy exploits the ability of CD8+ cytotoxic T-lymphocytes 
(CTLs) to recognise and kill tumour cells. Recognition depends on a 
three-component interaction: a peptide fragment derived from an intracellular 
protein, presented by an MHC Class I molecule on the cell surface, detected 
by a T-cell receptor. Tumour-specific somatic mutations can generate novel 
peptide sequences (neoantigens) not present in the normal proteome. If these 
neoantigens bind MHC-I with sufficient affinity (IC50 < 500 nM by convention), 
they can trigger a tumour-specific immune response.

Personalised cancer vaccines (currently in Phase III trials by Moderna/BioNTech 
for melanoma) work by identifying patient-specific neoantigens from tumour 
sequencing data and synthesising mRNA encoding those peptides. The bottleneck 
is computationally predicting which of thousands of candidate peptides will 
actually bind MHC — wet-lab validation of every candidate is infeasible.

### 1.2 HLA-A*02:01

HLA-A*02:01 is the most common MHC Class I allele in humans (~45% population 
frequency). It has the most curated binding data in IEDB and is the primary 
allele studied in neoantigen vaccine trials. It has known anchor preferences: 
hydrophobic residues (Leu, Met) at position 2, and hydrophobic residues 
(Val, Leu) at position 9 of 9-mer peptides.

### 1.3 Problem Statement

Given a 9-amino acid peptide sequence, predict whether it binds 
HLA-A*02:01 with IC50 < 500 nM (binary classification).

---

## 2. Data

### 2.1 Source
Immune Epitope Database (IEDB), queried for:
- MHC Restriction: HLA-A*02:01
- Assay type: IC50 (nM)
- Organism: Homo sapiens
- Downloaded: July 2026

### 2.2 Cleaning Pipeline

| Step | Action | Rows remaining |
|------|--------|----------------|
| Raw download | — | 30,146 |
| Filter nM IC50 assays | Remove non-IC50 assay types | 17,096 |
| Drop missing values | Remove NaN peptide or IC50 | 17,075 |
| Keep 9-mers | Filter by peptide length == 9 | 11,723 |
| Standard amino acids | Remove non-canonical residues | 11,723 |
| Deduplicate | Median IC50 per unique peptide | 8,346 |

### 2.3 Label Definition
- Binder: IC50 < 500 nM (positive class, label=1) — 3,112 peptides (37.3%)
- Non-binder: IC50 ≥ 500 nM (negative class, label=0) — 5,234 peptides (62.7%)

The 500 nM threshold is the standard used in NetMHCpan and the majority 
of MHC binding prediction literature.

---

## 3. Methods

### 3.1 Feature Engineering

**One-Hot Encoding (OHE):**  
Each peptide is represented as a 9×20 binary matrix where position (i,j)=1 
if amino acid j occupies position i. Flattened to a 180-dimensional vector. 
Treats all amino acids as equally dissimilar.

**BLOSUM62 Encoding:**  
Each amino acid is replaced by its row in the BLOSUM62 substitution matrix 
— a 20-dimensional vector of log-odds substitution scores derived from 
evolutionarily conserved protein blocks. This encodes physicochemical 
similarity: amino acids with similar biochemical properties receive similar 
vectors. Each 9-mer becomes a 9×20 matrix (180-dimensional vector after 
flattening), normalised to zero mean and unit variance before model training.

BLOSUM62 was chosen over OHE because MHC binding pockets select residues 
partly by size and hydrophobicity — properties that BLOSUM62 implicitly 
encodes through evolutionary conservation patterns.

### 3.2 Model

Support Vector Machine (SVM) with Radial Basis Function (RBF) kernel, 
implemented via scikit-learn's CalibratedClassifierCV wrapper for 
probability estimates.

Hyperparameters selected by 5-fold cross-validated grid search over:
- C ∈ {0.1, 1.0, 10.0, 100.0}
- γ ∈ {scale, 0.001, 0.01}

Optimal: C=10, γ=0.001.

Class imbalance handled by class_weight='balanced'.

### 3.3 Evaluation
- Primary metric: ROC-AUC (threshold-independent, appropriate for imbalanced data)
- Secondary: Accuracy, Precision, Recall, F1, Average Precision
- Train/test split: 80/20 stratified
- Cross-validation: 5-fold stratified

### 3.4 Interpretability
Permutation importance: each of the 180 features was randomly shuffled 
10 times; the mean decrease in ROC-AUC was recorded. Features were 
reshaped to a 9×20 matrix (position × amino acid) for biological 
interpretation.

---

## 4. Results

### 4.1 Model Performance

| Metric | Value |
|--------|-------|
| ROC-AUC | **0.9486** |
| Average Precision | 0.8968 |
| Accuracy | 88.6% |
| Non-binder Precision | 0.92 |
| Non-binder Recall | 0.90 |
| Binder Precision | 0.83 |
| Binder Recall | 0.87 |
| CV AUC (5-fold) | 0.934 ± 0.009 |

### 4.2 Biological Validation — Anchor Residues

Permutation importance analysis revealed:
- P2 is the single most important peptide position
- Within P2, Leucine (L) has the highest importance score of any 
  amino acid at any position in the entire 9×20 feature matrix
- P9 shows elevated importance relative to central positions

This precisely matches the experimentally determined anchor residue 
preferences of HLA-A*02:01: P2=Leucine/Methionine (primary anchor), 
P9=Valine/Leucine (secondary anchor). The model recovered this 
structural constraint from IC50 measurements alone.

### 4.3 Neoantigen Panel

| Peptide | Source | Score | Result |
|---------|--------|-------|--------|
| GILGFVFTL | Influenza M1 | 0.943 | ✓ |
| NLVPMVATV | CMV pp65 | 0.929 | ✓ |
| RMFPNAPYL | HPV / Cervical CA | 0.881 | ✓ |
| KILDGVQKL | Melanoma BRAF V600E | 0.840 | ✓ |
| YLQLQHFNL | Melanoma NRAS Q61R | 0.840 | ✓ |
| FLTPKKLQC | Negative control | 0.041 | ✓ |
| DDDPGIARI | Negative control | 0.002 | ✓ |
| GGGPGGGGG | Negative control | 0.000 | ✓ |
| KLGGALQAK | Influenza M1 (weak) | 0.053 | ✗* |

*KLGGALQAK is a documented weak binder with IC50 ~2,000 nM, above the 
500 nM classification threshold. The model's low score is consistent 
with this borderline status.

Panel accuracy: 8/9 (89%)

---

## 5. Discussion

### 5.1 Performance
AUC 0.949 substantially exceeds the project target of 0.85 and is 
competitive with published tools. NetMHCpan 4.1 reports AUC ~0.96 
on benchmark datasets, trained on significantly more data and using 
neural network architectures. The gap is consistent with our simpler 
model and smaller training set.

### 5.2 Why BLOSUM62 Works
The MHC binding groove selects peptides based on the physicochemical 
complementarity between residue side chains and groove pockets. 
BLOSUM62 encodes this implicitly: amino acids that substitute for 
each other in evolution tend to be biochemically similar. OHE treats 
all substitutions as equally wrong; BLOSUM62 treats similar amino acids 
as similar features, which is biologically appropriate.

### 5.3 Limitations
1. Allele-specific: model only applies to HLA-A*02:01
2. Binary threshold (500 nM) loses affinity information; regression 
   on log-IC50 would be more informative
3. No structural features (peptide flexibility, MHC groove geometry)
4. KLGGALQAK miss highlights boundary cases near the threshold

### 5.4 Future Directions
1. Pan-allele model: concatenate HLA pseudo-sequence with peptide encoding
2. Regression target: predict log-IC50 directly
3. Deep learning: 1D-CNN or Transformer encoder for positional motif learning
4. Application: screen TCGA somatic mutations to rank candidate neoantigens
5. Deployment: Streamlit web app for clinical/research use

---

## 6. Conclusion

Built a complete ML pipeline for MHC-I peptide binding prediction 
achieving AUC 0.949 on 8,346 HLA-A*02:01 peptides. The model's feature 
importance profile independently recovered the known P2 Leucine anchor 
preference — a structural constraint established by decades of biochemical 
and crystallographic research — demonstrating that BLOSUM62-encoded SVM 
models can extract meaningful biological signal from binding affinity data. 
This validates the approach as a foundation for more sophisticated 
neoantigen prioritisation pipelines.

---

## References

1. Jurtz V et al. NetMHCpan-4.0: Improved Peptide–MHC Class I Interaction 
   Predictions Integrating Eluted Ligand and Peptide Binding Affinity Data. 
   J Immunol. 2017.
2. Henikoff S, Henikoff JG. Amino acid substitution matrices from protein 
   blocks. PNAS. 1992.
3. Vita R et al. The Immune Epitope Database (IEDB): 2018 update. 
   Nucleic Acids Res. 2019.
4. Moderna/BioNTech mRNA-4157 Phase III trial. ClinicalTrials.gov NCT03897881.
