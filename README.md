# 🫁 Lung Cancer Classifier

A bioinformatics and machine learning pipeline for identifying genetic biomarkers of lung adenocarcinoma (LUAD) using RNA-seq data from the Cancer Genome Atlas (TCGA).

---

## 📖 Overview

Lung cancer is a leading cause of cancer-related mortality worldwide. The lifetime risk of developing it is approximately 1 in 16 for men and 1 in 17 for women. This project combines **differential gene expression analysis** with **supervised machine learning** to identify potential genetic biomarkers and build classifiers capable of distinguishing cancerous tissue from healthy tissue.

The pipeline was applied to **598 RNA-seq samples** (539 primary tumors, 59 normal solid tissue controls) from the TCGA-LUAD dataset, identifying **515 differentially expressed genes (DEGs)** and training models that achieve up to **100% classification accuracy**.

> 📄 Full thesis available in [`Full_Thesis.pdf`](./Full_Thesis.pdf)

---

## ✨ Key Results

| Model | Accuracy |
|---|---|
| Random Forest (RF) | 94% |
| Support Vector Machine (SVM) | 100% |
| Logistic Regression (GLM) | — |

**Top candidate biomarkers identified:**
- `ALPP` — downregulated
- `BEX1` — downregulated
- `HSD17B6` — upregulated
- `TCEAL2` — upregulated

---

## 🔬 Pipeline Summary

The analysis is implemented in a single R script (`Lung-Cancer-Classifier.R`) and follows these main steps:

### 1. Data Loading & Preprocessing
- Reads TCGA-LUAD RNA-seq `.tsv` files from a local dataset directory
- Filters for **protein-coding genes** only
- Aggregates duplicate gene symbols by taking the max expression value
- Aligns sample IDs with phenotype metadata from the TCGA sample sheet
- Excludes recurrent tumor samples to avoid data leakage

### 2. Differential Expression Analysis (DESeq2)
- Constructs a `DESeqDataSet` from raw read counts
- Runs DESeq2 normalization and differential expression
- Applies significance thresholds: `padj < 0.05` and `|log2FoldChange| > 1`
- Labels genes as **UP**, **DOWN**, or **Normal**

### 3. Exploratory Visualization
- **Heatmap** of the top 100 DEGs (ComplexHeatmap), annotated by sample type
- **Volcano plot** showing fold-change vs. significance (ggplot2)
- **2D and 3D PCA** plots colored by sample type (plotly)

### 4. Gene Set Enrichment Analysis (GSEA)
- Functional enrichment using `clusterProfiler`
- Pathway visualization via `pathview`
- Notable enriched pathways: **neuroactive ligand-receptor interaction**, **thrombophilia-associated pathways**

### 5. Machine Learning Classification (caret)
- Feature matrix: normalized, scaled expression profiles of DEGs
- 70/30 train/test split with stratification
- 3-repeat cross-validation (`repeatedcv`)
- Models trained: **Random Forest**, **SVM (linear kernel)**, **Logistic Regression**
- Parallel training using `doParallel`

### 6. Model Evaluation
- Accuracy, Precision, Recall, and F1 Score per model
- Confusion matrices with fourfold plots
- ROC curves for all three models on the same plot
- Variable importance analysis to rank top predictive genes

---

## 🗂️ Repository Structure

```
Lung-cancer-classifier-main/
├── Lung-Cancer-Classifier.R   # Full analysis pipeline
├── Full_Thesis.pdf            # Bachelor's thesis document
└── README.md                  # This file
```

> **Note:** The TCGA dataset is not included in this repository due to size and licensing. See the Data section below for access instructions.

---

## 📦 Requirements

**R version:** 4.2.0+

### Required R Packages

| Category | Packages |
|---|---|
| Data I/O | `readr` |
| DE Analysis | `DESeq2` |
| Visualization | `ggplot2`, `ComplexHeatmap`, `plotly`, `scatterplot3d`, `rgl` |
| Enrichment | `clusterProfiler`, `pathview`, `GOplot`, `org.Hs.eg.db` |
| Machine Learning | `caret`, `caretEnsemble`, `randomForest`, `randomForestSRC`, `e1071` |
| Evaluation | `pROC` |
| Utilities | `dplyr`, `ggfortify`, `doParallel` |

Install all packages via Bioconductor and CRAN:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "ComplexHeatmap", "clusterProfiler",
                        "pathview", "GOplot", "org.Hs.eg.db"))

install.packages(c("readr", "ggplot2", "plotly", "scatterplot3d", "rgl",
                   "dplyr", "ggfortify", "caret", "caretEnsemble",
                   "randomForest", "randomForestSRC", "pROC",
                   "e1071", "doParallel"))
```

---

## 🗄️ Data

This project uses publicly available data from the **Cancer Genome Atlas (TCGA)**:

- **Dataset:** TCGA-LUAD (Lung Adenocarcinoma)
- **Data type:** RNA-seq (raw read counts, TPM)
- **Access:** [GDC Data Portal](https://portal.gdc.cancer.gov/)

After downloading, organize the data as follows:

```
Dataset/
├── TCGA-LUAD-RNA_SEQ/      # Folder containing per-sample .tsv files
└── Sample_Sheet.tsv         # TCGA GDC sample sheet
```

Update the paths in the script if needed:

```r
data.path  <- "Dataset/TCGA-LUAD-RNA_SEQ"
pheno.path <- "Dataset/Sample_Sheet.tsv"
```

---

## 🚀 Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/youssefygj/Lung-cancer-classifier.git
   cd Lung-cancer-classifier
   ```

2. Download and organize the TCGA-LUAD dataset as described above.

3. Open `Lung-Cancer-Classifier.R` in RStudio or run it from the terminal:
   ```bash
   Rscript Lung-Cancer-Classifier.R
   ```

4. Key output files generated:
   - `res.txt` — Full DESeq2 results table
   - `res.degs.txt` — List of significant DEG gene names
   - `final.csv` — Annotated DEG results with fold-change labels
   - `dataForTraining.csv` — Processed feature matrix used for ML
   - `data.RDATA` — Saved R workspace with all major objects

---

## 📊 Sample Outputs

- Heatmap of top 100 differentially expressed genes
- Volcano plot (UP/DOWN/Normal gene labeling)
- Interactive 3D PCA (tumor vs. normal separation)
- ROC curves comparing RF, SVM, and Logistic Regression
- Grouped bar charts for Accuracy, Precision, Recall, and F1 Score

---

## 👤 Author

**Youssef Fahim**  
Bachelor's Project — 2023

---

## 📜 License

This project is for academic and research purposes. The TCGA data used in this analysis is subject to the [NIH GDC Data Access Policies](https://gdc.cancer.gov/access-data/data-access-policies).
