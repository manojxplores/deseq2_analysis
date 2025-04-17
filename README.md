# deseq2_analysis

**Title:** 🧬RNA-Seq Differential Expression Analysis of HD-Zip Transcription Factors in Soybean under Dehydration and Salt Stress

## 📌 Objective

To identify and characterize differentially expressed HD-Zip transcription factor genes in *Glycine max* (soybean) root samples under dehydration and salt stress conditions using RNA-Seq data.

---

## 📁 Dataset Overview

**Source:** NCBI GEO  
**Accession:** GSE57252 
**Experiment Type:** Expression profiling by high throughput sequencing ✅  
**Organism:** *Glycine max* (soybean)  
**Conditions:**  
- Control (0 hr)  
- Dehydration: 1 hr, 6 hr, 12 hr  
- Salt stress: 1 hr, 6 hr, 12 hr  
**Replicates:** 3 biological replicates per condition

---

## 🧪 Methods

### 🔹 Tools & Libraries
- R
- DESeq2
- RUVSeq
- pheatmap
- ggplot2
- EnhancedVolcano

### 🔹 Analysis Steps
1. **Preprocessing**: Import count data and metadata
2. **Normalization**: DESeq2 and optional RUVg batch correction
3. **Differential Expression Analysis**:
   - Comparisons:
     - Dehydration (1/6/12 hr) vs Control
     - Salt stress (1/6/12 hr) vs Control
4. **Filtering**: Genes with |log2FC| ≥ 1 and adjusted p-value < 0.1
5. **Visualization**: PCA, heatmap, volcano plot
6. **Interpretation**: Identify stress-responsive HD-Zip genes

---

## 📂 Repository Structure

```bash
.
├── datasets/                
├── results/             
├── scripts/  
└── README.md
