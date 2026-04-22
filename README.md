# Dengue-Transcriptomics-XAI
An Explainable AI approach for Rapid Multiplex RT-PCR Assay Design in Dengue Haemorrhagic Fever. 
# Translating Dengue Transcriptomics into Clinical Diagnostics

![GitHub last commit](https://img.shields.io/github/last-commit/A7MSDS/Dengue-Transcriptomics-XAI)
![R](https://img.shields.io/badge/Language-R-blue)
![Status](https://img.shields.io/badge/Status-Completed-success)

This repository contains the bioinformatics pipeline, Explainable AI (XAI) models, and primer design data for the master's thesis project: **"Translating Dengue Transcriptomics into Clinical Diagnostics: An Explainable AI Approach for Rapid Multiplex RT-PCR Assay Design."**

## 📌 Project Overview
Dengue Hemorrhagic Fever (DHF) is a severe, life-threatening complication of dengue infection. This project leverages Machine Learning (Random Forest) to classify patient samples based on transcriptomic data (GSE51808). Utilizing the Mean Decrease Gini index as a Global Explainer (XAI), we identified top predictive biomarkers (**NOTCH4, PHTF1, and YIPF1**). These findings were translated into a clinical setting by designing a strict, single-tube Multiplex RT-PCR assay.

## 📂 Repository Structure
* `02_Scripts/`: R scripts for data preprocessing, Random Forest training, and Gini Index extraction.
* `03_Clean_Data/`: Processed and normalized transcriptomic data ready for ML modeling.
* `04_Results/`: Output figures including Heatmaps and XAI Global Importance plots.
* `Multiplex_Assay/`: Primer sequences and PCR parameters for the designed multiplex assay.

## 🧬 Highlighted Biomarkers (Multiplex Targets)
1. **NOTCH4** (Product Size: 100 bp)
2. **PHTF1** (Product Size: 134 bp)
3. **YIPF1** (Product Size: 189 bp)
*Note: All primers are optimized for a unified melting temperature ($T_m$) of 60.0°C ± 0.4°C.*

## ⚙️ How to Reproduce
1. Clone this repository: `git clone https://github.com/A7MSDS/Dengue-Transcriptomics-XAI.git`
2. Open the `.Rproj` file in RStudio.
3. Run the scripts sequentially starting from `01_Data_Prep.R`.

## ✉️ Contact
For questions regarding the methodology or data, please refer to the main thesis document or contact the author.
