# ABCRpred: Prediction of Antibiotic Resistant Strains from Beta-Lactamase Proteins

Welcome to the official repository for **ABCRpred**, a machine learning-based computational tool developed for predicting ceftazidime-resistant and ceftazidime-sensitive bacterial strains using beta-lactamase protein sequences.

Web Server: https://webs.iiitd.edu.in/raghava/abcrpred/

---

## Citation

Maryam, L., Dhall, A., Patiyal, S., Usmani, S. S., Sharma, N., & Raghava, G. P. S. (2021).

**Prediction of antibiotic resistant strains of bacteria from their beta-lactamases protein**

bioRxiv.

https://doi.org/10.1101/2021.06.26.450028

---

## About the Study

Antimicrobial resistance (AMR) is one of the major global public health threats. Beta-lactamases are enzymes produced by bacteria that deactivate beta-lactam antibiotics such as ceftazidime.

This study presents a computational framework to predict whether a bacterial strain is:

- Resistant to ceftazidime
- Sensitive to ceftazidime

using only the beta-lactamase protein sequence.

The developed system uses machine learning techniques and sequence-derived protein descriptors for classification.

---

## Dataset Information

The dataset was collected from the β-lactamase database (BLDB).

### Dataset Statistics

- 199 beta-lactamase protein sequences
- 87 ceftazidime-sensitive proteins
- 112 ceftazidime-resistant proteins

### External Validation

- 22 resistant beta-lactamase sequences from the RGI database

---

## Feature Generation

Protein descriptors were generated using the **Pfeature** package.

A total of **9149 composition-based features** were calculated, including:

- Amino Acid Composition (AAC)
- Dipeptide Composition (DPC)
- Tripeptide Composition (TPC)
- Shannon Entropy
- Pseudo Amino Acid Composition (PAAC)
- Quasi Sequence Order (QSO)
- Distance Distribution of Residues (DDOR)

---

## Feature Selection

The following feature selection methods were used:

- Variance Threshold Method
- SVC-L1 Regularization
- Feature Ranking

Final models were developed using:

- Top 10 features
- Top 20 features
- Top 33 selected features

---

## Machine Learning Algorithms

The following classifiers were implemented:

- Random Forest (RF)
- Support Vector Classifier (SVC)
- Logistic Regression (LR)
- XGBoost (XGB)
- K-Nearest Neighbor (KNN)
- Decision Tree (DT)
- Gaussian Naive Bayes (GNB)

---

## Best Model Performance

### Random Forest Model (33 Features)

| Dataset | AUROC | MCC | Accuracy |
|----------|--------|------|-----------|
| Training | 0.80 | 0.48 | 74.10% |
| Validation | 0.79 | 0.46 | 73.33% |

### External Validation

- Correctly predicted 19 out of 22 resistant sequences

---

## Biological Insights

The analysis revealed that:

### Resistant Beta-Lactamases

Enriched in amino acids with non-polar side chains:

- Alanine (A)
- Glycine (G)
- Leucine (L)
- Proline (P)
- Arginine (R)

### Sensitive Beta-Lactamases

Enriched in polar and charged residues:

- Aspartic Acid (D)
- Isoleucine (I)
- Lysine (K)
- Asparagine (N)
- Threonine (T)
- Tyrosine (Y)

---

## Webserver Features

### Predict Module

Users can:

- Submit protein sequences in FASTA format
- Upload multiple sequences
- Predict resistance/susceptibility
- Obtain Random Forest prediction scores

### Standalone Package

A standalone Linux/Unix package is also available for offline prediction.

Standalone:
https://webs.iiitd.edu.in/raghava/abcrpred/stand.php

---

## Applications

- Antibiotic resistance prediction
- AMR surveillance
- Beta-lactamase analysis
- Personalized antibiotic therapy
- Computational microbiology research

---

## Authors
Gajendra Pal Singh Raghava
raghava@iiitd.ac.in
Department of Computational Biology  
IIIT Delhi, India

---

## License

This work is intended for academic and research purposes.

---

## Reference

bioRxiv:
https://doi.org/10.1101/2021.06.26.450028

Webserver:
https://webs.iiitd.edu.in/raghava/abcrpred/
