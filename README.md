# PCR-bias

This repository contains the code and data associated with the study:

> **Deep learning uncovers sequence-specific amplification bias in multi-template PCR**

It includes scripts for deep learning model validation and motif discovery, aiming to analysis sequence-specific biases that arise during PCR amplification.

---

## Contents

- **CluMo.py**  
  A Python script that implements the **motif discovery** approach introduced in the study.  

- **InternalValidation.py**  
  Performs **5-fold internal validation** on the selected dataset.  

- **ExternalValidation.py**  
  Performs **external validation** on the selected dataset and evaluate all other datasets to measure generalization.  

- **analysis/**  
  Contains additional scripts/notebooks for analyzing results and generating figures for the manuscript.  

- **utils/**  
  Utility functions for data loading, preprocessing, model construction, and training.

- **Data/**  
  DNA sequence dataset with binarization of PCR efficiency under different thresholds.

---

## Dependencies
The software is implemented using Python 3.9.7.

All major dependencies can be found in `requirements.txt`.  

You can install these packages by running:

```bash
pip install -r requirements.txt
```
---

## Quick Demo

Running motif discovery
```bash
python CluMo.py --filename dataset
```
The results will be saved under `CNN/motifs/{dataset}/{threshold}/`

Running internal and external validation
```bash
python InternalValidation(ExternalValidation).py --filename dataset
```
The results will be saved under `CNN/results/interal(external)/{dataset}/{threshold}/`

The `dataset` should be specified as one of the 7 datasets used in this study: 
`
- "Choi_et_al",
- "Erlich_et_al",
- "Gao_et_al",
- "GCall",
- "GCfix",
- "Koch_et_al",
- "Song_et_al"
`

---

## Citation

Please use the following to cite our work:

```bibtex
@article{gimpel2024deep,
  title={Deep learning uncovers sequence-specific amplification bias in multi-template PCR},
  author={Gimpel, Andreas L and Fan, Bowen and Chen, Dexiong and W{\"o}lfle, Laetitia OD and Horn, Max and Meng-Papaxanthos, Laetitia and Antkowiak, Philipp L and Stark, Wendelin J and Christen, Beat and Borgwardt, Karsten and others},
  journal={bioRxiv},
  pages={2024--09},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```
