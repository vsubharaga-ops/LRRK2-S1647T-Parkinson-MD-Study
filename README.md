# LRRK2 S1647T Mutation — Parkinson's Disease Study

## Integrative WES and Molecular Dynamics Analysis

**Institution:** Bversity School of Life Sciences  
**Company:** Insilicomics, Ooty  
**Authors:** Nandini Solanki | V Subharaga  
**Date:** April 2026  

---

## Project Overview

This project identifies and characterises the **S1647T missense mutation** in the 
**LRRK2 (PARK8)** gene associated with Parkinson's disease using an integrative 
computational pipeline combining:

- Whole-Exome Sequencing (WES) variant discovery
- All-atom Molecular Dynamics (MD) simulations

---

## Repository Structure
---

## WES Pipeline Summary

| Step | Tool | Purpose |
|------|------|---------|
| QC | FastQC | Raw read quality assessment |
| Trim | Trimmomatic | Adapter and quality trimming |
| Align | BWA-MEM → GRCh38 | Read alignment |
| Post-process | SAMtools + Picard | BAM processing |
| Variant call | GATK HaplotypeCaller | SNP/INDEL calling |
| Filter | GATK Hard filter | Quality filtering |
| Annotate | ANNOVAR | ClinVar/dbSNP annotation |
| Prioritise | LRRK2 filter | S1647T selected |

**Dataset:** SRR19500742 (NCBI SRA)  
**Result:** S1647T identified from 3 LRRK2 missense variants in 87,000 raw calls

---

## MD Simulation Summary

| Parameter | Value |
|-----------|-------|
| Software | GROMACS 2023.3 |
| Force field | AMBER99SB-ILDN |
| Water model | TIP3P |
| Box type | Dodecahedron, 1.0 nm padding |
| Temperature | 300 K (V-rescale) |
| Pressure | 1 bar (Parrinello-Rahman) |
| Production MD | 100 ps (50,000 steps, dt=2 fs) |
| System size | ~355,807 atoms |

---

## Key Results

| Parameter | WT | S1647T | Δ |
|-----------|-----|--------|---|
| RMSD (nm) | 0.218 ± 0.058 | 0.181 ± 0.062 | −0.037 |
| Rg (nm) | 4.505 ± 0.012 | 4.525 ± 0.028 | +0.020 |
| H-bonds | 1206 ± 15.3 | 1211 ± 13.9 | +5.3 |
| SASA (nm²) | 882.34 ± 13.8 | 866.98 ± 17.7 | −15.36 |

**Conclusion:** S1647T introduces subtle structural stabilisation of the LRRK2 
COR domain with allosteric propagation to distal domains.

---

## Requirements

```bash
# GROMACS
sudo apt install gromacs

# Python packages
pip install matplotlib numpy biopython pdbfixer
```

---

## How to Reproduce

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/LRRK2-S1647T-Parkinson-MD-Study

# 2. Introduce S1647T mutation
python3 scripts/mutate_S1647T.py

# 3. Run MD (follow MDP files in mutant_md/ and wt_md/)
gmx pdb2gmx -f data/LRRK2_WT_SER1647_clean.pdb -ff amber99sb-ildn -water tip3p

# 4. Generate comparison plots
cd analysis/
python3 ../scripts/final_plot.py
```

---

## Citation

If you use this work, please cite:
> Solanki N, Subharaga V. (2026). Integrative WES and MD study of S1647T 
> mutation in LRRK2 for Parkinson's disease. Insilicomics Mini-Project CP003.

---

## License

MIT License — see LICENSE file for details.
