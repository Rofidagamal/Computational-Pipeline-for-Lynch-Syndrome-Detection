# Lynch Syndrome Variant Detection Pipeline — Streamlit GUI

A no-code web interface for detecting Lynch Syndrome-associated germline variants
from RNA-seq or VCF data. Built for clinical biologists without programming experience.

## Quick Start

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then open http://localhost:8501 in your browser.

## Features

- **Upload & Run** — Upload FASTQ / VCF / BAM and run the full 10-step pipeline
- **Variant Results** — Filterable variant table (gene, ClinVar, consequence) with CSV export
- **Statistics** — Sensitivity/specificity gauges, confusion matrix, FASTQ QC metrics, Plotly charts
- **Pipeline Steps** — Expandable documentation for all 10 steps with parameters and tips
- **Help & Docs** — Installation guide, file format specs, FAQ

## Accepted File Formats

| Format | Extension |
|--------|-----------|
| Raw reads | `.fastq`, `.fq`, `.fastq.gz` |
| Variants | `.vcf`, `.vcf.gz` |
| Aligned reads | `.bam` |

## MMR Genes Analyzed (GRCh38)

| Gene  | Chromosome | Start     | End       |
|-------|-----------|-----------|-----------|
| MLH1  | chr3      | 36,993,336| 37,050,918|
| MSH2  | chr2      | 47,403,067| 47,709,873|
| MSH6  | chr2      | 47,783,192| 47,806,321|
| PMS2  | chr7      | 5,970,838 | 6,009,982 |
| EPCAM | chr2      | 47,291,953| 47,344,088|

## QC Thresholds (Step 9)

- QUAL > 200
- DP ≥ 10×
- GQ ≥ 20

## Citation

If you use this pipeline, please cite:
> Gamal R. et al. "A Reproducible Computational Pipeline for Detecting Lynch Syndrome
> Variants from RNA Sequencing Data." IEEE NILES 2025.

GitHub: https://github.com/Rofidagamal/Computational-Pipeline-for-Lynch-Syndrome-Detection
