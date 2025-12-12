# Image Comparison under Geometric Misalignments

This repository contains the full R implementation for the study:

**“Image Comparison in the Presence of Various Geometric Misalignments”**  
**Authors:** Anik Roy, Sagar Ghosh, Partha Sarathi Mukherjee

---

## Overview

This project provides tools and scripts for **robust image comparison** under various geometric misalignments. It includes:

- Simulation experiments
- Landmark-based shape comparisons
- TRS-invariant metrics
- JLC-based edge extraction
- Real tumor growth analysis

The code is modularized:

- Helper functions in `R/`
- Main analysis in `code/Misaligned_Image_Comparison.Rmd`
- Data in `data/`
- Manuscript drafts in `paper/`
- Figures and outputs can be saved locally (e.g., `results/`)

---

## Repository Structure

.Rhistory
.RData
.Rproj.user/
.DS_Store
results/*
data/*.jpg
data/*.jpeg



---

## Dependencies

Install the required R packages:

```r
install.packages(c(
  "raster", "sp", "magick", "OpenImageR", "foreach",
  "doParallel", "twosamples", "magic", "pracma", "jpeg"
))


