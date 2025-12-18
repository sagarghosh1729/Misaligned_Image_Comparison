# A Landmark-based Framework for Image Comparison in the Presence of Various Geometric Misalignments

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
- Location of Maximal Change Detection

The code is modularized:

- Helper functions in `R/`
- Main analysis in `Misaligned_Image_Comparison.Rmd`

---

## Dependencies

Install the required R packages:

```r
install.packages(c(
  "raster", "sp", "magick", "OpenImageR", "foreach",
  "doParallel", "twosamples", "magic", "pracma", "jpeg"
))


