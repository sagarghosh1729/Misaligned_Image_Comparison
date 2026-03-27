# A Landmark-based Framework for Image Comparison in the Presence of Various Geometric Misalignments

This repository contains the full R implementation for the study:

**A Landmark-based Framework for Image Comparison in the Presence of Various Geometric Misalignments**  
**Authors:** Anik Roy, Sagar Ghosh, Partha Sarathi Mukherjee

---

## Overview

This project provides tools and scripts for **A Landmark-based Framework for robust image comparison** under various geometric misalignments. It includes:

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

## Required Libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EBImage")


library(EBImage)
library(raster)
library(sp)
library(magick)
library(OpenImageR)
library(foreach)
library(doParallel)
library(DRIP)
library(twosamples)
library(OpenImageR)
library(magic)
library(pracma)
library(jpeg)
library(SAFARI)
#library(imager)



## Operating System and R Studio Environment
Hardware Overview:

      Model Name: MacBook Pro
      Model Identifier: Mac16,6
      Chip: Apple M4 Max
      Total Number of Cores: 14 (10 Performance and 4 Efficiency)
      Memory: 36 GB
      System Firmware Version: 18000.101.7
      OS Loader Version: 18000.101.7

R Studio Overview:

      RStudio 2026.01.1+403 "Apple Blossom" Release (0e924abb984501b0d66b204ea06b60fc7813275a, 2026-02-04) for         macOS




