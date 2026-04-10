# A Landmark-based Framework for Image Comparison in the Presence of Various Geometric Misalignments

This repository contains the full R implementation for the study:

**A Landmark-based Framework for Image Comparison in the Presence of Various Geometric Misalignments**  
**Authors:** Anik Roy, Sagar Ghosh, Partha Sarathi Mukherjee

---

## Overview

This project provides tools and scripts for **A Landmark-based Framework for robust image comparison** under various geometric misalignments. It includes:

- Simulation experiments : The Simulations folder contains three simulated experiments on power computations: One use cases each from Amoeba, Ellipse and Polygon. This folder also contains one use case of the local change detection simulation on two lake images. 
  
- R: This Folder contains all the functions and an end-to-end function named as "Power_Computation.R" which intakes two image paths and produces the power output considering the first image as null. This folder also contains the file "local_change_detection_fucntion.R" which has an end-to-end function taking two image paths as input and producing the location on the second image where the maximal change is observed relative to the first image
  
- Datasets: This folder contains all the **lake images** and a link to the **medical images** which were used to detect the position of the maximal local change. 
  
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
## Required Libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))  
#  install.packages("BiocManager")

#BiocManager::install("EBImage")


library(EBImage)
#library(raster)
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

      RStudio 2026.01.1+403 "Apple Blossom" Release (0e924abb984501b0d66b204ea06b60fc7813275a, 2026-02-04) for  macOS




