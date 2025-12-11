# Image Comparison in the Presence of various Geometric Misalignments  
Code for the paper:  
**“Image Comparison in the Presence of Various Geometric Misalignment”**  
Authors: Anik Roy, Sagar Ghosh, Partha Sarathi Mukherjee  

---

##  Overview  
This repository contains the full R implementation used in our study on robust 
image comparison under geometric misalignment. It includes simulation experiments, 
landmark-based shape comparisons, TRS-invariant metrics, JLC-based edge extraction, 
and real tumor growth analysis.

The code is modularized: helper functions reside in `R/`, the full pipeline is in 
`code/misalignment-analysis.Rmd`, and output figures are stored in `results/`.

---

## Repository Structure

|--code/
|   |-Misaligned_Image_Comparison.Rmd # Main R Markdown
|--R/
|   |-Centering_Matrix.R
|   |-Frobenius_Norm.R
|   |-JLC_Sample_Edge.R
|   |-Landmark_Selection.R
|   |-Plot_Landmark_Metrix.R
|   |-TRS_Inv_Metric.R
|--data/
|   |-Ellipse
|   |-Lake
|   |-Tumor
|--paper/
|--LICENSE
|--README.md


--

##  Dependencies  

Install R packages:

```r
install.packages(c(
  "raster", "sp", "magick", "OpenImageR", "foreach",
  "doParallel", "twosamples", "magic", "pracma", "jpeg"
))
