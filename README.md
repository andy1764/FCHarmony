# FCHarmony
### Methods for harmonization of functional connectivity matrices

--------
**Maintainer**: Andrew Chen, andrewac@pennmedicine.upenn.edu

**License**: Artistic License 2.0

**References**: If you are using FC-ComBat, please cite the following paper

|       | Citation     | Paper Link
| -------------  | -------------  | -------------  |
| ComBat for functional connectivity    | Yu, Meichen, Kristin A. Linn, Philip A. Cook, Mary L. Phillips, Melvin McInnis, Maurizio Fava, Madhukar H. Trivedi, Myrna M. Weissman, Russell T. Shinohara, and Yvette I. Sheline. 2018. “Statistical Harmonization Corrects Site Effects in Functional Connectivity Measurements from Multi-Site FMRI Data.” Human Brain Mapping 39 (11): 4213–27. https://doi.org/10.1002/hbm.24241. | [Link](https://doi.org/10.1002/hbm.24241)|

## Table of content
- [1. Installation](#id-section1)
- [2. Background](#id-section2)
- [3. Software](#id-section3)

<div id='id-section1'/>

## 1. Installation
The R package and its Github dependencies can be installed via devtools by 
running the following code

```
# install.packages("devtools")
devtools::install_github("andy1764/CovBat_Harmonization/R")
devtools::install_github("andy1764/FCHarmony")
```

Then, you can load this package via

```
library(FCHarmony)
```

The R package provides the `fcComBat`, `blComBat`, and `fcCovBat` functions for harmonization of functional connectivity. We include `test_regress` for regression-based site effect evaluation and `plot_fc` for plotting results from `test_regress` or plotting functional connectivity matrices.

<div id='id-section2'/>

## 2. Background
Community detection on graphs constructed from functional magnetic resonance imaging (fMRI) data has led to important insights into brain functional organization. Large studies of brain community structure often include images acquired on multiple scanners across different studies. Differences in scanner can introduce variability into the downstream results, and these differences are often referred to as scanner effects. Such effects have been previously shown to significantly impact common network metrics. 

In our preprint, we identify scanner effects in data-driven community detection results and related network metrics. We assess the performance of [FC-ComBat](https://doi.org/10.1002/hbm.24241) and introduce two new methods for harmonizing function connectivity:

1. **Block-ComBat (Bl-ComBat)**, which leverages existing knowledge about network structure
2. **FC-CovBat**, which harmonizes patterns of covariance in the data using [CovBat](https://github.com/andy1764/CovBat_Harmonization).

We demonstrate that our new methods reduce scanner effects in community structure and network metrics. Our results highlight scanner effects in studies of brain functional organization and provide additional tools to address these unwanted effects. `FCHarmony` implements our harmonization methods and provides tools to evaluate harmonization performance. 

<div id='id-section3'/>

## 3. Software
The R and Python implementations of CovBat are available [here](https://github.com/andy1764/CovBat_Harmonization) and are extensions of ComBat implemented by Jean-Phillipe Fortin in R, Python, and Matlab in the  [ComBat](https://github.com/Jfortin1/ComBatHarmonization) package maintained by Jean-Philippe Fortin. Network analyses in our paper largely utilize the [Brain Connectivity Toolbox version Version 2019-03-03](https://sites.google.com/site/bctnet/).
