# Mendelian-Randomization-and-Multi-Omics-Analysis

## single-cell RNA sequencing analysis 
* Cellsubtype discovery.R assists in identifying cellular targets through single-cell RNA sequencing analysis. The process includes the following steps.
  ### Data reading
  ### Quality control
  ### Normalize
  ### Dimensionality reduction and clustering
  ### Cell type annotation
* The focus was on CD8+ T cells, and a similar analysis process was conducted to further identify CD8+ Teff cells. The key difference from the previous step lies in the cell type annotation method: here, manual annotation was employed instead of using singleR.
  ### T cells trajectory analysis
  * T cell trajectory analysis.R helps reveal the developmental relationships among different T cell subtypes, aiding in comprehensive subtype recognition.
  ### Cell chat analysis
  * cellchat.R facilitates the exploration of intercellular communication between CD8+ Teff cells and other cell types.

## MR analysis
* MR analysis.R aids in identifying SNPs within CD8+ Teff cells that may contribute to the pathogenesis of SSc-ILD.
* MR analysis with validation test.R modifies the outcome dataset to validate SNPs within CD8+ Teff cells implicated in the pathogenesis of SSc-ILD.
* biMR analysis.R further confirmed that SNPs cause SSc-ILD rather than the reverse. The SNPs within the N4BP2L1 gene have been confirmed to contribute to the pathogenesis of SSc-ILD.
* regional association plot.R helps reveal the genetic association analysis of SNPs within N4BP2L1 with traits related to SSc-ILD.
* steiger filtering.R computes the amount of variance each SNP explains in the exposure and in the outcome variable.

## Downstream analysis of the gene N4BP2L1 regarding its impact on the pathogenesis of SSc-ILD and IPF


### 三级标题

* 要点一
* 要点二
