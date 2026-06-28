# Adjusted Positive Agreement (APA): A Chance-Corrected Agreement Measure for Unconstrained Pile Sorts

## Publication Information

**Title:** Adjusted Positive Agreement: A Chance-Corrected Agreement Measure for Unconstrained Pile Sorts

**Author:** Jonas Idevall Hagren

**Journal:** Field Methods

**Publication Year:** 2026

**DOI:** 10.1177/1525822X261464653

## Overview

This directory contains the simulation study code and analysis scripts for the peer-reviewed publication on Adjusted Positive Agreement (APA). The code demonstrates a comprehensive comparison of three respondent-by-respondent agreement measures for unconstrained pile sorts:

- **APA**: Adjusted Positive Agreement (novel measure)
- **Pearson Correlation**: On binary co-placement vectors
- **ARI**: Adjusted Rand Index

## Simulation Study

The simulation evaluates agreement measures under two experimental conditions:

### Condition 1: Single Culture Scenario
All respondents sort according to the same taxonomic hierarchy (H), but differ in classificatory resolution:
- Low resolution (5-8 piles)
- Intermediate resolution (9-11 piles) 
- High resolution (12-16 piles)

### Condition 2: Multiple Cultures Scenario
- 75% of respondents sort according to taxonomic hierarchy (H)
- 25% sort according to an ecological organization (H*) that cross-cuts the taxonomy

## Main Code

**File:** `APA_FM_2026_v3.R`

The R script includes:

1. **Setup and data generation:** Animal taxonomy and trait-based ecological organization
2. **Latent item organizations:** Defines two distinct categorization systems (H and H*)
3. **Visualizations:** Generates figures for hierarchical and cross-cutting organizations
4. **Respondent partition generation:** Simulates respondents with varying resolution levels
5. **Agreement calculations:** Implements APA, Pearson, and ARI measures
6. **Simulation framework:** Runs 1000 simulation replicates across both conditions
7. **Eigenanalysis:** Performs spectral decomposition and eta-squared recovery analysis
8. **Output generation:** Produces figures and data tables for publication

## Requirements

### R Libraries
```r
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(igraph)
library(ggraph)
library(grid)
library(patchwork)
```

### System Requirements
- R version 3.6 or higher
- Approximately 2GB of RAM (for 1000 simulations)
- Runtime: ~2-3 hours on typical desktop hardware

## Running the Analysis

1. Ensure all required libraries are installed
2. Set `SAVE_FIGURES` and `SAVE_TABLES` flags as needed
3. Run the script from start to finish:
   ```r
   source("APA_FM_2026_v3.R")
   ```

## Output

### Figures Generated
- `Hstar_size_habitat_grid.png` - Ecological organization visualization
- `H_taxonomic_hierarchy_circular.png` - Taxonomic hierarchy diagram
- `Figure_2_latent_organisations_H_and_Hstar.png` - Combined organization figure
- `sim_heatmap_single_culture_APA_Pearson_ARI.png` - Single-culture agreement comparison
- `sim_eigenvector_separation_single_culture_APA_Pearson_ARI.png` - Eigenvector recovery (single culture)
- `Figure_single_culture_heatmap_and_resolution_separation.png` - Combined single-culture figure
- `sim_eigenvalue_spectra_multiple_cultures_APA_Pearson_ARI.png` - Eigenvalue spectra
- `sim_eigenvector_separation_multiple_cultures_APA_Pearson_ARI.png` - Eigenvector recovery (multiple cultures)
- `Figure_eigenvalue_spectra_and_eigenvector_separation.png` - Combined eigenanalysis figure

### Data Tables Generated
- `sim_summary.csv` - Overall summary statistics by simulation
- `sim_pairwise.csv` - Pairwise respondent agreement values
- `sim_metadata.csv` - Respondent metadata and characteristics
- `sim_structure_pair_summary.csv` - Agreement by latent organization pairing
- `sim_eigenvalues_APA_Pearson_ARI.csv` - Eigenvalue decomposition results
- `sim_eigenvectors_APA_Pearson_ARI.csv` - Eigenvector coordinates for respondents
- `sim_eigenvector_separation_by_sim.csv` - Eta-squared recovery by simulation
- `sim_eigenvector_separation_summary.csv` - Summary of eigenvector separation
- `sim_eigenvalue_ratio_multiple_summary.csv` - Eigenvalue ratio analysis for multiple-culture condition

## Key Findings

The simulation demonstrates that APA provides a chance-corrected measure of respondent agreement that:
- Preserves agreement across respondents sharing the same organizational framework but differing in resolution
- Distinguishes between different latent categorization systems
- Shows interpretable eigenvalue spectra that reflect the underlying organizational structure

## Citation

If you use this code or findings in your research, please cite:

Idevall Hagren, J. (2026). Adjusted Positive Agreement: A Chance-Corrected Agreement Measure for Unconstrained Pile Sorts. *Field Methods*, [Volume(Issue)], pages. https://doi.org/10.1177/1525822X261464653

## License

This code is provided for research and educational purposes. Please refer to the journal publication for methodological details and theoretical background.

## Questions or Issues

For questions about the implementation or methodology, please refer to the published article or contact the author through the journal.
