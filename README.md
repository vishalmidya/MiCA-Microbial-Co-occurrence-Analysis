# Microbial-Co-occurrence-Analysis (MiCA)
## Vishal Midya, Chris Gennings, and Shoshannah Eggers

This repository contains a vignette of how to apply the Microbial-Co-occurrence-Analysis (MiCA) analysis to discover gut microbial cliques associated with any health outcome. MiCA is conducted in two stages to identify microbial cliques associated with chemical exposures. The first part of this algorithm uses a machine learning-based prediction framework to discover microbial cliques predictive of the health outcome. The next stage restores the directionality and dives into estimating the association between the joint-relative abundance of the discovered cliques and the health outcome using a causal inference (or simply classical association-based) framework. We also include a CSV file containing simulated relative abundance data for 60 Taxa for about 500 people. It also contains four covariates and a continuous outcome. This fictitious dataset demonstrates why and how to use the Microbial-Co-occurrence-Analysis (MiCA) analysis. Please read the `MiCA-vignette.md` tab for further discussion. 


## Reference

Midya V, Lane JM, Gennings C, Torres-Olascoaga LA, Wright RO, Arora M, Tellez-Rojo MM, Eggers S. Prenatal Pb exposure is associated with reduced abundance of beneficial gut microbial cliques in late childhood: an investigation using Microbial Co-occurrence Analysis (MiCA). medRxiv. 2023 May. doi: https://doi.org/10.1101/2023.05.18.23290127.
