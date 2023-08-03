# How to implement the MiCA algorithm
## Vishal Midya, Chris Gennings, and Shoshannah Eggers

This article presents a step-by-step guide and intuitions to implement the MiCA algorithm.  MiCA is conducted in two stages to identify microbial cliques associated with the outcome of interest. The first part of this algorithm uses a machine learning-based prediction framework to discover microbial cliques predictive of the outcome. The next stage restores the directionality and is divided into estimating the association between the outcome and the joint-relative abundance of the discovered cliques using a causal inference (or classical association-based) framework. The microbial cliques are searched using repeated hold-out signed-iterated Random Forest (rh-SiRF), where the predictors are relative abundances of the selected taxa. The SiRF (Signed Iterative Random Forest) algorithm combined a state-of-the-art predictive tool called "Iterative Random Forests" with Random Intersection Trees (RIT) to search for combinations of taxa predictive of the outcome <sup>1</sup>. On top of the SIRF algorithm, we introduced a repeated hold-out step that randomly partitions the data in training and testing sets for better generalizability. 




`library(kableExtra)`

## References

1. Basu, S.; Kumbier, K.; Brown, J. B.; Yu, B. Iterative Random Forests to Discover Predictive and Stable High-Order Interactions. Proc. Natl. Acad. Sci. 2018, 115 (8), 1943–1948.
2. 
