# How to implement the MiCA algorithm
## Vishal Midya, Chris Gennings, and Shoshannah Eggers

This article presents a step-by-step guide and intuitions to implement the MiCA algorithm. MiCA is conducted in two stages to identify microbial cliques associated with the outcome of interest. The first part of this algorithm uses a machine learning-based prediction framework to discover microbial cliques predictive of the outcome. The next stage restores the directionality and is divided into estimating the association between the outcome and the joint-relative abundance of the discovered cliques using a causal inference (or classical association-based) framework. The microbial cliques are searched using repeated hold-out signed-iterated Random Forest (rh-SiRF), where the predictors are relative abundances of the selected taxa. The SiRF (Signed Iterative Random Forest) algorithm combined a state-of-the-art predictive tool called "Iterative Random Forests" with "Random Intersection Trees" to search for microbial cliques predictive of the outcome <sup>1</sup>. On top of the SIRF algorithm, we introduced a repeated hold-out step that randomly partitions the data in training and testing sets for better generalizability. See more technical details in <sup>2</sup>.

##  Simulated Taxa data

We first present simulated data on around `500` hypothetical participants with `60` simulated Taxa. The dataset, named `data.simulated.csv`, is uploaded as a part of the demonstration.

## Required `R` packages

Kindly install the following `R` packages
1. `mvtnorm`
2. `iRF`

Find the instructions [here](https://github.com/sumbose/iRF/tree/master) to install the `iRF` package.


## Creating a 3<sup>rd</sup> order microbial clique and an outcome

1. __Run the following chuck of functions:__

`require(mvtnorm)`

```{}
make.SIGMA <- function(rho,dimension){
  SIGMA = matrix(NA,dimension,dimension)
  for(i in 1:dimension){
    for(j in 1:dimension){
      if(i != j){
        a <- sign(rnorm(1,0,1))*rnorm(1,0.1, 0.01)
        if(a>0) {SIGMA[i,j] = rho + a}
        else if (a <0)  {SIGMA[i,j] = rho}
      }
      else if (i == j){
        SIGMA[i,j] = 1
      }
    }
  }
  
  (SIGMA + t(SIGMA))/2
}
```

```{}
make.X0 <- function(n, rho, p0){
  X0 <- (mvtnorm::rmvnorm(n, mean = rep(0, nrow(make.SIGMA(rho, p0))), sigma = make.SIGMA(rho, p0)))
}
```

2. __Next, we create four covariates that are moderately correlated to each other__

```{}
n <- dim(sample.data.simulated)[1]
set.seed(12456)
covariates <- make.X0(n,0.1,4)
```
The set of four covariates is already supplied in the dataset.  

3. __Create a 3<sup>rd</sup> order microbial clique__

Without loss of generality, we choose `Taxa.1`, `Taxa.3`, and `Taxa.11` to form a clique. The clique is present only in those individuals where (1) the relative abundance of `Taxa.1` is _less_ than its total sample median, (2) the relative abundance of `Taxa.3` is _greater_ than its total sample median, and lastly, (3) the relative abundance of `Taxa.11` is _greater_ than its total sample median. In terms of code,

```{}
three_clique <- as.numeric(data.simulated$Taxa.1 <= quantile(data.simulated$Taxa.1, 0.5)
                             & data.simulated$Taxa.3 >= quantile(data.simulated$Taxa.3, 0.5)
                                & data.simulated$Taxa.11 >= quantile(data.simulated$Taxa.11, 0.5))
```

4. __Create the outcome__

We create a Gaussian outcome composed of the indicator for microbial clique, covariates, and random Gaussian error.

```{}
set.seed(45667)
outcome <-  1 * three_clique  + covariates %*% rep(0.1, 4) + rnorm(n, 0, 0.5)
```
This outcome is already supplied in the dataset. However, note that the effect size of the `three_clique` is set at `1`.




## References

1. Basu, S.; Kumbier, K.; Brown, J. B.; Yu, B. Iterative Random Forests to Discover Predictive and Stable High-Order Interactions. Proc. Natl. Acad. Sci. 2018, 115 (8), 1943–1948.
2. Midya V, Lane JM, Gennings C, Torres-Olascoaga LA, Wright RO, Arora M, Tellez-Rojo MM, Eggers S. Prenatal Pb exposure is associated with reduced abundance of beneficial gut microbial cliques in late childhood: an investigation using Microbial Co-occurrence Analysis (MiCA). medRxiv. 2023 May. doi: https://doi.org/10.1101/2023.05.18.23290127.

