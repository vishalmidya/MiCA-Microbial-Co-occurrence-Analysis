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

## The aim of the MiCA algorithm

1. To recover the 3<sup>rd</sup> order microbial clique of `Taxa.1`, `Taxa.3`, and `Taxa.11`.
2. To recover the thresholds for the construction of the clique.
3. To estimate the effect of the clique on the outcome.

Note that, there are a total of `1770` two-ordered, and `34220` three-ordered combinations to choose from `60` Taxa. Further, as the number of Taxa increases, the total number of possible combinations to test from exponentially increases. Therefore, we used the signed-iterated Random Forest (SiRF) algorithm to search for the optimal combinations. Further, we introduced a repeated-holdout stage on the SiRF algorithm to keep the false positives and negatives in check.

### Finding the optimal combination of Taxa


Run the following function that finds the most frequently occurring combinations of Taxa:

```{}
require("iRF")
clique.finder <- function(exposures, outcome, iterations, validation, seed.value, n.bootstrap, min.prevalence, min.stability, data){
  n <- dim(data)[1]
  if(n == 0){
    stop("Provide a dataset in data.frame format") 
  }
  if(length(exposures) < 2){
    stop("Need at least two exposures to find a clique") 
  }
  if(validation == 0 | validation >=1){
    stop("Validation must be a positive fraction within 0 and 1") 
  }
  train.sample <- 1 - validation
  if(iterations < 50){
    warning("Use more than 50 iterations for reliable replication")
  }
  if(n.bootstrap < 100){
    warning("Use more than 100 bootstraps for reliable replication")
  }
  tab <- data.frame( int = NA_character_,  prevalence = NA_real_, precision = NA_real_, cpe = NA_real_,
                     sta.cpe = NA_real_,  fsd = NA_real_, sta.fsd = NA_real_, mip = NA_real_,
                     sta.mip = NA_real_, stability =NA_real_, occurance.freq = NA_character_ , occurance.freq.sum = NA_real_)
  for(i in 1:iterations){
    set.seed(runif(1, 0, (seed.value + 10)))
    train.id <- sample(seq(1,n), ceiling(n*train.sample))
    test.id <- setdiff(1:n, train.id)
    fit <- iRF(x=data.simulated[train.id, exposures], 
               y=data.simulated[train.id,outcome], 
               xtest=data.simulated[test.id, exposures],
               ytest=data.simulated[test.id,outcome],
               n.iter=10, 
               n.core=3,
               select.iter = T,
               n.bootstrap=n.bootstrap
    )
    SiRF_table <- as.data.frame(fit$interaction[fit$interaction$stability > min.stability,])
    if(dim(SiRF_table)[1] != 0){
      
      SiRF_table_subset <- SiRF_table
      SiRF_table_subset_prev <- SiRF_table_subset[SiRF_table_subset$prevalence > min.prevalence,]
      
      if(dim(SiRF_table_subset_prev)[1] != 0){
        nam <- NA_character_
        for(i in 1:nrow(SiRF_table_subset_prev)){nam <- c(nam, strsplit(SiRF_table_subset_prev$int[i],"_")[[1]])}
        yt<-as.data.frame(table(nam))
        SiRF_table_subset_prev$occurance.freq <- rep(NA_character_,nrow(SiRF_table_subset_prev))
        for(i in 1:nrow(SiRF_table_subset_prev)){
          SiRF_table_subset_prev$occurance.freq[i] <- as.character((paste0(yt$Freq[which(yt$nam == strsplit(SiRF_table_subset_prev$int[i],"_")[[1]][1])],"/",yt$Freq[which(yt$nam == strsplit(SiRF_table_subset_prev$int[i],"_")[[1]][2])])))
        }
        SiRF_table_subset_prev$occurance.freq.sum <- rep(NA_character_,nrow(SiRF_table_subset_prev))
        for(i in 1:nrow(SiRF_table_subset_prev)){
          SiRF_table_subset_prev$occurance.freq.sum[i] <- sum(strsplit(SiRF_table_subset_prev$occurance.freq[i],"/")[[1]] %in% "1")
        }
        tab <- rbind(tab, SiRF_table_subset_prev)
      } else {
        tab <- rbind(tab, SiRF_table_subset_prev)
      }
    } else {
      tab <- rbind(tab, SiRF_table)
    }
  }
  clique <- tab
  if(dim(tab)[1] == 1){
    stop("No clique was found, decrease the value of min.stability and increase the number of bootstraps")
  }
  clique <- clique[-1,]
  for(i in 1:nrow(clique)){
    if(sum(strsplit(clique$int[i],"")[[1]] %in% "+") != 0 & sum(strsplit(clique$int[i],"")[[1]] %in% "-") == 0){
      x <- strsplit(clique$int[i],"+")[[1]]
      clique$int[i] <- paste0(x[x!= "+"], collapse = "")  
    } else if(sum(strsplit(clique$int[i],"")[[1]] %in% "-") != 0 & sum(strsplit(clique$int[i],"")[[1]] %in% "+") == 0){
      x <- strsplit(clique$int[i],"-")[[1]]
      clique$int[i] <- paste0(x[x!= "-"], collapse = "")  
    } else if (sum(strsplit(clique$int[i],"")[[1]] %in% "-") != 0 & sum(strsplit(clique$int[i],"")[[1]] %in% "+") != 0){
      x <- strsplit(clique$int[i],"")[[1]]
      clique$int[i] <- paste0(x[x!= "-" & x!= "+"], collapse = "")  
    }
  }
  rtf <- as.data.frame(table(clique$int))
  rtf <- rtf[order(rtf$Freq, decreasing = T),]
  return(as.data.frame(rtf))
}
```

Finally, run the `function` called `clique.finder`. Below we discuss each argument for this function and what they entail.

```{}
clique.finder(exposures = paste0("Taxa.", seq(1,60)), outcome = "outcome",  iterations = 500, validation = 0.4, 
              seed.value = 1234, n.bootstrap = 200, min.prevalence = 0.05, min.stability = 0.25, data = data.simulated)
```

1. `exposures`: a vector of all possible Taxa names (among which one intends to find the combinations)
2. `outcome`: name of the outcome variable
3. `iterations`: the number of repeated holdouts (should be more than 100)
4. `validation`: the proportion of the dataset which is set aside for validation at each repeated holdout iteration 
5. `seed.value`: random initial seed value for partitioning the dataset
6. `n.bootstrap`: the number of bootstrap iterations employed at each of the repeated holdouts on the training dataset (should be more than 100)
7. `min.prevalence`: the minimum proportion (lower bound) of the sample that has the clique. Here we chose `5%` as the lower bound of the prevalence. 
8. `min.stability`: the stability implies the proportion of times the clique was recovered across bootstrap replicates. The `min.stability` is the lower bound. Here we chose `25%` as the lower bound.
9. `data`: name of the dataset

Note that, lowering the values of `min.prevalence` and `min.stability` finds more combinations of Taxa; however, due to the implementation of the repeated holdout technique, lowering these bounds do not have any significant effect on the most stable combinations. In this example, we used `1500` repeated holdouts. One can utilize the [parallel R package](https://www.rdocumentation.org/packages/parallel/versions/3.6.2) for fast parallel computation. 



### References

1. Basu, S.; Kumbier, K.; Brown, J. B.; Yu, B. Iterative Random Forests to Discover Predictive and Stable High-Order Interactions. Proc. Natl. Acad. Sci. 2018, 115 (8), 1943–1948.
2. Midya V, Lane JM, Gennings C, Torres-Olascoaga LA, Wright RO, Arora M, Tellez-Rojo MM, Eggers S. Prenatal Pb exposure is associated with reduced abundance of beneficial gut microbial cliques in late childhood: an investigation using Microbial Co-occurrence Analysis (MiCA). medRxiv. 2023 May. doi: https://doi.org/10.1101/2023.05.18.23290127.


### Acknowledgments

This method was developed at the Dept. of Environmental Medicine and Public Health, Icahn School of Medicine at Mount Sinai, NYC, with funding and support from the National Institute of Environmental Health Sciences (K99ES032884, P30ES023515, and U2C ES026555-01).
