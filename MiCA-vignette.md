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


1. __Run the following chunk of functions:__

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
This outcome is already supplied in the dataset. However, note that the effect size of the `three_clique` is set at `1`. Further, note that individually `(outcome ~ Taxa + cov1 + cov2 + cov3 + cov4)`, the relative abundances of the three Taxa are _not_ significantly associated with the outcome. Only through the three-ordered clique, is there a significant statistical association.  

```{}
             Estimate Std. Error t value Pr(>|t|)
Taxa.1      -0.009965   0.006567  -1.517 0.129816    
Taxa.3       0.007614   0.006652   1.145  0.25297    
Taxa.11      0.008184   0.006479   1.263 0.207165    
clique        1.24607    0.06482  19.223  < 2e-16 ***

```

## The aim of the MiCA algorithm

1. To recover the 3<sup>rd</sup> order microbial clique of `Taxa.1`, `Taxa.3`, and `Taxa.11`.
2. To recover the thresholds for the construction of the clique.
3. To estimate the effect of the clique on the outcome.

Note that, there are a total of `1770` two-ordered, and `34220` three-ordered combinations to choose from `60` Taxa. Further, as the number of Taxa increases, the total number of possible combinations to test from exponentially increases. Therefore, we used the signed-iterated Random Forest (SiRF) algorithm to search for the optimal combinations. Further, we introduced a repeated-holdout stage on the SiRF algorithm to keep the false positives and negatives in check.

## Finding the optimal combination of Taxa

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

Note that lowering the values of `min.prevalence` and `min.stability` finds more combinations of Taxa; however, due to the implementation of the repeated holdout technique, lowering these bounds does not significantly affect the most stable combinations. In this example, we used `1500` repeated holdouts. One can utilize the [parallel R package](https://www.rdocumentation.org/packages/parallel/versions/3.6.2) for fast parallel computation. The same code can be used for binary outcomes as well.

Here is the result of the top 10 combinations obtained from the simulated dataset

```{}
                       Var1    Freq
             Taxa.11_Taxa.3 11.9250
             Taxa.1_Taxa.11 11.5353
              Taxa.1_Taxa.3 10.2261
      Taxa.1_Taxa.11_Taxa.3  3.1224
            Taxa.11_Taxa.17  3.1154
            Taxa.11_Taxa.41  2.6534
             Taxa.1_Taxa.41  2.6510
             Taxa.1_Taxa.17  2.6440
            Taxa.11_Taxa.33  2.0443
             Taxa.17_Taxa.3  2.0303
```
The first column, `Var1`, denotes all possible combinations picked up by the algorithm, and the `Freq` denotes the percentage of their occurrence over all possible detected combinations and all bootstrap and repeated holdout combinations. The total number of detected combinations is around `400`. Note the top three combinations, `Taxa.11_Taxa.3`, `Taxa.1_Taxa.11`, and `Taxa.1_Taxa.3` occurred the most (more than 10%) among all possible detected combinations. Moreover, a three-ordered `Taxa.1_Taxa.11_Taxa.3` combination was among the most occurring. The three-ordered clique of `Taxa.1_Taxa.11_Taxa.3` induced down-stream of further two-ordered cliques. With smaller sample sizes, this tool effectively detects lower-ordered cliques. However, as the sample size increases, the frequency of the true higher-ordered combinations also increases. Although this is very subjective, a good practice is to choose the top (first three or first five) most frequently occurring combinations as long as they form a clique.  Although we found the microbial combination, how it is associated with the outcome or the directionality is unknown. In the next stage, once the combination of the Taxa is found, we estimate the thresholds of the microbial clique and its association with the outcome. 

## Estimating the thresholds of the microbial clique and its association with the outcome

Run the following function that finds the thresholds for relative abundances of `Taxa.1`, `Taxa.3`, and `Taxa.11`. Each Taxa's directionality is chosen based on its univariate association. The following code _should only be used_ based on the output from the `clique.finder` function; otherwise, overfitting is possible. 

```{}
clique.tba <- function(clique.names, outcome, covariates, grid.quantile, min.prevalence,  data, family = "gaussian"){
  len <- length(clique.names)
  if(len < 2){
    stop("Need at least two exposures to form a meaningful clique") 
  }
  n <- dim(data)[1]
  if(n == 0){
    stop("Please provide a dataset in data.frame format") 
  }
  beta.data <- data.frame(Exposure = rep(NA_character_, len), effect_size = rep(NA_real_, len))
  beta.data$Exposure <- clique.names
  for(i in 1:len){
    g.out <- data[,outcome]
    if(family == "gaussian"){
      fit <- summary(lm(g.out ~ as.matrix(data[, c(beta.data$Exposure[i], covariates)]), data = data))  
    }
    if(family == "binomial"){
      fit <- summary(glm(g.out ~ as.matrix(data[, c(beta.data$Exposure[i], covariates)]), data = data , family = family))  
    }
    if(family == "poisson"){
      fit <- summary(glm(g.out ~ as.matrix(data[, c(beta.data$Exposure[i], covariates)]), data = data , family = family))  
    }
    beta.data$effect_size[i] <- fit$coefficients[2,1]
  }
  x <- grid.quantile
  if(length(grid.quantile) == 1){
    stop("Please increase the number of possible thresholds")
  }
  if(sum(grid.quantile >= 1) != 0){
    stop("All the values provided must be less than 1")
  }
  if(sum(grid.quantile <= 0) != 0){
    stop("All the values provided must be more than 0")
  }
  d1 <- do.call(expand.grid, replicate(len, x, simplify = F))
  d1$min.prevalence <- rep(NA_real_, dim(d1)[1])
  d1$effect_size <- rep(NA_real_, dim(d1)[1])
  d1$se <- rep(NA_real_, dim(d1)[1])
  d1$pvalue <- rep(NA_real_, dim(d1)[1])
  for(i in 1:nrow(d1)){
    mat.len <- as.data.frame(matrix(NA_real_, ncol = len, nrow = nrow(data)))
    for(j in 1:len){
      if(sign(beta.data$effect_size)[j] < 0){
        mat.len[,j] <- as.numeric(data[,clique.names[j]] <= quantile(data[,clique.names[j]], d1[i,j] )) 
      }else{
        mat.len[,j] <- as.numeric(data[,clique.names[j]] >= quantile(data[,clique.names[j]], d1[i,j] ))
      }
    }
    clique.int <- apply(mat.len, 1, function(x){prod(x)})
    
    if(sum(clique.int)!= 0){
      
      d1$min.prevalence[i] <- as.numeric(table(clique.int)/sum(table(clique.int)))[2]
      
      if(d1$min.prevalence[i] >= min.prevalence){
        
        data$clique.int <- clique.int
        g.out <- data[,outcome]
        s <- summary(lm(g.out ~ as.matrix(data[, c("clique.int", covariates)]), data = data))
        d1$effect_size[i] = s$coefficients[2,1]
        d1$se[i] = s$coefficients[2,2]
        d1$pvalue[i] = s$coefficients[2,4]
      }
    }
  }
  d2 <- na.omit(d1)
  if(dim(d2)[1] == 0){
    stop("Please decrease the min.prevalence, but keep it more than 5% for reliability")
  }
  d2 <- d2[d2$min.prevalence > 0.1,]
  d2 <- d2[order(abs(d2$effect_size), decreasing = T),]
  out <- d2[1,]
  
  cuts <-colnames(out)[!(colnames(out) %in% c("min.prevalence", "effect_size", "se",  "pvalue" ))]
  colnames(out)[1:length(cuts)] <- paste0(clique.names, ":Threshold")
  for(i in 1:nrow(beta.data)){
    if(beta.data[i,"effect_size"] < 0){
      out[,i] <- paste0("<=", out[,i]*100,"th Percentile")
    } else {out[,i] <- paste0(">=", out[,i]*100,"th Percentile")}
  }
  return(out)
}
```
Finally, run the `function` called `clique.tba`. Below we discuss each argument for this function and what they entail.

```{}
clique.tba(clique.names = c("Taxa.1", "Taxa.3", "Taxa.11"), outcome= "outcome", covariates = paste0("cov",seq(1,4)),
            grid.quantile = seq(0.2, 0.8, 0.1), min.prevalence = 0.1, family = "gaussian", data = data.simulated)
```
1. `clique.names`: a vector containing the names of most frequently occurring Taxa that form a "clique". We chose `Taxa.1`, `Taxa.3`, and `Taxa.11` based on the output from `clique.finder`.
2. `outcome`: name of the outcome variable
3. `covariates`: a vector containing the names of the covariates
4. `grid.quantile`: choices of the quantiles to search for the thresholds. We removed the lower and upper 20<sup>th</sup> quantiles for stable results. 
5. `min.prevalence`: the minimum proportion (lower bound) of the sample that has the clique. Here we chose `10%` as the lower bound of the prevalence. 
6. `family`: choice of `glm` family of distributions
7. `data`: name of the dataset

I'm sharing below the final output from the simulated example.

```{}
  Taxa.1:Threshold   Taxa.3:Threshold  Taxa.11:Threshold min.prevalence effect_size         se       pvalue
<= 50th Percentile >= 50th Percentile >= 50th Percentile      0.1300813    1.246074 0.06482274 1.150153e-61
```

The `Taxa.1:Threshold`, `Taxa.3:Threshold`, and `Taxa.11:Threshold` denote the estimated thresholds for `Taxa.1`, `Taxa.3`, and `Taxa.11`, respectively. Therefore, this microbial clique is formed in those having (1) Taxa.1 less than 50<sup>th</sup> percentile of the sample, (2) Taxa.3 greater than 50<sup>th</sup> percentile of the sample, and lastly (3) Taxa.11 more than 50<sup>th</sup> percentile of the sample. The recovered estimated effect size is `1.2`, and the estimated prevalence of this microbial clique is almost `13%`. 

### References

1. Basu, S.; Kumbier, K.; Brown, J. B.; Yu, B. Iterative Random Forests to Discover Predictive and Stable High-Order Interactions. Proc. Natl. Acad. Sci. 2018, 115 (8), 1943–1948.
2. Midya V, Lane JM, Gennings C, Torres-Olascoaga LA, Wright RO, Arora M, Tellez-Rojo MM, Eggers S. Prenatal Pb exposure is associated with reduced abundance of beneficial gut microbial cliques in late childhood: an investigation using Microbial Co-occurrence Analysis (MiCA). medRxiv. 2023 May. doi: https://doi.org/10.1101/2023.05.18.23290127.


### Acknowledgments

This method was developed at the Dept. of Environmental Medicine and Public Health, Icahn School of Medicine at Mount Sinai, NYC, with funding and support from the National Institute of Environmental Health Sciences (K99ES032884, P30ES023515, and U2C ES026555-01).

### Contact

Please email Vishal Midya (vishal.midya@mssm.edu) for reporting typos, general questions, and feedbacks.  
