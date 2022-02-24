# Feature-selection-and-metaheuristics
Feature selection using diferent search strategies. 

### Notes
This is a first version of an in-progress work. more documentations, references, and improvements will be added later.

## Version
1.0

## Required packages

```{r}
install.packages(GA)
install.packages(dplyr)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biobase")
```

### Example

```{r}
# Example on the Butterfly dataset available in IDmining R library
df <- IDmining::Butterfly(1000)[,-9]
res <- Ga_UfsCov(df, nBits=NULL, pmutation=0.02, maxiter=10, popSize=100, 
                      pcrossover=0.8)
res$BestSolution
```
