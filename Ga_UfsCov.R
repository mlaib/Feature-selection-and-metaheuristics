#' Ga_UfsCov: Unsupervised Feature Selection Via Space-Filling Design and Genetic Algorithm 
#'
#' Applies feature selection algorithm based on the space filling concept,
#' and using Genetic Algorithm (GA) as search strategy.
#' @usage Ga_UfsCov(df, nBits=NULL, pmutation=0.02, maxiter=50, popSize=50, pcrossover=0.8)
#' @param df Data of class: \code{matrix} or \code{data.frame}.
#' @param nBits An integer, number of variable (parameter from GA library).
#' @param pmutation Probability of mutation in a parent chromosome (parameter from GA library).
#' @param maxiter Maximum number of iterations  (parameter from GA library).
#' @param popSize Population size  (parameter from GA library).
#' @param pcrossover Probability of crossover between pairs of chromosomes (parameter from GA library).
#'
#' @return A list of two elements:
#'  \itemize{
#'   \item \code{BestSolution} a binary vector containing the selected features.
#'   \item \code{GA} A ga-class (See more information on the GA package).
#'   }
#' @author Mohamed Laib \email{Mohamed.Laib@@list.
#'
#' @note The work is presented in a Conference. ROADEF 2022: 23ème édition du congrès annuel de 
#' la Société Française de Recherche Opérationnelle et d'Aide à la Décision. In this work, 
#' the algorithm was applied to steel manufacturing data (provided by ArcelorMittal). The work is
#' part of the PAX project (FNR-Bridge).
#' 
#' @details The algorithm does not deal with missing values and constant
#' features. Please make sure to remove them before applying it.
#' 
#' Since the algorithm is based on pairwise distances, and
#' according to the computing power of your machine, large number of
#' data points can take much time and needs more memory.
#' 
#' The evaluation measure (Coverage measure) does not need to tune parameters. the used parameters in
#' this function are needed mainly for the search strategy based on genetic algorithm (provided by
#' GA R library).
#'
#' @references
#' Riad Aggoune, Mohamed Laib, A Genetic Algorithm for Feature Selection Applied to Data 
#' From Multiples Sources: Application to Manufacturing Data, Recherche Opérationnelle et 
#' d'Aide à la Décision, Roadef 2022.
#' 
#' M. Laib, M. Kanevski, 
#' \href{https://www.elen.ucl.ac.be/Proceedings/esann/esannpdf/es2018-57.pdf}{A novel 
#' filter algorithm for unsupervised feature selection based on a space filling measure}. 
#' Proceedings of the 26rd European Symposium on Artificial Neural Networks, Computational 
#' Intelligence and Machine Learning (ESANN), pp. 485-490, Bruges (Belgium), 2018.
#' 
#' M. Laib and M. Kanevski, A new algorithm for redundancy minimisation in 
#' geo-environmental data, 2019.
#' \href{https://www.sciencedirect.com/science/article/pii/S0098300418310975}{Computers & 
#' Geosciences, 133 104328}.
#' 
#' Scrucca L. (2013). GA: A Package for Genetic Algorithms in R. Journal of Statistical 
#' Software, 53(4), 1-37, doi: 10.18637/jss.v053.i04.
#' 

### Loading libraries ####
library(Biobase) 
library(GA)  
library(dplyr) 

### Internal function: evaluation measure ####
EvalMeasure <- function (design) {
  X <- as.matrix(design)
  n <- nrow(X)
  Dmin <- as.matrix(matchpt(X)[, 2])
  gammabar <- mean(Dmin)
  s <- sum(apply(Dmin, 1, FUN = function(a) ((a - gammabar)^2)))
  cov <- (1/gammabar) * ((1/n) * s)^(1/2)
  return(cov)
}

### Main function #####
Ga_UfsCov <- function(df, nBits=NULL, pmutation=0.02, maxiter=50, popSize=50, 
                      pcrossover=0.8){
  EvalFunc <- function(x) {
    cols <- colnames(df)[x == 1]
    new_df <- df %>% select(one_of(cols))
    NewDim <- -EvalMeasure(new_df) 
    return(NewDim)
  }
  nBits = ifelse(is.null(nBits), ncol(df), nBits)
  GA <- ga(type = "binary", fitness = EvalFunc, nBits = nBits, pmutation = pmutation, 
           maxiter = maxiter, popSize = popSize, pcrossover = pcrossover)
  a <- GA@solution
  colnames(a) <- colnames(df)
  
  return(list(BestSolution=a, GA=GA))
}
