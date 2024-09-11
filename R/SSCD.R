#' Sample Specific Compartment Deconvolution
#'
#'
#' \code{SSCD} will perform sample specific gene expression deconvolution on a matrix or data frame of gene expression data using FaStaNMF and gene-specific sample-specific optimization.
#'
#' @param data Gene expression target data, a matrix-like object. The rows should represent genes, and each row must have a unique row name. Each column should represent a different sample.
#'
#' @param compartments_n The factorization rank (number of factors) to be used during NMF. This function argument should be a positive integer value.
#'
#' @param nrun The desired number of FaStaNMF runs to be run
#'
#' @param nmf_seed The desired seed to be used for NMF
#'
#' @param mvg A numerical argument determining how many of the most variable genes to look at during the first steps of FaStaNMF.
#'
#' @param ... Other arguments to be passed to FaStaNMF.
#'
#' @return A list containing W matrix for every sample; in each matrix, rows are genes and columns are deconvolved compartments/factors/tissues types 
#' 
#' @export
#'
#'
#' @import aged 
#'

SSCD <- function(data, compartments_n = 3, nrun = 200, nmf_seed = 124578, mvg = 1000, parallel_n = 2, ...) {

  data_nmf <- aged::fastanmf(data, rank = compartments_n, nrun = nrun, nmf_seed = nmf_seed, mvg = mvg, ...)
  
  
  ##############
  # get h and w
  w = basis(data_nmf)
  h = coef(data_nmf)
  ##############
  
  ##############
  # normalization
  
  # normalize columns of learned W
  median_w = apply(w, 2, median)
  lambda = 100/median_w
  B = diag(lambda)
  B.inv = diag(1/lambda)
  w_new = w %*% B
  h_new = B.inv %*% h
  
  # normalize to have H colSums close to 1
  h_sums <- colSums(h_new)
  mean_sum <- mean(h_sums)
  w <- w_new * mean_sum
  h <- h_new / mean_sum
  ##############
  
  ##############
  # Creates list of Ws for every sample (for NMF, all Ws are the same)
  resultW <- lapply(seq_len(ncol(h)), function(X) w)
  resultH <- apply(h, 2, FUN = sum_to_1)
  ##############
  
  
  ##############################################
  # sets the variables
  samples_n <- ncol(h)
  genes_n <- nrow(data)
  ##############
  
  improvedW <- run_compartment_optimization(data, compartments_n = compartments_n, samples_n = samples_n, genes_n = genes_n, resultH = resultH, resultW = resultW, parallel_n = parallel_n)
  return(improvedW)
}

# helper function                    
sum_to_1 = function(col)
{
  return (col/sum(col))
}
                    

                    
