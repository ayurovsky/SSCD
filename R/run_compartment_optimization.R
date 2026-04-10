#' Run optimization for individual samples
#'
#'
#' \code{run_compartment_optimization} will perform gene-specific sample-specific optimization on the output of FastaNFM or other deconvolution technique
#'
#' @param data Gene expression target data, a matrix-like object. The rows should represent genes, and each row must have a unique row name. Each column should represent a different sample.
#'
#' @param compartments_n The factorization rank (number of factors) to be used during NMF. This function argument should be a positive integer value.
#'
#' @param samples_n The number of samples in the dataset
#'
#' @param genes_n The number of genes in the dataset
#'
#' @param resultH The H matrix from (FaStaNMF) deconvolution
#'
#' @param parallel_n The number of cores available for parallel optimization - will drastically impact runtime
#'
#' @return A list containing W matrix for every sample; in each matrix, rows are genes and columns are deconvolved compartments/factors/tissues types 
#' 
#'
#' @export
#'
#'
#' @import foreach
#' @import doParallel
#'
#'


run_compartment_optimization <-  function(data, compartments_n, samples_n, genes_n, resultH, resultW, parallel_n=2) {

  registerDoParallel(parallel_n)

  print("Starting Sample Specific Compartment Optimization...")
  start_time <- Sys.time()

  closed_form_W <-  foreach (sample_n=1:samples_n) %dopar% {
    sample_w <- vector()
    for (gene_n in 1:genes_n) {
      mixed <- mixture[gene_n,sample_n]
      nmf_h <- resultH[,sample_n]
      nmf_w <- resultW[[sample_n]][gene_n,]

      # initialize not-to-use-indexes
      not_to_use <- integer(0)
      iterate <- TRUE

      while(iterate) {
	iterate <- FALSE # ideally this loop will run just once
	# calculate the closed form solution
	fh_sum <- 0
	fh_sum_sq <- 0
	for (i in 1:compartments_n) {
	  if (!(i %in% not_to_use)) { # ignore compartments set to zero in previous iteration(s)
	    fh_sum <- fh_sum + nmf_h[i]*nmf_w[i]
	    fh_sum_sq <- fh_sum_sq + nmf_h[i]*nmf_w[i]*nmf_h[i]*nmf_w[i]
	  }
	}
	common <- (mixed - fh_sum)/fh_sum_sq
	closed_w <- c()
	for (i in 1:compartments_n){
	  if (i %in% not_to_use) { # ignore compartments set to zero in previous iteration(s)
	    closed_w <- c(closed_w, 0.0)
	  } else {
	    new_val <- nmf_w[i] + nmf_h[i]*nmf_w[i]*nmf_w[i]*common
	    if (new_val < 0.0) {
	      not_to_use <- c(not_to_use, i)
	      iterate <- TRUE # negative, set to zero, will need to iterate
	      closed_w <- c(closed_w, 0.0)
	    } else {
	      closed_w <- c(closed_w, new_val) # adding the solved value for this compartment
	    }
	  }
	}
      }
      sample_w <- rbind(sample_w, closed_w)
    }
    sample_w
  }
  end_time <- Sys.time()
  print(end_time - start_time)

  return(closed_form_W)

}


