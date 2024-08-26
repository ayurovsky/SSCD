#'
#' @export
#'
#'
#' @import ROI
#' @import ROI.plugin.alabama
#' @import foreach
#' @import doParallel
#'
#'


run_compartment_optimization <-  function(trueMixed, compartments_n, samples_n, genes_n, resultH, resultW, parallel_n=2) {

  registerDoParallel(parallel_n)

  print("Starting Sample Specific Compartment Optimization...")
  start_time <- Sys.time()
  improvedW <- foreach (sample_n=1:samples_n) %dopar% {
    sample_w <- vector()
    for (gene_n in 1:genes_n) {
      mixed <- trueMixed[gene_n,sample_n]
      nmf_h <- resultH[,sample_n]
      nmf_w <- resultW[[sample_n]][gene_n,]

      # build objective function
      eval_f <- "function(x) { return ( "
      for (i in 1:(compartments_n-1)){
         eval_f <- paste0(eval_f, "((x[", i ,"] - ", nmf_w[i], ")/", nmf_w[i], ")^2 + ")
      }
      eval_f <- paste0(eval_f,"((x[",i+1,"] - ", nmf_w[i+1], ")/", nmf_w[i+1], ")^2 ) }")
      eval_f <- eval(parse(text = eval_f)) #noquote(eval_f)

      # equality constraints
      eval_g_eq <- "function(x) { return ( "
      for (i in 1:(compartments_n-1)){
        eval_g_eq <- paste0(eval_g_eq, "x[", i ,"] * ", nmf_h[i], " + ")
      }
      eval_g_eq <-  paste0(eval_g_eq, "x[", i+1 ,"] * ", nmf_h[i+1], " - ", mixed, ") }")
      eval_g_eq <- eval(parse(text = eval_g_eq)) #noquote(eval_g_eq)

      # lower bound constraints

      #initial values
      x0 <- as.vector(nmf_w)

      nlp <- OP(objective = F_objective(eval_f, n=compartments_n),
              constraints = F_constraint(eval_g_eq, dir = "==", rhs = 0), bounds = V_bound(lb = rep(0.0001, compartments_n)))

      sol <- ROI_solve(nlp, solver = "alabama", start = x0)

      neww <- as.vector(sol$solution)

      if (sol$status$msg$symbol != "SUCCESS") {
        neww <- nmf_w
      }

      sample_w <- rbind(sample_w, neww)
      rownames(sample_w)[gene_n] <- rownames(trueMixed)[gene_n]
    }
    sample_w
  }
  end_time <- Sys.time()
  print(end_time - start_time)

  return(improvedW)

}


