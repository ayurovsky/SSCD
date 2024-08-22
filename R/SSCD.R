#' 
#' @export
#'
#'
#' @import aged 
#'

SSCD <- function() {
  result <- tryCatch(
  {
		aged::aged(NULL)
  },
  	error = function(e) {
    print("An error occurred: ", e$message)
    NULL 
  	}
	)
	return("After the error")
}
