#Initializing functions

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to the CytoTRACE R package, a tool for the unbiased prediction of differentiation states in scRNA-seq data. For more information about this method, visit https://cytotrace.stanford.edu.\n")

  warn_if_no_scanoramaCT <- function() {
    have_scanoramaCT <<- reticulate::py_module_available("scanoramaCT")
    have_numpy <<- reticulate::py_module_available("numpy")
    if (!have_scanoramaCT)
      warning("The ScanoramaCT python module is not accessible. The iCytoTRACE function for integration across multiple datasets will be disabled. Please follow the instructions in https://github.com/gunsagargulati/CytoTRACE to install the necessary Python packages for this application.", call. = FALSE)
    if(!have_numpy)
      warning("The numpy python module is not accessible. The iCytoTRACE function for integration across multiple datasets will be disabled. Please follow the instructions in https://github.com/gunsagargulati/CytoTRACE to install the necessary Python packages for this application.", .call = FALSE)
  }
  warn_if_no_scanoramaCT()
}






