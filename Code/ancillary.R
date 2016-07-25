#' This files provide basic functions to load multiple
#' libraries, source all the files in directory, ...


#' `sourceDir()` sources all the `.R` files in a directory
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[Rr]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

#' `try_and_install()` testa for the presence of a package in the library
#' and, if not present, installa it.
try_and_install <- function(package_names){
  installed_packages <- installed.packages()
  installed <- character(0)
  for(Name in package_names){
    if(!Name %in% installed_packages){
      installed <- c(installed,Name)
      install.packages(Name,
                       verbose = F)
    }
  }
  if(length(installed) == 0){return("Nothing to install.")}
  installed_text <- paste(installed,
                          sep="\n")
  return(paste("All packages installed:",installed_text))
}