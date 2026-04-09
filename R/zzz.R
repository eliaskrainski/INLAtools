.onAttach <- function(libname, pkgname) {
  if(.Platform$OS.type != "windows") {
    bINLA <- paste0(find.package("INLA"), "/bin/", apINLAbin())
    if((!is.na(file.info("/Library")$isdir)) &&
       (!is.na(file.info("/Applications")$isdir))) {
    ##  Sys.getenv("DYLD_LIBRARY_PATH")
      Sys.setenv(DYLD_LIBRARY_PATH=paste0(bINLA,
                                          ":",
                                          Sys.getenv("DYLD_LIBRARY_PATH"))
                 )
    ##  Sys.getenv("DYLD_LIBRARY_PATH")
    } else {
##      print(Sys.getenv("LD_LIBRARY_PATH"))
      Sys.setenv(LD_LIBRARY_PATH = paste0(bINLA,
                                          ":",
                                          Sys.getenv("LD_LIBRARY_PATH"))
                 )
  ##    Sys.getenv("LD_LIBRARY_PATH")
    }
  }
}

