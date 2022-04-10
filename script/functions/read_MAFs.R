#' Integrate MAFs into one file
#'
#' @param samples a character vector contains sample names
#' @param path_prefix the path prefix in front of 'sap'.hg38.annovar.txt
#'
#' @return an integrated MAF file
#' @export
#'
#' @examples
read_MAFs = function(samples, 
                     path_prefix, 
                     multiThread=FALSE, 
                     cores = 10, 
                     verbose = TRUE) {
  
  ## using multi threads
  if(multiThread){
    
    library(foreach)
    library(doParallel)
    
    cl = makeCluster(cores)
    registerDoParallel(cl)
    
    maf = foreach(
      sap = samples,
      .combine = plyr::rbind.fill,
      .packages = c(
        "plyr",
        "maftools",
        "rlang"
      ),
      .errorhandling = 'stop',
      .verbose = verbose
    ) %dopar% {
      path = paste0(path_prefix, '/', sap, ".hg38.annovar.txt")
      
      if (!file.exists(path))
        stop("Error, Please check MAF file path!")
      
      a_maf = maftools::annovarToMaf(path)
      a_maf$Tumor_Sample_Barcode = sap
      
      return(a_maf)
    }
    
    stopImplicitCluster()
    stopCluster(cl)
    
   
    
  }else {
    ## don't use multi threads
    maf = NULL
    for (sap in samples) {
      path = paste0(path_prefix, '/', sap, ".hg38.annovar.txt")
      
      if (!file.exists(path))
        stop("Error, Please check MAF file path!")
      
      a_maf = maftools::annovarToMaf(path)
      a_maf$Tumor_Sample_Barcode = sap
      maf = rbind(maf, a_maf)
    }
  }
  
  return(maf)
}