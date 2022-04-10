#' Title
#'
#' @param panel a panel to capture
#' @param type data source
#' @param cores threading cores number
#' @param verbose bool
#' @param exons the exons in panel used. Default NULL means all
#'
#' @return a data.frame contains phenotype data of patients.
#' @export
#'
#' @examples

# # for test
# library(dplyr)
# library(magrittr)
# library(GenomicRanges)
# library(stringr)
# source("script/functions/Map2Selector.R")
# source("script/functions/Polisher.R")
# a = 1
# panel="Our"
# prt = ifelse(a, "ACC", "ALL-US")
# cores = 10
# type = ifelse(a, "TCGA", "ICGC")
# verbose = TRUE
# exons = NULL
phenoDataPrep = function(panel = c("Our", "Newman", "Burgener", "MSK", "Our_ori"),
                         type = c('TCGA', "ICGC"),
                         cores = 10,
                         verbose = TRUE,
                         exons = NULL) {
  library(foreach)
  library(doParallel)
  
  # param prep
  panel = match.arg(panel)
  type = match.arg(type)
  cl = makeCluster(cores)
  registerDoParallel(cl)
  
  
  load("Rdata/project_details.df.Rdata")
  projects = pro_df$Project[pro_df$Source == type]
  
  # do parallel  ------------------------------------------------------------
  
  results = foreach(
    prt = projects,
    .export = exons,
    .combine = plyr::rbind.fill,
    .packages = c(
      "plyr",
      "stringr",
      "magrittr",
      "maftools",
      "rlang",
      "GenomicRanges",
      "dplyr"
    ),
    .errorhandling = 'stop',
    .verbose = verbose
  ) %dopar% {
    # data prep ---------------------------------------------------------------
    source("script/functions/Map2Selector.R")
    source("script/functions/Polisher.R")
    library(dplyr)
    
    ###########################
    ## for TCGA data loading ##
    ###########################
    ### set ref version
    grHg = "hg38"
    
    ### load os data
    os_dat = data.table::fread(paste0(
      "D:/Data/xena/TCGA-",
      prt,
      "/phenotype/",
      prt,
      "_survival.txt.gz"
    )) %>%
      dplyr::select(contains('patient'), OS, OS.time, PFI, PFI.time) %>%
      magrittr::set_colnames(c("PatientID", "OS", "OS.time", "PFS", "PFS.time")) %>%
      dplyr::distinct()
    
    ### if has, load TNM etc., else, donot load it.
    namecol = c("Age",
                "Sex",
                "T.stage",
                "N.stage",
                "M.stage",
                "Stage")
    phe_dat = data.table::fread(paste0(
      "D:/Data/xena/TCGA-",
      prt,
      "/phenotype/",
      prt,
      "_clinicalMatrix"
    ))
    try({
      os_dat = dplyr::left_join(
        os_dat,
        phe_dat %>% dplyr::select(PatientID = `_PATIENT`,
                                  Age = age_at_initial_pathologic_diagnosis),
        by = "PatientID"
      )
    })
    try({
      os_dat = dplyr::left_join(os_dat,
                                phe_dat %>% dplyr::select(PatientID = `_PATIENT`,
                                                          Sex = gender),
                                by = "PatientID")
    })
    try({
      os_dat = dplyr::left_join(
        os_dat,
        phe_dat %>% dplyr::select(PatientID = `_PATIENT`,
                                  T.stage = ends_with("_T")[1]),
        by = "PatientID"
      )
    })
    try({
      os_dat = dplyr::left_join(
        os_dat,
        phe_dat %>% dplyr::select(PatientID = `_PATIENT`,
                                  N.stage = ends_with("_N")[1]),
        by = "PatientID"
      )
    })
    try({
      os_dat = dplyr::left_join(
        os_dat,
        phe_dat %>% dplyr::select(PatientID = `_PATIENT`,
                                  M.stage = ends_with("_M")[1]),
        by = "PatientID"
      )
    })
    try({
      os_dat = dplyr::left_join(
        os_dat,
        phe_dat %>% dplyr::select(PatientID = `_PATIENT`,
                                  Stage = ends_with("_stage")[1]),
        by = "PatientID"
      )
    })
    ### add NAs if subset is not success
    out = namecol[!(namecol %in% colnames(os_dat))]
    if (length(out) != 0) {
      out_df = data.frame(matrix(NA, nrow = nrow(os_dat), ncol = length(out))) %>%
        magrittr::set_colnames(out)
      os_dat = cbind(os_dat, out_df) %>%
        dplyr::select("PatientID",
                      "OS",
                      "OS.time",
                      "PFS",
                      "PFS.time",
                      all_of(namecol))
    }
    
    ### return NULL if os_dat is NULL, for stop this thresholds
    if (nrow(os_dat) == 0)
      return(NULL)
    
    ### load maf
    maf_path = paste0("Rdata/Rmaf/", type, '-', prt, '_maf.Rdata')
    if (file.exists(maf_path)) {
      load(maf_path)
    } else {
      maf_dat = maftools::read.maf(paste0("D:/Data/TCGA_maf/TCGA.", prt, ".somatic.maf.gz"),
                                   verbose = F)@data
      save(maf_dat, file = maf_path)
    }
    ### format and subset maf
    maf_dat = maf_dat %>% as.data.frame() %>%
      dplyr::mutate(VAF = t_alt_count / t_depth * 100) %>%
      dplyr::select(contains(
        c(
          'Tumor_Sample_Barcode',
          "chr",
          "start",
          "end",
          "Tumor_Seq_Allele2",
          "VAF"
        ),
        ignore.case = T
      )) %>%
      magrittr::set_colnames(c("Sample", "Chr", "Start", "End", "Allele", "VAF")) %>%
      dplyr::mutate(PatientID = str_sub(Sample, 1, 12)) %>%
      dplyr::distinct() #%>% filter(PatientID %in% os_dat$PatientID)
    ### format to gr
    maf_gr = GRanges(
      seqnames = maf_dat$Chr,
      ranges = IRanges(start = maf_dat$Start, end = maf_dat$End),
      PatientID = maf_dat$PatientID,
      ALT = maf_dat$Allele,
      VAF = maf_dat$VAF
    )
    
    
    
    
    # Polishing for our panel only -------------------------------------------
    if (panel == "Our") {
      index = Polisher(maf_gr, hg = grHg)
      maf_gr = maf_gr[index]
    }
    
    
    # get capture data -------------------------------------------------------
    maf_pat = data.frame(PatientID = maf_gr$PatientID %>% unique)
    temp_dat = Map2Selector(maf_gr,
                            panel = panel,
                            grHg = grHg,
                            exons = exons) %>%
      as.data.frame()
    
    
    # Wed Dec 01 20:48:25 2021 ------------------------------
    
    ## concat it to other variables
    cap_num = temp_dat %>%
      dplyr::group_by(PatientID) %>%
      dplyr::summarise(
        SNVs = length(Region),
        NumMutRegion = length(Region %>% unique),
        mVAF = mean(VAF)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::right_join(maf_pat, by = "PatientID")
    
    
    cap_num[is.na(cap_num)] = 0
    
    prt_df = dplyr::right_join(os_dat, cap_num, by = "PatientID")
    
    ## add project id
    prt_df$Project = prt
    
    return(prt_df)
    
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  #print(results)
  return(results)
}