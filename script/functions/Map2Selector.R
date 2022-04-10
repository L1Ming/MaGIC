#' Map GRanges to different panels
#'
#' @param gr GRanges object to map selector
#' @param grHg The reference version, if it isn't hg38, then transform it to hg38
#' @param panel The panel of selector to use
#' @param exons The exons of a panel to use, default all
#'
#' @return GRanges object, mapped gr
#' @export
#'
#' @examples
Map2Selector = function(gr = NULL,
                        grHg = c("hg38", "hg19"),
                        panel = c("Our", "Newman", "Burgener", "MSK", "Our_ori"),
                        exons = NULL) {
  library(GenomicRanges)
  library(dplyr)
  library(stringr)
  
  ## input data check ---------------------------
  
  grHg = match.arg(grHg)
  panel = match.arg(panel)
  
  # load panel ----------------------------
  if (panel == "Our_ori") {
    # primary selector
    selector = openxlsx::read.xlsx("panels/Selector.xlsx") %>%
      dplyr::transmute(Chr = Chr,
                       Start = Start.bp,
                       End = End.bp)
    # selector %>% as.data.frame() %>%
    #   mutate(a = ".", b = ".", c = ".") %>%
    #   write.table('data/Enricher.txt', quote = F,
    #               sep = '\t', row.names = F, col.names = F)
    
    ## hg38
    Ref = 'hg38'
  } else if (panel == "Our"){
    ## optimal selector
    selector = data.table::fread("panels/Catcher.csv") %>%
      dplyr::transmute(Chr = seqnames,
                       Start = start,
                       End = end)
    
    # selector %>% as.data.frame() %>%
    #   mutate(a = ".", b = ".", c = ".") %>%
    #   write.table('data/Suppressor.txt', quote = F,
    #               sep = '\t', row.names = F, col.names = F)
    ## hg38
    Ref = 'hg38'
    
  } else if (panel == "Newman") {
    ## hg19
    selector = openxlsx::read.xlsx("panels/Newman.xlsx", 2)
    Ref = 'hg19'
  } else if (panel == "Burgener") {
    ## hg19
    selector = openxlsx::read.xlsx("panels/Burgener.xlsx", 5)
    Ref = 'hg19'
  } else if (panel == "MSK") {
    ## hg19
    selector = openxlsx::read.xlsx("panels/MSK.xlsx", startRow = 4) %>% 
      dplyr::rename(Start = start, End = stop) %>% 
      mutate(Chr = paste0('chr', Chr))
    Ref = "hg19"
  }
  
  # transform to hg38 ----------------------------
  ## define function
  source('script/functions/transHg38.R')

  
  ## transform
  se_gr = GRanges(
    seqnames = selector$Chr,
    ranges = IRanges(start = selector$Start, end = selector$End),
    Region = paste0("Exon_", 1:nrow(selector))
  ) %>% transHg38(from = Ref)
  
  
  if(is.null(gr)) {
    return(se_gr)
  } else {
    
    if ((str_sub(seqnames(gr)[1] %>% as.character(), 1, 1) != "c"))
      stop("The chrnosome must start with 'chr'.")
    
    gr = transHg38(gr, from = grHg)
    ## subset by indicated exons ----------------------
    if (!is.null(exons)) {
      se_gr = se_gr[se_gr$Region %in% exons]
    }
    
    # overlap and subset -------------------
    overlap = findOverlaps(query = se_gr, subject = gr) %>% as.data.frame()
    grCap = gr[overlap$subjectHits]
    grCap$Region = se_gr$Region[overlap$queryHits]
    grCap$RegionWidth = width(se_gr[overlap$queryHits])
    
    ## cat
    message(
      "captured SNVs is ",
      length(grCap),
      " with used regions in Selector=",
      length(se_gr),
      "\n"
    )
    
    return(grCap)
  }
}