# this script was designed to count the coverage rate in pancncer level
# Date: 2021-9-10
# Author: Ming Li

library(dplyr)
library(magrittr)
library(stringr)
library(purrr)
library(ggplot2)
library(dplyr)

# get data  ---------------------------------------------------------------

# Time difference of 3.851677 mins

if (1) {
  rm(list = ls())
  start = Sys.time()
  print(start)
  source("script/functions/phenoDataPrep.R", encoding = 'utf-8')
  
  ## loop start -------
  for (panel in c("Our_ori",  "Our",  "Newman", "Burgener", "MSK")) {
    for (Source in c("TCGA")) { 
      message(Source, ' of ', panel)
      resTemp = phenoDataPrep(
        panel = panel,
        type = Source,
        cores = 15,
        verbose = F
      )
      resTemp$Panel = panel
      resTemp$Source = Source
      save(resTemp, file = paste0("Rdata/", panel, "_panel_", Source, "_capRate.Rdata"))
    }
  }
  
  ## loop end ------
  end = Sys.time()
  print(end)
  print(end - start)
  ## turn off the computer
  # system('shutdown -s')
}


