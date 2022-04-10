library(openxlsx)

# Sup Table 1. Patients characteristics
load("Rdata/NPC_cedat.Rdata")
ce_dat %>% dplyr::select(1, 3:13, "Response") %>% 
  mutate(Sex = ifelse(Sex == 0, 'female', 'male')#,
         # SampleName = paste0(Patient, "_", Group) %>% str_replace("-", '_'),
         # filename1 = paste0(Sample, "_R1.fq.gz"),
         # filename2 = paste0(Sample, "_R2.fq.gz")
         ) %>% 
  write.xlsx("supTable/SupTable1.xlsx", overwrite = T)

# Sup Table 2. The genomic regions that MaGIC target
Enricher = openxlsx::read.xlsx("panels/Selector.xlsx") %>% 
  mutate(Region = paste0('Region_', 1:nrow(.)))
Suppressor = data.table::fread("panels/Catcher.csv") %>% 
  select(1:3) %>% 
  mutate(Region = paste0('Region_', 1:nrow(.)))

openxlsx::write.xlsx(x = list(Enricher, Suppressor),
                     file = 'supTable/SupTable2.xlsx', 
                     sheetName = c("MaGIC-Enricher", "MaGIC-Suppressor"),
                     overwrite = T)

# Sup Table 3 . Sequencing quality statistics

file.copy('data/stat/xlsx/reads.stat.xlsx',
          'supTable/SupTable3.xlsx', overwrite = T)



# Sup Table 4. values in NPC ROC plot
openxlsx::read.xlsx("data/NPC_Our_tss.xlsx") %>% 
  openxlsx::write.xlsx("supTable/SupTable4.xlsx")
