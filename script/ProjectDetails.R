library(dplyr)
library(stringr)

world_cancer = c("Lung", "Colorectal", "Liver", "Stomach", "Breast", 
                 "Esophagus", "Pancreas", "Bladder", "Uterus", "Blood")
women_cancer = c("Breast", "Cervix","Ovary", "Uterus")

pro_df = openxlsx::read.xlsx("data/proj_details.xlsx") %>% as_tibble() %>% 
  mutate(World = if_else(Primary.Site %in% world_cancer, "Y", "N"),
         China = if_else(Country == "China", "Y", "N"),
         Women = if_else(Primary.Site %in% women_cancer, "Y", "N")) %>% 
  dplyr::filter(str_detect(Project.Name, "TCGA", negate = T))
save(pro_df, file = "Rdata/project_details.df.Rdata")