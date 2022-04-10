library(dplyr)

ce_dat = openxlsx::read.xlsx('data/NPC_clin.xlsx') %>% 
  mutate(
    Patient = str_replace(Patient, pattern = '_', replacement = ''),
    Group = if_else(Group == "Before", "Pre-therapy", "Post-therapy"),
    Response = if_else(Response == "PR",  "Responder", "Nonresponder")
  ) %>%
  dplyr::mutate(
    Group = factor(
      Group,
      levels = c("Pre-therapy", "Post-therapy"),
      ordered = T
    ),
    Response = factor(
      Response,
      levels = c("Responder", "Nonresponder"),
      ordered = T
    )
  ) %>% as_tibble()

save(ce_dat, file = 'Rdata/NPC_cedat.Rdata')