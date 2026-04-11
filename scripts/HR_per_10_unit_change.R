rm(list = ls())

set.seed(123)

library(docxtractr)
library(dplyr)
library(tidyr)
library(officer)
library(flextable)

# Set global defaults for all following flextables
set_flextable_defaults(
  font.family = "Times New Roman",
  line_spacing = 1,
  padding.top = 0,
  padding.bottom = 0
)

#  word_sect_properties
sect_properties <- prop_section(
  page_size = page_size(
    orient = "landscape",
    height = 11.7,
    width = 8.3
  ),
  page_margins = page_mar(
    top = 0.5,
    bottom = 0.5,
    left = 0.5,
    right = 0.5,
    header = 0.3,
    footer = 0.3,
    gutter = 0
  ),
  type = "continuous"
)

tenUnitHR <- function(HR){
  beta <- log(HR)
  beta <- beta*10
  return(exp(beta))
}

doc_filename <- "Supplementary_tables .docx"

  
filename <- file.path(
  "/Users",
  "debda22",
  "Projects",
  "T1D_projects",
  "REVISION_10_04_2026",
  doc_filename
)

# Read the Word document
doc <- docxtractr::read_docx(path = filename)

# List all tables in the document
docx_describe_tbls(doc)

# Table Name	Description	Data
# Table T1	Baseline statistics. P-values from chi-squared test.	
# 
# 
# Full follow-up
# Table T2	Scaled Schoenfeld residuals with time	
# Table T3	Cox proportional hazards models with Firth's Penalized Likelihood: conc. 25(OH)D	
# Table T4	Cox proportional hazards models with Firth's Penalized Likelihood: DBP	
# Table T5	Conditional regression: stratified interaction model with child covariate balancing	
# Table T6	Conditional regression: stratified interaction model with maternal covariate balancing	
# 
# Table T7	Cox proportional hazards models with Firth's Penalized Likelihood: conc. 25(OH)D	Up to age ~5 years
# Table T8	Cox proportional hazards models with Firth's Penalized Likelihood: DBP	
# 
# Table T9	Cox proportional hazards models with Firth's Penalized Likelihood: conc. 25(OH)D	Up to age ~10 years
# Table T10	Cox proportional hazards models with Firth's Penalized Likelihood: DBP	
# 
# Table T11	Cox proportional hazards models with Firth's Penalized Likelihood: conc. 25(OH)D	Up to age ~15 years
# Table T12	Cox proportional hazards models with Firth's Penalized Likelihood: DBP	
# 
# Table T13	Cox proportional hazards models with Firth's Penalized Likelihood: conc. 25(OH)D	Boy-only
# (full follow-up)
# Table T14	Cox proportional hazards models with Firth's Penalized Likelihood: DBP	
# 
# Table T15	Cox proportional hazards models with Firth's Penalized Likelihood: conc. 25(OH)D	Girl-only
# (full follow-up)
# Table T16	Cox proportional hazards models with Firth's Penalized Likelihood: DBP	

table_number <- 15 + 1

# Extract a specific table (e.g., the first table)
curr_table <- docx_extract_tbl(doc, tbl_number = table_number) 

curr_table <- curr_table %>%
  tidyr::separate_wider_regex(
    cols = `HR..95..CI.`,
    patterns = c(
      HR = ".*",              # Get the Hazard Ratio
      " \\(",                 # Match the space and opening parenthesis (not saved)
      CI_low = ".*",          # Get the lower bound
      "-",                    # Match the dash (not saved)
      CI_high = ".*",         # Get the upper bound
      "\\)"                   # Match the closing parenthesis (not saved)
    )
  )  %>%
  mutate(across(c("HR","CI_low","CI_high"), as.numeric))

curr_table <- curr_table %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    HR_new = ifelse(Level == "+1", tenUnitHR(HR), HR),
    CI_low_new = ifelse(Level == "+1", tenUnitHR(CI_low), CI_low),
    CI_high_new = ifelse(Level == "+1", tenUnitHR(CI_high), CI_high)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(`NEW HR (95% CI)` = paste0(
    round(HR_new, 3),
    " (",
    round(CI_low_new, 3),
    "-",
    round(CI_high_new, 3),
    ")"
  ))

curr_table <- curr_table[,c(1:3,11)] %>%
  dplyr::rename(`HR (95% CI)` = `NEW HR (95% CI)` ) %>%
  dplyr::mutate(Level = gsub("\\+1", "per 10 nmol/L increase", Level))

outDir <- "/Users/debda22/Projects/T1D_projects/REVISION_10_04_2026"

if(!dir.exists(outDir)){
  dir.create(outDir, recursive = TRUE)
} else {
  message("Directory exists!")
}

output_doc_filename <- "Supplementary Table T15.docx"
  
# Save the table
flextable::save_as_docx(
  curr_table %>%
    flextable::flextable() %>% 
    flextable::autofit() %>% 
    flextable::set_table_properties(width = 1, layout = "autofit"),
  path = file.path(outDir,output_doc_filename),
  pr_section = sect_properties
)
