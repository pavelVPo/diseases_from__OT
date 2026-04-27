library(tidyverse)
library(arrow)
library(RCurl)

## Input
# Input local, SEE: Table S2: Manual Mapping Review in https://www.mdpi.com/1467-3045/45/4/223#app1-cimb-45-00223
manual_evds <- read_tsv("C:/.../manual_evds.tsv") |>
					separate_longer_delim(ID, delim = "_")
# Input FTP
raw_evds <- getURL("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/21.04/output/etl/json/evidence/sourceId%3Deuropepmc/",
	                 verbose=TRUE,ftp.use.epsv=FALSE, dirlistonly = TRUE, crlf = TRUE) |>
					str_match_all("part.*json") |>
					unlist() |>
					unique() |>
					str_replace('".*', "")
# Raw FTP data for manual_evds
for (i in seq(1:length(raw_evds))) {
	assign(str_glue("data_{i}"), read_json_arrow(str_glue("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/21.04/output/etl/json/evidence/sourceId%3Deuropepmc/{raw_evds[i]}")) |> semi_join(manual_evds, by = c("id" = "ID")))
}

## Process
# select data to work with
data_names <- tibble(var = ls.str(envir = .GlobalEnv)) |> filter(str_starts(var, "data_[0-9]")) |> pull(var)
data_expr  <- str_glue("data_all <- bind_rows({str_c(data_names, collapse = ', ')})") 
eval(str2expression(data_expr))
# filter
data_all <- data_all |> inner_join(manual_evds, by = c("id" = "ID")) |> filter(`Does the text fragment allow to make a judgement on the category of experimental method?` == 1 &
								`Is dictionary-defined level study mentioned in the text fragment?` == 1 &
								(`Category by the Dictionary` == "Nucleic Acids" | `Category by the Dictionary` == "Proteins"))
disease_nucl <- data_all |> filter(`Category by the Dictionary` == "Nucleic Acids") |> select(diseaseLabel) |> mutate(nucl = 1) |>
								group_by(diseaseLabel) |>
								summarize(nucl = sum(nucl)) |>
								ungroup()
disease_prot <- data_all |> filter(`Category by the Dictionary` == "Proteins") |> select(diseaseLabel) |> mutate(prot = 1) |>
								group_by(diseaseLabel) |>
								summarize(prot = sum(prot)) |>
								ungroup()
disease = inner_join(disease_nucl, disease_prot, by = "diseaseLabel") |>
					rowwise() |>
					mutate(total = nucl + prot, diff = max(nucl, prot) - min(nucl, prot)) |>
					ungroup() |>
					arrange(desc(diff))

# Glance:
head(disease)