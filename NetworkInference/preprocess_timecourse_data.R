## Preprocessing of time-course data
## Using comparison to unexposed over time
library(readr)
library(here)
library(dplyr)


## Original time course (NB: Not included in repository)
timecourse_data = read_tsv(here('Data/Combined log fold changes all vs unexposed.txt'))

timecourse_matrix = as.matrix(timecourse_data[,2:11])
observation_times = c(0, 0.5, 1, 2, 4, 8, 12, 24, 48, 72)
gene_names = timecourse_data[,1]

saveRDS(list(gex_fc=timecourse_matrix, times=observation_times,
             gene_names=gene_names), file=here('Data/gex_fc_timecourse.rds'))

## Pre-processing for ageing time course
ageing_times = c(8, 12, 24, 48, 72)
ageing_files = paste0(here('Data/'), ageing_times, 'h - Unexposed.txt')

# Get probe names
ageing_probes = (read_tsv(ageing_files[1]) %>%
  arrange(SystematicName))$SystematicName

# Extract log fold change data
ageing_matrix = sapply(ageing_files, function(file) {
    ageing_data = read_tsv(file) %>%
      arrange(SystematicName)
    
    logfc_data = matrix(ageing_data$logFC)
    
    return(logfc_data)
  })

colnames(ageing_matrix) = paste0(ageing_times, 'hourFC')

saveRDS(list(gex_fc=ageing_matrix, times=ageing_times,
             gene_names=ageing_probes), file=here('Data/ageing_fc_timecourse.rds'))

