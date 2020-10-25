# Simple script to concatenate the 1500 files from the network inference.
library(readr)
library(here)
library(dplyr)

file_names = list.files(here('Results/'), 'bn_marg_probs_bipartite_',
                        full.names=TRUE)

# Load in results
results_table_list = lapply(file_names, function(file) {
                                          read_csv(file) %>% 
                                            mutate(smoothing=ifelse(grepl('loess', file), 
                                                                    'loess', 'none'),
                                                   age_correct=ifelse(grepl('age_correction', file), 
                                                                    'corrected', 'uncorrected'),
                                                   circ_correct=ifelse(grepl('no_circ', file), 
                                                                      'uncorrected', 'corrected'))})

results_table = do.call('rbind', results_table_list)

write_csv(results_table, here('Results/full_results_bn_marg_probs.csv'))