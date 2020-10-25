# The following code applies the EDISON package to infer networks from a subset 
# of genes in an insecticide resistance timecourse. In this case, the subset is 
# selected from the transcription factors with large interaction effects in a 
# graphical Gaussian model of a separate, spatial gene expression dataset.
library(here)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(EDISON)

source(here('Code/util.R'))

loess_flag = TRUE # Use LOESS interpolation to obtain equally-spaced observations
circadian_correction = TRUE # Use circadian correction to avoid picking up 
# circadian rhythms
ageing_flag = TRUE # Use data from an experiment without insecticide exposure 
# to avoid picking up changes due to ageing

time_vec = c(0,0.5,1,2,4,8,12,24,48,72) # Time points (in hours)
ageing_time_vec = c(8,12,24,48,72) # Time points for ageing experiment

# Read in data (note, data not currently included in repository)
dataset = readRDS(here('Data/gex_fc_timecourse.rds'))

# List of transcription factors and target genes to include
full_names_tf = read_table(here('Data/TF_list.txt'), col_names=FALSE)[[1]]
full_names_targets = read_tsv(here('Data/Target_list.txt'), col_names=TRUE)[[1]]


# Loop over target genes (NB: If running on a multi-core machine, we recommend
# to parallelise this loop)
for(i in 1: length(full_names_target)) {
  full_names = c(full_names_tf, full_names_targets[i])
  
  gene_types = gsub('.*-', '', full_names)
  gene_indices = dataset$gene_names[[1]] %in% full_names
  
  # Reformat data for pre-processing
  input_frame = as.data.frame(dataset$gex_fc[gene_indices,]) %>%
    mutate(gene_name = (dataset$gene_names[[1]])[gene_indices],
           gene_type = gsub('.*-', '', dataset$gene_names[[1]])[gene_indices]) %>%
    gather('Sample', 'Gex', ImmediateFC:`72hourFC`) %>%
    mutate(Time=rep(time_vec,each=sum(gene_indices)))
  
  ## Average replicate genes (multiple probes)
  input_matrix = average_replicate_genes(input_frame)
  
  # Ageing correction
  if(ageing_flag) {
    ageing_data = readRDS(here('Data/ageing_fc_timecourse.rds'))
    ageing_indices = ageing_data$gene_names %in% full_names
    ageing_frame = as.data.frame(ageing_data$gex_fc[ageing_indices,]) %>%
      mutate(gene_name = ageing_data$gene_names[ageing_indices],
             gene_type = gsub('.*-', '', ageing_data$gene_names)[ageing_indices]) %>%
      gather('Sample', 'Gex', `8hourFC`:`72hourFC`) %>%
      mutate(Time=rep(ageing_time_vec,each=sum(ageing_indices)))
    
    ## Average replicate genes
    ageing_matrix = average_replicate_genes(ageing_frame)[,-1] 
  } else {
    ageing_matrix = NULL
  }
  
  ## Correct for 24-hour cycle effects
  if(circadian_correction) {
    input_matrix = make_cyc_columns(input_matrix)
  } 
  
  input_matrix = input_matrix[,-1]
  
  # Combine matrix of ageing time course with insecticide exposure time course
  if(ageing_flag) {
    colnames(ageing_matrix) = paste0(colnames(ageing_matrix), '_ageing')
    
    # Pad with zeros and shift by one to make sure prediction happens at the same time
    # point
    ageing_matrix = rbind(matrix(0, dim(input_matrix)[1]-dim(ageing_matrix)[1]-1, 
                                 dim(ageing_matrix)[2]), 
                          ageing_matrix)
    ageing_matrix = rbind(ageing_matrix, matrix(0, 1, dim(ageing_matrix)[2]))
    input_matrix = cbind(input_matrix, ageing_matrix)
  }
  
  # Set options for EDISON network inference run
  options = defaultOptions()
  options$burnin = TRUE
  options$cp.fixed=TRUE
  options$maxTF=3
  
  # Specify target variables (target genes)
  target_vars = grep(full_names_targets[i], colnames(input_matrix))
  target = target_vars[1]
  ageing_target = target_vars[2]
  
  ## Specify predictor variables (transcription factors)
  ## select original predictors, plus ageing value for target,
  ## but exclude target itself
  pred_vars = (1:length(full_names))[-target]
  
  ## LOESS imputation to ensure equi-distant time points
  if(loess_flag) {
    new_matrix = loess_smooth_predictors(time_vec, input_matrix, pred_vars)
  } else {
    new_matrix = input_matrix
  }
  
  if(ageing_flag) {
    # Add ageing version of target gene at time t as predictor of 
    # exposed version at time t.
    pred_vars = c(pred_vars, ageing_target)
  }
  
  if(circadian_correction) {
    # Add in sine and cosine predictors for cyclical trend
    pred_vars = c(pred_vars, length(full_names)+1, length(full_names)+1)
  }
  
  ## Apply EDISON (two independent runs to check convergence)
  test_1 = EDISON.run(new_matrix, num.iter=500000, options=options, 
                      target.vars=target, pred.vars=pred_vars)
  test_2 = EDISON.run(new_matrix, num.iter=500000, options=options, 
                      target.vars=target, pred.vars=pred_vars)
  
  ## Calculate edge posterior probabilities
  edge_probs_1 = calculateEdgeProbabilities(test_1, num.preds=length(pred_vars))
  edge_probs_2 = calculateEdgeProbabilities(test_2, num.preds=length(pred_vars))
  
  marg_probs_1 = edge_probs_1$probs.segs[[1]]
  marg_probs_2 = edge_probs_2$probs.segs[[1]]
  
  tfs = rep(colnames(input_matrix)[pred_vars],2)
  
  ## Compile results
  result_tibble = tibble(marg_probs=c(marg_probs_1, marg_probs_2), 
                         target=colnames(input_matrix)[target], 
                         tfs=tfs,
                         reps=rep(1:2, each=length(pred_vars)))
  
  suffix = ''
  
  if(loess_flag) {
    suffix = paste0(suffix, 'loess_')
  } 
  if(!circadian_correction) {
    suffix = paste0(suffix, 'no_circ_')
  }
  if(ageing_flag) {
    suffix = paste0(suffix, 'age_correction_')
  }
  if(nchar(gene_choice) > 0) {
    suffix = paste0(suffix, gene_choice, '_')
  }
  
  ## Save result for this target gene
  write_csv(result_tibble, path=here(paste0('Results/bn_marg_probs_bipartite_', suffix, i, '.csv')))
  
}