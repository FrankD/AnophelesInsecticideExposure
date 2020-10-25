# Utility functions

# Function for doing a LOESS smoothing of predictors
# in order to obtain t-0.5 values. Variables not in pred_vars are
# left unchanged from the input matrix.
loess_smooth_predictors <- function(time, data_matrix, pred_vars, interval=0.5) {
  new_matrix = data_matrix
  
  for(pred_var in pred_vars) {
    loess_mod = loess(data_matrix[,pred_var]~time)
    # Remove -0.5 as loess cannot estimate that. Add T (maximum time) as the 
    # last data point to obtain a complete matrix (but this is not used in estimation).
    new_matrix[,pred_var] = predict(loess_mod, c(time-interval, time[length(time)]))[-1]
  }  
  
  return(new_matrix)
}

## Function to average genes with multiple probes 
average_replicate_genes <- function(input_frame) {
  # average
  gene_data = input_frame %>%
    group_by(Sample, Time, gene_name, gene_type) %>%
    summarize(Gex=mean(Gex)) %>%
    mutate(ID=gene_name)
  
  # create matrix
  preprocessed_data = gene_data %>% 
    ungroup() %>% 
    dplyr::select(ID,Gex,Time) %>% 
    spread(ID, `Gex`)
  
  input_matrix = as.matrix(preprocessed_data)
  
  return(input_matrix)
}

# Function for correcting for circadian rhythms
# Generate two additional columns to capture possible 
# cycles
make_cyc_columns <- function(X, period=1/3){
  max_time = max(X[,1])
  
  X[,1] = X[,1]/max_time
  
  x_cyc_cols = cbind(sin(2*pi*X[,1]/period),
                     cos(2*pi*X[,1]/period))
  colnames(x_cyc_cols) = c('sin_cyc', 'cos_cyc')
  
  X_new = cbind(X, x_cyc_cols)
  
  return(X_new)
}