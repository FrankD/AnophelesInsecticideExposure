# This codes runs the de-sparsified node-wise lasso implementation in the R package SILGGM
# on fold change transcriptomic data from 31 microarrays of resistant Anopheles gambiae s.l populations
# to build a gaussian graphical model, or gene regulatory network.
# 3 datasets are remove due to very weak phenotype (Hai, Dar, Muheza)

library(tidyverse)
library(SILGGM)

fold_changes <- read.delim("Fold Changes.txt", 
                           header=TRUE, 
                           sep="\t", 
                           as.is=TRUE, 
                           row.names = NULL)

############# Take the mean of duplicate probes ############
# deselect metadata
metadata <- fold_changes[1:4,]
fold_changes <- fold_changes[5:14915,]
fold_changes$Detoxification.Classification<- NULL
fold_changes$TranscriptType <- NULL
fold_changes <- fold_changes[,1:32] #select expression values only, not q values

#make numeric 
for (i in 2:32){
  fold_changes[[i]] <- as.numeric(fold_changes[[i]])
}

#merge duplicate probes, take average
final_fc <- fold_changes %>% 
  group_by(SystematicName) %>% 
  summarise_all(funs(mean))

#gene names to rownames
final_fc <- final_fc %>% column_to_rownames(var='SystematicName')

#select only resistant populations
final_resistant <- final_fc %>% select(-'HaiFC', -'DarFC', -'MuhezaFC')

#run SILGGM - desparsified lasso
genes.SILGGM <- SILGGM(as.matrix(t(final_resistant)), 
                       cytoscape_format = TRUE)

#save to csv file
write.csv(genes.SILGGM[[1]], "all_transcripts.dslasso.csv")

#filter to 0.1 p value to reduce file size
genes.SILGGM <- genes.SILGGM %>% filter("p_precision" <= 0.1)


#read gene names in
names1 <- read.csv("names1.csv")
names <- read.csv("names.csv")

#left_join to keep all rows of X
genes.SILGGM2 <- left_join(genes.SILGGM, names1, by='gene1')
genes.SILGGM <- left_join(genes.SILGGM2, names, by='gene2')

#write to file
write.csv(genes.SILGGM, "all_transcripts.dslasso.0.1.csv", 
          row.names = FALSE)



