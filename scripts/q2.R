#This script calculates the Q2 among two vectors. These vectors can be, for example, germline gene or diversity profiles.

library(forecast) #needed for CV function


germline_gene_frequencies_df <- data.frame(apply(replicate(1, runif(100, min = 0, max = 1)), 2, function(x) x/sum(x))) #mock germline gene frequency vector normalized to 1
colnames(germline_gene_frequencies_df) <- c("V1") 

germline_gene_frequencies_df[,2] <- germline_gene_frequencies_df[,1] + 0.05 * runif(length(germline_gene_frequencies_df[,1]), min = 0, max = 1 ) # generate second correlated germline gene vector

germline_gene_frequencies_df[,2] <- germline_gene_frequencies_df[,2]/sum(germline_gene_frequencies_df[,2]) #normalize second germline gene vector to 1

q2 <- function(germline_gene_frequencies_df) {
    
    germline_gene_frequencies_df_scaled <- scale(germline_gene_frequencies_df) #mean centering, scaling to unit variance
    
    fit <- lm( germline_gene_frequencies_df_scaled[,2] ~  germline_gene_frequencies_df_scaled[,1]) #linear regression model
    
    unname((1-CV(fit)[1]/mean(germline_gene_frequencies_df_scaled[,2]^2))*100) #calculation of Q2
  }


######
#Output using above mock germline gene frequency vectors
######

q2(germline_gene_frequencies_df) 
