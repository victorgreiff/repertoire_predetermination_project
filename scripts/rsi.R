#This script calculates the repertoire similarity within a given antibody CDR3 repertoire. It can easily be expanded to calculate the similarity among any two repertoires.

library(stringdist) #necessary for R function stringdistmatrix()

cdr3s <- c("AAA", "AAB", "AAC") #mock antibody CDR3 sequences

rsi <- function(cdr3s) {
    
    nchar_cdr3s <- nchar(cdr3s) #calculate length of CDR3s
    
    normalized_distance_matrix <- stringdistmatrix(cdr3s, cdr3s, method = "lv")/ #computes Levenshtein distance among CDR3s (dissimilarity matrix)
    sapply(nchar_cdr3s, function(k) ifelse(k > nchar_cdr3s, k, nchar_cdr3s)) #computes length-based normalization so that similarity ranges between 0 and 100%

   median(1- normalized_distance_matrix[upper.tri( normalized_distance_matrix)])*100 # (i) convert dissimilarity to similarity matrix by substracting from 1, (ii) multiply by 100 so that RSI ranges between 0 and 100.  
    
} 

######
#Output using above mock antibody CDR3 sequences
######
rsi(cdr3s) #66.66667
