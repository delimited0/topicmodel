lda <- function(dtm, K, alpha=0.01) {
  if (K %% 1 != 0) stop("K must be an integer")
  if (K <= 0) stop("K must be positive")
  
  if (class(dtm)[1] == "DocumentTermMatrix") stop("dtm must be of class DocumentTermMatrix")
  
}