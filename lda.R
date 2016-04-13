lda <- function(dtm, K, alpha = 0.1, eta = 0.01, gam_tol = 1e-3, em_tol = 1e-3, 
                em_max_iter = 50, doc_max_iter = 100, method = "vb", burnin = 0, 
                iter = 1000) {
  if (K %% 1 != 0) stop("K must be an integer")
  if (K <= 0) stop("K must be positive")
  
  if (class(dtm)[1] != "DocumentTermMatrix") stop("dtm must be of class DocumentTermMatrix")
  
  if (method == "vb") {
    result <- lda_vb(as.matrix(dtm), K, alpha, eta, gam_tol, em_tol, em_max_iter, doc_max_iter)
    result$log_liks <- result$log_liks[result$log_liks != 0]
  }
  else if (method == "gibbs") {
    result <- lda_gibbs(as.matrix(dtm), K, alpha, eta, burnin, iter)
  }
  else {
    stop("method must be one of vb or gibbs")
  }
  result$vocab <- colnames(dtm)
  return(result)
}