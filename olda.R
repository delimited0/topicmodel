olda <- function(dtm, K, alpha = 0.01, eta = 0.01, kappa = .75, tau = 20, 
                 gam_tol = 1e-3, doc_max_iter = 100) {
  if (K %% 1 != 0) stop("K must be an integer")
  if (K <= 0) stop("K must be positive")
  if (kappa <= 0.5 | kappa > 1) stop("kappa must be in (0.5, 1]")
  if (tau < 0) stop("tau must be >= 0")
  if (class(dtm)[1] != "DocumentTermMatrix") stop("dtm must be of class DocumentTermMatrix")
  
  result <- lda_online(as.matrix(dtm), K, alpha, eta, kappa, tau, gam_tol, doc_max_iter)
  result$vocab <- colnames(dtm)
  class(result) <- c("lda", "online")
  return(result)
}
  