per_term_topic_prob <- function(lambda) {
  return(lambda / rowSums(lambda))
}

per_doc_topic_prop <- function(gamma) {
  return (gamma / rowSums(gamma))
}

term_score <- function(beta) {
  beta * log(beta / apply(beta, 1, prod)^(1/nrow(beta)))          
}

get_topics <- function(lda, k) {
# show k topic words for all topics, in order of decreasing term-score (Blei and Lafferty 2009)
  if (class(lda)[1] != "lda") stop("lda must be of class lda")
  lambda <- lda$lambda
  betas <- per_term_topic_prob(lambda)
  term_scores <- term_score(betas)
  vocab_mat <- matrix(data=lda$vocab, nrow=length(lda$vocab), ncol=nrow(lda$lambda))
  apply(term_scores, 2, function(x) {
    vocab_mat[order(x, decreasing=TRUE)[1:k]]
  })
}