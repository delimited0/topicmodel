get_topics <- function(lda, K) {
  term_score <- term_score(lda, K)
  vocab_mat <- matrix(data=lda$vocab, nrow=length(lda$vocab), K)
  apply(term_score, 2, function(x) {
    vocab_mat[order(x)[1:K]]
  })
}