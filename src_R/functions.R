accuracy <- function(matrix) {
  # True positive
  tp <- matrix[2, 2]
  # true negative
  tn <- matrix[1, 1]
  return((tn + tp) / sum(matrix))
}

precision <- function(matrix) {
  # True positive
  tp <- matrix[2, 2]
  # false positive
  fp <- matrix[1, 2]
  return(tp / (tp + fp))
}

recall <- function(matrix) {
  # true positive
  tp <- matrix[2, 2]
  # false positive
  fn <- matrix[2, 1]
  return(tp / (tp + fn))
}

f1 <- function(matrix) {
  prec <- precision(matrix)
  rec <- recall(matrix)
  2 * ((prec * rec) / (prec + rec))
}

"%ni%" <- Negate("%in%") # define 'not in' func
