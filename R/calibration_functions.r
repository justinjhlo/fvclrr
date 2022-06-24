#' Fuses multiple LR systems and provides single calibrated output
#'
#' Combines LR scores from (multiple) test systems to provide a single set of fused LR scores, based on a logistric-regression model trained with the training data input. Self-calibration based on training with the same set of data is possible, in which case the training data argument should be omitted.
#' An R implementation of \code{calibration_fusion.m} by Hughes (2015) and \code{train_llr_fusion_robust.m} by Morrison (2009) (http://geoff-morrison.net/#TrainFus).
#'
#' @param test_list A list of LR matrices to be calibrated and fused. The input must be a LIST. Any data frame or matrix will work, as long as it follows the same format as the output in MVKD_loop (a row for each suspect, a column for each offender, LR value not transformed in any way). The LR matrices should all be of the same dimensions.
#' @param dev_list (Optional) A list of LR matrices to be trained on. The list must be of the same length as test_list and matrices within must correspond with those in test_list. If left blank, training will be based on test_list.
#' @param log Boolean. LRs already log-transformed (TRUE; default) or not (FALSE)?
#' @return A named list following the exact format as MVKD_loop: likelihood_ratio_matrix, cllr and eer.
#' @export

calibration_fusion <- function(test_list, dev_list = NULL, log = TRUE, ...){
  dev_list <- dev_list %||% test_list
  d <- length(test_list)
  speakers <- colnames(test_list[[1]])

  dev_list <- lapply(dev_list, as.matrix)
  test_list <- lapply(test_list, as.matrix)
  if(!log){
    dev_list <- lapply(dev_list, log)
    test_list <- lapply(test_list, log)
  }

  # Add matrix of ones for offset
  test_list[[d + 1]] <- test_list[[1]] * 0 + 1

  # Split dev LR matrices into SS and DS for training
  train_targets <- t(sapply(dev_list, diag))
  train_non_targets <- t(sapply(dev_list, non_diag))

  weights <- train_llr_fusion_robust(train_targets, train_non_targets, ...)

  # Fuse and calibrate test LR matrices
  fused_LR <- Reduce('+', mapply("*", weights, test_list, SIMPLIFY = FALSE))

  calibrated_SS <- diag(fused_LR)
  calibrated_DS <- non_diag(fused_LR)

  fused_LR <- setNames(as.data.frame(fused_LR, row.names = speakers), speakers)
  fused_cllr <- cllr(calibrated_SS, calibrated_DS, log = TRUE)
  fused_EER <- EER_linear(calibrated_SS, calibrated_DS, log = TRUE)

  return(list(likelihood_ratio_matrix = fused_LR, cllr = fused_cllr, eer = fused_EER))
}

train_llr_fusion_robust <- function(targets, non_targets, prior = 0.5, robust_weight = 0){
  d <- nrow(targets)
  nt <- ncol(targets)
  nn <- ncol(non_targets)

  x <- rbind(cbind(targets, non_targets, -non_targets, -targets),
             c(rep_len(1, nn + nt), rep_len(-1, nn + nt)))

  weights <- c(rep_len(prior * (nn / nt + 1), nt),
               rep_len(prior * (nt / nn + 1) * robust_weight, nn),
               rep_len((1 - prior) * (nt / nn + 1), nn),
               rep_len((1 - prior) * (nn / nt + 1) * robust_weight, nt))

  offset <- log(prior / (1 - prior)) * c(rep_len(1, nt + nn), rep_len(-1, nt + nn))

  w <- numeric(d + 1)
  gp <- w
  up <- NULL
  repeat{
    train <- update_w(w, x, weights, offset, gp, up)
    w <- w + train$term
    if(max(abs(train$term)) < 1e-5) break
    gp <- train$gp
    up <- train$up
  }

  return(w)
}

update_w <- function(w, x, weights, offset, gp, up){
  s1 <- 1 / (1 + exp(as.vector(w %*% x) + offset))
  g <- as.vector(x %*% (s1 * weights))
  if(is.null(up)) u <- g
  else u <- cg_dir(up, g, gp)

  ug <- as.numeric(u %*% g)
  ux <- (as.vector(u %*% x)) ^ 2
  a <- weights * s1 * (1 - s1)
  uhu <- as.numeric(ux %*% a)

  if(uhu < .Machine$double.xmin) uhu <- .Machine$double.xmin

  return(list(term = u * ug / uhu, gp = g, up = u))
}

cg_dir <- function(old_dir, grad, old_grad){
  delta <- grad - old_grad
  den <- as.numeric(old_dir %*% delta)
  if(den == 0) return(grad * 0)
  beta <- as.numeric(grad %*% delta) / den
  return(grad - beta * old_dir)
}
