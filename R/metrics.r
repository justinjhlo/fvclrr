#' Calculates the log likelihood ratio cost
#' 
#' @param LLR_SS A vector of LRs from same-speaker comparisons.
#' @param LLR_DS A vector of LRs from different-speaker comparisons.
#' @param log (optional) Boolean. Specifies whether LRs supplied are already in (natural) log form (default TRUE).
#' @return A single numerical Cllr.
#' @md
#' @export
cllr <- function(LLR_SS, LLR_DS, log = TRUE){
  # More numerically robust implementation to work better with log LRs (using Maechler, 2012)
  # See: Maechler, M. (2012). Accurately Computing log(1-exp(-|a|)) Assessed by the Rmpfr package.
  
  # (Updated 2025-03-18) generates warning if contains NA
  if(anyNA(LLR_SS)) warning("cllr: LLR_SS contains NA values")
  if(anyNA(LLR_DS)) warning("cllr: LLR_DS contains NA values")
  
  # (Updated 2025-03-18) removes any NA if present
  LLR_SS <- LLR_SS[!is.na(LLR_SS)]
  LLR_DS <- LLR_DS[!is.na(LLR_DS)]

  if(!log){
    LLR_SS <- log(LLR_SS)
    LLR_DS <- log(LLR_DS)
  }

  cllr_SS <- mean(log1pexp(-LLR_SS)) / log(2)
  cllr_DS <- mean(log1pexp(LLR_DS)) / log(2)

  # Cllr
  c <- (cllr_SS + cllr_DS) / 2

  return(c)
}

#' Calculates the equal error rate
#' 
#' @param LLR_SS A vector of LRs from same-speaker comparisons.
#' @param LLR_DS A vector of LRs from different-speaker comparisons.
#' @param log (optional) Boolean. Specifies whether LRs supplied are already in (natural) log form (default FALSE).
#' @return A single numerical EER between 0 and 1 (inclusive).
#' @md
#' @export
EER_linear <- function(LR_SS, LR_DS, log = FALSE){
  # (Updated 2025-03-18) generates warning if contains NA
  if(anyNA(LR_SS)) warning("EER: LLR_SS contains NA values")
  if(anyNA(LR_DS)) warning("EER: LLR_DS contains NA values")
  
  # (Updated 2025-03-18) removes any NA if present
  LR_SS <- LR_SS[!is.na(LR_SS)]
  LR_DS <- LR_DS[!is.na(LR_DS)]
  
  # Convert to log LR if not already logged
  if(!log){
    LR_SS <- log(LR_SS)
    LR_DS <- log(LR_DS)
  }

  num_thresholds <- 5000
  counter <- 0
  min_threshold <- min(LR_SS, LR_DS)
  max_threshold <- max(LR_SS, LR_DS)

  if(max_threshold == Inf) max_threshold <- log(.Machine$double.xmax)
  if(min_threshold == -Inf) min_threshold <- log(.Machine$double.xmin)

  SS_corr <- numeric(num_thresholds)
  DS_corr <- SS_corr

  for(threshold in seq(from = min_threshold, to = max_threshold, length.out = num_thresholds)){
    counter <- counter + 1
    SS_corr[counter] <- sum(LR_SS > threshold)
    DS_corr[counter] <- sum(LR_DS <= threshold)
  }

  SS_corr <- SS_corr / length(LR_SS)
  DS_corr <- DS_corr / length(LR_DS)

  diffind <- which(abs(SS_corr - DS_corr) == min(abs(SS_corr - DS_corr)))

  EER <- 1 - min(SS_corr[diffind] + DS_corr[diffind]) / 2

  return(EER)
}

EER_alt <- function(LR_SS, LR_DS){
  SS_total <- length(LR_SS)
  DS_total <- length(LR_DS)

  LR_all <- as.data.frame(matrix(c(LR_SS, LR_DS, rep(c(TRUE, FALSE), c(SS_total, DS_total))), ncol = 2))
  LR_all <- LR_all[order(LR_all[, 1], LR_all[, 2]), ]

  LR_all$SS <- cumsum(LR_all[, 2]) / SS_total
  LR_all$DS <- 1 - cumsum(!LR_all[, 2]) / DS_total
  LR_all$ER <- (LR_all$SS + LR_all$DS) / 2

  EER <- min(LR_all$ER)
}
