#' Calculates likelihood ratio based on MVKD
#'
#' This function is in part an adaptation of \code{Morrison, G.S. (2007). Matlab implementation of Aitken & Lucy\'s (2004) forensic likelihood-ratio software using multivariate-kernel-density estimation. [Software]}. (http://geoff-morrison.net/#MVKD)
#' For details of how MVKD works, see \code{Aitken, C.G.G. & Lucy, D. (2004). Evaluation of trace evidence in the form of multivariate data. Applied Statistics, 54, 109-122.}
#'
#' @param sus_data A matrix containing all suspect data to be compared
#' @param off_data A matrix containing all offender data to be compared
#' @param bg A data frame that contains all data in the background model. The first column identifies the speaker.
#' @param bg_stat (Optional) Data frame of pre-calculated means and covariance matrices. Avoid providing your own bg_stat at all cost, because it will be generated within if left blank. It exists primarily for use in conjunction with MVKD_loop(), which will generate it in prior to save computation time.
#' @return Numeric. The log-likelihood ratio.
#' @export

multivar_kernel_LR <- function(sus_data, off_data, bg, bg_stat = NULL){
  # Background setup
  num_obs <- nrow(bg)
  num_var <- ncol(bg) - 1

  # Calculate speaker means and covariance matrices, if not already supplied
  if(is.null(bg_stat)){
    # Split full background data by speaker
    background_data_by_speaker <- split_by_speaker(bg_data)

    # Calculate speaker means and covariance matrices in advance for multivar_kernel_LR to extract
    background_mean <- colMeans(subset(bg_data, select = -speaker))
    bg_stat <- list(
      bg_means = lapply(background_data_by_speaker, colMeans),
      bg_within_covar = mapply(function(x, y) nrow(x) * y,
                               background_data_by_speaker,
                               lapply(background_data_by_speaker, sse),
                               SIMPLIFY = FALSE)
    )
    bg_stat[["bg_between"]] <- mapply(function(x, y) nrow(x) * tcrossprod(y - background_mean),
                                      background_data_by_speaker,
                                      bg_stat$bg_means,
                                      SIMPLIFY = FALSE)
  }

  # Number of speakers in the background
  num_speakers <- length(bg_stat[[1]])

  # Extract speaker means and covariance matrices
  bg_means <- do.call(rbind, bg_stat$bg_means)
  bg_within_covar <- Reduce('+', bg_stat$bg_within_covar) * num_speakers / (num_obs * (num_obs - num_speakers))
  bg_between <- Reduce('+', bg_stat$bg_between)
  bg_within_covar_inv <- solve(bg_within_covar)

  # Suspect & Offender
  n_susdata <- nrow(sus_data)
  sus_mean <- colMeans(sus_data)
  sus_covar_inv <- bg_within_covar_inv * n_susdata

  n_offdata <- nrow(off_data)
  off_mean <- colMeans(off_data)
  off_covar_inv <- bg_within_covar_inv * n_offdata

  sus_off_mean_diff <- off_mean - sus_mean
  sus_off_mean_typicality <- solve(off_covar_inv + sus_covar_inv,
                                   off_covar_inv %*% off_mean + sus_covar_inv %*% sus_mean)

  # Kernel
  smooth_power <- 1 / (num_var + 4)
  smooth_param <- (4 / (num_speakers * (2 * num_var + 1))) ^ smooth_power
  kernel <- smooth_param^2 * (bg_between / (num_speakers - 1) - bg_within_covar) * num_speakers / num_obs
  kernel_inv <- solve(kernel)

  off_kern_inv <- solve(bg_within_covar / n_offdata + kernel)
  sus_kern_inv <- solve(bg_within_covar / n_susdata + kernel)
  off_sus_covar_inv <- bg_within_covar_inv * n_susdata * n_offdata / (n_susdata + n_offdata)
  off_sus_kern_inv <- solve(bg_within_covar / (n_susdata + n_offdata) + kernel)

  kernden_at_typicality <- sum(apply(mean_distance(bg_means, sus_off_mean_typicality), 1, kernel_build, kern_inv = off_sus_kern_inv))
  dist_bg_to_sus <- sum(apply(mean_distance(bg_means, sus_mean), 1, kernel_build, kern_inv = sus_kern_inv))
  dist_bg_to_off <- sum(apply(mean_distance(bg_means, off_mean), 1, kernel_build, kern_inv = off_kern_inv))

  # Likelihood ratio: equation simplified to reduce computation
  llr <- 0.5 * (determinant(kernel)$modulus +
                 determinant(off_covar_inv + kernel_inv)$modulus +
                 determinant(sus_covar_inv + kernel_inv)$modulus -
                 determinant(off_covar_inv + sus_covar_inv + kernel_inv)$modulus) +
    log(kernel_build(sus_off_mean_diff, off_sus_covar_inv)) +
    log(kernden_at_typicality) +
    log(num_speakers) -
    log(dist_bg_to_off) -
    log(dist_bg_to_sus)

  return(llr)
}

kernel_build <- function(x, kern_inv) as.numeric(exp(-0.5 * x %*% (kern_inv %*% x)))
