#' Performs GMM-UBM-based LR calculation on a set of speakers
#'
#' Models reference population using a single GMM (UBM) and models suspect by adapting from UBM. Called by LR_test when the mode GMM-UBM is chosen. Procedure is based on \code{Reynolds, D.A., Quatieri, T.F. & Dunn, R.B. (2000). Speaker Verification Using Gaussian Mixture Models. Digital Signal Processing, 10(1-3), 19-41.}
#'
#' @param bg_data A data frame of background data. The first column identifies speakers. All other columns contain data.
#' @param by_speaker_data A named list of 2 sub-lists: \code{suspect_data} and \code{offender_data}. Contains data in the \code{test_speakers} set, divided by speaker. Each is a named list of data frames, with speaker IDs as the names.
#' @param test_speakers A vector of test speakers.
#' @param G Number of components in the Gaussian Mixture Model.
#' @param r Relevance factor for speaker adaptation.
#' @param background_model (Optional) Model pre-fitted on bg_data. Useful when the same model is used for different test data sets.
#' @return A named list of 3 items:
#' * \code{likelihood_ratio_matrix}: A data frame. Rows and columns are named after the speaker identifiers. Each row and column represents a speaker as suspect and offender respectively, and each cell contains a single logLR score.
#' * \code{cllr}: Numeric. Reports the logLR cost.
#' * \code{eer}: Numeric. Reports the equal error rate (between 0 and 1).
#' @md
#' @importFrom mclust Mclust
#' @importFrom mvnfast dmvn
#' @export

GMM_UBM_loop <- function(bg_data, by_speaker_data, test_speakers, G = 8, r = 16, background_model = NULL){
  likelihood_ratio_matrix <- matrix(0, nrow = length(test_speakers),
                                    ncol = length(test_speakers),
                                    dimnames = list(as.character(test_speakers),
                                                    as.character(test_speakers)))

  data_dim <- ncol(bg_data) - 1

  UBM <- background_model %||% Mclust(bg_data[, -1], G = G)

  bg_pro <- UBM$parameters$pro
  bg_mean <- t(UBM$parameters$mean)
  if(data_dim == 1) bg_mean <- t(bg_mean)
  bg_variance <- array_to_list(UBM$parameters$variance$sigma, data_dim)

  # Iterate through test-speaker set to perform same- and different-speaker comparisons
  for(sus in test_speakers){
    suspect_data <- as.matrix(by_speaker_data$suspect_data[[sus]])
    suspect_alignment <- predict(UBM, suspect_data)$z

    # Compute sufficient statistics
    suff_weight <- colSums(suspect_alignment)
    suff_mean <- (t(suspect_alignment) %*% suspect_data) / suff_weight
    suff_variance <- (t(suspect_alignment) %*% (suspect_data * suspect_data)) / suff_weight

    # Adapted parameters: proportions, means and variances
    adapt_coeff <- suff_weight / (suff_weight + r)
    sus_pro <- adapt_coeff * suff_weight / nrow(suspect_data) + (1 - adapt_coeff) * bg_pro
    sus_pro <- sus_pro / sum(sus_pro)

    sus_mean <- adapt_coeff * suff_mean + (1 - adapt_coeff) * bg_mean
    if(data_dim == 1){
      sus_variance <- mapply("+",
                             mapply("*", bg_variance, 1 - adapt_coeff, SIMPLIFY = FALSE),
                             split(suff_variance * adapt_coeff +
                                     bg_mean * bg_mean * (1 - adapt_coeff) -
                                     sus_mean * sus_mean,
                                   seq(G)),
                             SIMPLIFY = FALSE)
    } else {
      sus_variance <- mapply("+",
                             mapply("*", bg_variance, 1 - adapt_coeff, SIMPLIFY = FALSE),
                             lapply(split(suff_variance * adapt_coeff +
                                            bg_mean * bg_mean * (1 - adapt_coeff) -
                                            sus_mean * sus_mean,
                                          seq(G)),
                                    diag),
                             SIMPLIFY = FALSE)
    }

    for(off in test_speakers){
      offender_data <- as.matrix(by_speaker_data$offender_data[[off]])

      # Calculate log-likelihood ratio, using average LLR to compensate for dependencies
      llr <- mean(log(mapply(dmvn,
                             mu = matrix_to_rows(sus_mean),
                             sigma = sus_variance,
                             MoreArgs = list(X = offender_data)
                             ) %*% sus_pro) -
                    log(mapply(dmvn,
                               mu = matrix_to_rows(bg_mean),
                               sigma = bg_variance,
                               MoreArgs = list(X = offender_data)
                               ) %*% bg_pro))

      likelihood_ratio_matrix[as.character(sus), as.character(off)] <- llr
    }
  }

  LR_SS <- diag(likelihood_ratio_matrix)
  LR_DS <- non_diag(likelihood_ratio_matrix)
  likelihood_ratio_matrix <- as.data.frame(likelihood_ratio_matrix)

  c <- cllr(LR_SS, LR_DS, log = TRUE)
  eer <- EER_linear(LR_SS, LR_DS, log = TRUE)

  return(list(likelihood_ratio_matrix = likelihood_ratio_matrix, cllr = c, eer = eer))
}
