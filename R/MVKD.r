#' Performs MVKD-based LR calculation on a set of speakers
#'
#' Models reference population using a multivariate kernel density and models suspect by normal density
#'
#' @param bg_data A data frame. The first column identifies speakers.
#' @param by_speaker_data A named list of 2 sub-lists: \code{suspect_data} and \code{offender_data}. Contains data in the \code{test_speakers} set, divided by speaker. Each is a named list of data frames, with speaker IDs as the names.
#' @param test_speakers (optional) A vector. Specifies the set of speakers for which LR calculations should be derived. Takes all speakers from \code{data} if unspecified.
#' @param bg_speakers (optional) A vector. Specifies the set of speakers to form the background model for LR calculations. Equal to \code{test_speakers} if unspecified.
#' @return A named list of 3 items:
#' * \code{likelihood_ratio_matrix}: A data frame. Rows and columns are named after the speaker identifiers. Each row and column represents a speaker as suspect and offender respectively, and each cell contains a single LR score.
#' * \code{cllr}: Numeric. Reports the logLR cost.
#' * \code{eer}: Numeric. Reports the equal error rate (between 0 and 1).
#' @md
#' @export

MVKD <- function(bg_data, by_speaker_data, test_speakers, bg_speakers){
  likelihood_ratio_matrix <- matrix(0, nrow = length(test_speakers),
                                    ncol = length(test_speakers),
                                    dimnames = list(as.character(test_speakers),
                                                    as.character(test_speakers)))

  # Check if test speakers are part of background
  test_in_bg <- length(setdiff(test_speakers, bg_speakers)) == 0

  background_data_by_speaker <- split_by_speaker(bg_data)
  n_by_speaker <- sapply(background_data_by_speaker, nrow)

  # Calculate speaker means and covariance matrices in advance for multivar_kernel_LR to extract
  background_mean <- colMeans(subset(bg_data, select = -speaker))
  background_stat <- list(
    bg_means = lapply(background_data_by_speaker, colMeans),
    bg_within_covar = mapply("*",
                             n_by_speaker,
                             lapply(background_data_by_speaker, sse),
                             SIMPLIFY = FALSE)
  )
  background_stat[["bg_between"]] <- mapply(function(x, y) x * tcrossprod(y - background_mean),
                                            n_by_speaker,
                                            background_stat$bg_means,
                                            SIMPLIFY = FALSE)

  # Iterate through test speaker set to perform same- and different-speaker LR comparisons
  for(sus in test_speakers){
    suspect_data <- as.matrix(by_speaker_data$suspect_data[[sus]])

    for(off in test_speakers){
      offender_data <- as.matrix(by_speaker_data$offender_data[[off]])

      # Create background data and background index
      # If test speakers are part of background, remove suspect and offender/random speaker from background
      if(test_in_bg){
        if(sus == off) removed_speaker <- sample(setdiff(bg_speakers, sus), 1)
        else removed_speaker <- off
      }

      llr <- tryCatch({
        if(!test_in_bg) multivar_kernel_LR(suspect_data, offender_data, bg_data, background_stat)
        else multivar_kernel_LR(suspect_data, offender_data,
                                       bg_data[-which(bg_data$speaker %in% c(sus, removed_speaker)), ],
                                       background_stat_remove_speaker(background_stat, c(sus, removed_speaker)))
      },
      error = function(e){
        print(paste("Error in comparison:", sus, off))
        print(e)
        NA
      })

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

background_stat_remove_speaker <- function(background_stat, speakers) sapply(names(background_stat), function(x) speaker_filter(background_stat[[x]], speakers, remove = TRUE), simplify = FALSE, USE.NAMES = TRUE)
