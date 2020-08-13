#' Performs LR calculation on a set of speakers
#'
#' Iterates through all pairings in a set of test speakers and calculates the likelihood ratio for each speaker-pair. \code{GMM-UBM} and \code{MVKD}-based approaches available for use.
#'
#' @param data A data frame. The first column identifies speakers.
#' @param test_speakers (optional) A vector. Specifies the set of speakers for which LR calculations should be derived. Takes all speakers from \code{data} if unspecified.
#' @param bg_speakers (optional) A vector. Specifies the set of speakers to form the background model for LR calculations. Equal to \code{test_speakers} if unspecified.
#' @param data_col (optional) A vector. Specifies which columns in \code{data} to take into account in LR calculation. Uses all available columns if unspecified (excludes first column, which identifies speakers).
#' @param test_data (optional) A data frame in the same format as \code{data}. If specified, offender/disputed-speaker data will be taken from here. Suspect/known-speaker data and background data will be taken from \code{data}.
#' @param test_data_col (optional) A vector. Specifies which columns in \code{data} to take into account in LR calculation. Follows \code{data_col} if unspecified.
#' @param cross_full (optional) Boolean. When \code{test_data} and \code{test_data_col} are specified, determines if the full set of data from each speaker are used in the comparison. Default TRUE. If FALSE, suspect data come from first half of \code{data} and questioned data come from second half of \code{test_data}, to simulate the scenario when \code{test_data} is not specified. Useful for keeping number of data points constant, in the case of comparisons involving non-contemporaneous recordings or mismatched conditions and the results of which are compared with contemporaneous comparisons.
#' @param mode "gmm_ubm" (default) or "mvkd".
#' @param ... Additional arguments, e.g. \code{G} (default 8) to specify how many components to use in GMM-UBM and \code{r} (default 16) to specify relevance factor for the speaker-adaptation step in GMM-UBM.
#' @return A named list of 3 items:
#' * \code{likelihood_ratio_matrix}: A data frame. Rows and columns are named after the speaker identifiers. Each row and column represents a speaker as suspect and offender respectively, and each cell contains a single LR score.
#' * \code{cllr}: Numeric. Reports the logLR cost.
#' * \code{eer}: Numeric. Reports the equal error rate (between 0 and 1).
#' @md
#' @export
LR_test <- function(data, test_speakers = NULL, bg_speakers = NULL, data_col = NULL, test_data = NULL, test_data_col = NULL, cross_full = TRUE, mode = c("gmm_ubm", "mvkd"), ...){
  mode <- match.arg(mode)

  colnames(data)[1] <- "speaker"
  if(!is.null(test_data)) colnames(test_data)[1] <- "speaker"

  test_speakers <- unique(test_speakers) %||% unique(data[["speaker"]])
  bg_speakers <- unique(bg_speakers) %||% test_speakers

  # Error handling: not all speakers contained in (test_)data
  stopifnot(all(test_speakers %in% data[["speaker"]]),
            all(bg_speakers %in% data[["speaker"]]))
  if(!is.null(test_data)){
    stopifnot(all(test_speakers %in% test_data[["speaker"]]),
              all(bg_speakers %in% test_data[["speaker"]]))
  }

  ## Data preparation

  # Ignore test_data_col if test_data is not specified
  if(is.null(test_data) & !is.null(test_data_col)) warning("Warning: test_data not specified; data_col used instead of test_data_col")

  if(!is.null(data_col)) data <- data[, c(1, data_col)]

  data <- data[complete.cases(data), ]

  if(!is.null(test_data)){
    test_data_col <- test_data_col %||% data_col
    if(!is.null(test_data_col)) test_data <- test_data[, c(1, test_data_col)]
    test_data <- test_data[complete.cases(test_data), ]
  }

  background_data_full <- data[which(data$speaker %in% bg_speakers), ]

  # Retrieve suspect and offender data for each test speaker from data (and test_data, if applicable)
  # If test_data not specified: suspect data come from first half of data, offender data come from second half of data
  if(is.null(test_data)){
    by_speaker_data <- speaker_filter(split_by_speaker(data), test_speakers)
    by_speaker_data <- list(
      suspect_data = lapply(by_speaker_data, extract_data_half, half = "first"),
      offender_data = lapply(by_speaker_data, extract_data_half, half = "second")
      )
  } else {
    # (Updated 2020-08-10) If test_data specified and cross_full is TRUE: suspect data come from test_data, offender data come from data, background is derived from data
    if(cross_full){
      by_speaker_data <- list(
        suspect_data = speaker_filter(split_by_speaker(test_data), test_speakers),
        offender_data = speaker_filter(split_by_speaker(data), test_speakers)
        )
    } else {
      # (Updated 2020-08-10) If test_data specified and cross_full is FALSE: suspect data come from first half of test_data, offender data come from second half of data, background is derived from data
      by_speaker_data <- list(
        suspect_data = speaker_filter(lapply(split_by_speaker(test_data),
                                             extract_data_half, half = "first"),
                                      test_speakers),
        offender_data = speaker_filter(lapply(split_by_speaker(data),
                                              extract_data_half, half = "second"),
                                       test_speakers)
        )
      }
  }

  if(mode == "gmm_ubm") return(GMM_UBM_loop(background_data_full, by_speaker_data, test_speakers, ...))
  else return(MVKD(background_data_full, by_speaker_data, test_speakers, bg_speakers))
}
