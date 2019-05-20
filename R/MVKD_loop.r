#' Performs LR calculation on a set of speakers
#'
#' Iterates through all pairings in a set of test speakers and calculates the likelihood ratio for each speaker-pair
#'
#' @param data A data frame. The first column identifies speakers.
#' @param test_speakers (optional) A vector. Specifies the set of speakers for which LR calculations should be derived. Takes all speakers from \code{data} if unspecified.
#' @param bg_speakers (optional) A vector. Specifies the set of speakers to form the background model for LR calculations. Equal to \code{test_speakers} if unspecified.
#' @param data_col (optional) A vector. Specifies which columns in \code{data} to take into account in LR calculation. Uses all available columns if unspecified (excludes first column, which identifies speakers).
#' @param test_data (optional) A data frame in the same format as \code{data}. If specified, offender/disputed-speaker data will be taken from here. Suspect/known-speaker data and background data will be taken from \code{data}.
#' @param test_data_col (optional) A vector. Specifies which columns in \code{data} to take into account in LR calculation. Follows \code{data_col} if unspecified.
#' @param cross_full (optional) Boolean. When \code{test_data} and \code{test_data_col} are specified, determines if the full set of data from each speaker are used in the comparison. Default TRUE. If FALSE, suspect data come from first half of \code{data} and questioned data come from second half of \code{test_data}, to simulate the scenario when \code{test_data} is not specified. Useful for keeping number of data points constant, in the case of comparisons involving non-contemporaneous recordings or mismatched conditions and the results of which are compared with contemporaneous comparisons.
#' @return A named list of 3 items:
#' * \code{likelihood_ratio_matrix}: A data frame. Rows and columns are named after the speaker identifiers. Each row and column represents a speaker as suspect and offender respectively, and each cell contains a single LR score.
#' * \code{cllr}: Numeric. Reports the logLR cost.
#' * \code{eer}: Numeric. Reports the equal error rate (between 0 and 1).
#' @md
#' @importFrom dplyr group_by mutate select summarise summarise_all transmute
#' @importFrom purrr map map_int map2
#' @importFrom tidyr nest unnest
#' @export

MVKD_loop <- function(data, test_speakers = NULL, bg_speakers = NULL, data_col = NULL, test_data = NULL, test_data_col = NULL, cross_full = TRUE){
  ## Data preparation

  # Ignore test_data_col if test_data is not specified
  if(is.null(test_data) & !is.null(test_data_col)) warning("Warning: test_data not specified; data_col used instead of test_data_col")

  # Discard unused data columns if specified
  if(!is.null(data_col)) data <- data[, c(1, data_col)]

  # Rename speaker column to ensure uniformity
  colnames(data)[1] <- "speaker"

  # Remove rows with NA or NaN in data
  data <- data[complete.cases(data), ]

  # Do the same to test_data if specified
  if(!is.null(test_data)){
    # If test_data_col is not specified, follow data_col
    if(is.null(test_data_col)) test_data_col <- data_col

    if(!is.null(test_data_col)) test_data <- test_data[, c(1, test_data_col)]
    colnames(test_data)[1] <- "speaker"
    test_data <- test_data[complete.cases(test_data), ]
  }

  # Use all speakers as test_speakers if not specified
  if(is.null(test_speakers)) test_speakers <- as.vector(unique(data$speaker))

  # Use test speakers as background if not specified
  if(is.null(bg_speakers)) bg_speakers <- test_speakers

  # End program if not all speakers contained in (test_)data
  if(length(setdiff(test_speakers, data$speaker)) > 0) stop("Error: test_speakers not all contained in data")
  if(length(setdiff(bg_speakers, data$speaker)) > 0) stop("Error: bg_speakers not all contained in data")
  if(!is.null(test_data)){
    if(length(setdiff(test_speakers, test_data$speaker)) > 0) stop("Error: test_speakers not all contained in test_data")
    if(length(setdiff(bg_speakers, test_data$speaker)) > 0) stop("Error: bg_speakers not all contained in test_data")
  }

  ## Likelihood ratios

  # Count number of test speakers
  num_test_speakers <- length(test_speakers)

  # Set up LR matrix, SS matrix and DS matrix
  likelihood_ratio_matrix <- data.frame(matrix(0, nrow = num_test_speakers, ncol = num_test_speakers))
  colnames(likelihood_ratio_matrix) <- as.character(test_speakers)
  rownames(likelihood_ratio_matrix) <- as.character(test_speakers)
  LR_SS <- numeric()
  LR_DS <- numeric()

  # Check if test speakers are part of background
  test_in_bg <- length(setdiff(test_speakers, bg_speakers)) == 0

  # Create full background data
  background_data_full_nested <- nest(group_by(data, speaker))
  background_data_full_nested <- background_data_full_nested[which(background_data_full_nested$speaker %in% bg_speakers), ]

  # Calculate speaker means and covariance matrices in advance for multivar_kernel_LR to extract
  background_mean <- unlist(summarise_all(select(unnest(background_data_full_nested), -speaker), mean))
  background_stat <- select(mutate(background_data_full_nested,
                                   n = map_int(data, nrow),
                                   sse = map(data, sse),
                                   bg_means = map(data, colMeans),
                                   bg_within_covar = map2(n, sse, ~ .x * .y),
                                   bg_between = map2(n, bg_means, ~ .x * tcrossprod(.y - background_mean))),
                            speaker, bg_means, bg_within_covar, bg_between)

  # If test speakers are not part of background, set up background for LR comparison
  if(!test_in_bg){
    background_data <- unnest(background_data_full_nested)
    bg_stat <- background_stat[, -1]
  }

  # Retrieve suspect and offender data for each test speaker from data (and test_data, if applicable)
  # If test_data not specified: suspect data come from first half of data, offender data come from second half of data
  if(is.null(test_data)){
    by_speaker_data <- select(mutate(nest(group_by(data, speaker)),
                                     ntoken = map_int(data, nrow),
                                     suspect_data = map2(data, ntoken, ~ .x[1:ceiling(.y/2), ]),
                                     offender_data = map2(data, ntoken, ~ .x[(ceiling(.y/2) + 1):.y, ])),
                              speaker, suspect_data, offender_data)
    by_speaker_data <- by_speaker_data[which(by_speaker_data$speaker %in% test_speakers), ]
  } else {
    # If test_data specified and cross_full is TRUE: suspect data come from data, offender data come from test_data
    if(cross_full){
      by_speaker_suspect_data <- nest(group_by(data, speaker), .key = "suspect_data")
      by_speaker_suspect_data <- by_speaker_suspect_data[which(by_speaker_suspect_data$speaker %in% test_speakers), ]
      by_speaker_offender_data <- nest(group_by(test_data, speaker), .key = "offender_data")
      by_speaker_offender_data <- by_speaker_offender_data[which(by_speaker_offender_data$speaker %in% test_speakers), ]
      by_speaker_data <- merge(by_speaker_suspect_data, by_speaker_offender_data)
    } else {
      # If test_data specified and cross_full is FALSE: suspect data come from first half of data, offender data come from second half of test_data
      by_speaker_suspect_data <- select(mutate(nest(group_by(data, speaker)),
                                               ntoken = map_int(data, nrow),
                                               suspect_data = map2(data, ntoken, ~ .x[1:ceiling(.y/2), ])),
                                        speaker, suspect_data)
      by_speaker_suspect_data <- by_speaker_suspect_data[which(by_speaker_suspect_data$speaker %in% test_speakers), ]
      by_speaker_offender_data <- select(mutate(nest(group_by(test_data, speaker)),
                                                ntoken = map_int(data, nrow),
                                                offender_data = map2(data, ntoken, ~ .x[(ceiling(.y/2) + 1):.y, ])),
                                         speaker, offender_data)
      by_speaker_offender_data <- by_speaker_offender_data[which(by_speaker_offender_data$speaker %in% test_speakers), ]
      by_speaker_data <- merge(by_speaker_suspect_data, by_speaker_offender_data)
    }
  }

  # Iterate through test speaker set to perform same- and different-speaker LR comparisons
  for(sus in test_speakers){
    # Assign suspect data
    suspect_data <- as.matrix(by_speaker_data$suspect_data[which(by_speaker_data$speaker == sus)][[1]])

    for(off in test_speakers){
      # Assign offender data
      offender_data <- as.matrix(by_speaker_data$offender_data[which(by_speaker_data$speaker == off)][[1]])

      # Create background data and background index
      # If test speakers are part of background, remove suspect and offender/random speaker from background
      if(test_in_bg){
        if(sus == off){
          # Generate random speaker to be removed
          removed_speaker <- sample(setdiff(bg_speakers, sus), 1)
        } else removed_speaker <- off

        background_data <- unnest(background_data_full_nested[-which(background_data_full_nested$speaker %in% c(sus, removed_speaker)), ])
        bg_stat <- background_stat[-which(background_stat$speaker %in% c(sus, removed_speaker)), -1]
      }

      # Run LR on reduced/full background set
      lr <- multivar_kernel_LR(suspect_data, offender_data, background_data, bg_stat)

      # Store LR in LR matrix
      likelihood_ratio_matrix[as.character(sus), as.character(off)] <- lr

      # If same-speaker comparison, store LR in SS
      if(sus == off) LR_SS[length(LR_SS) + 1] <- lr
      # If different-speaker comparison, store LR in DS
      else LR_DS[length(LR_DS) + 1] <- lr
    }
  }

  # Cllr
  c <- cllr(LR_SS, LR_DS, log = TRUE)

  # EER
  eer <- EER_linear(LR_SS, LR_DS)

  # Output: LR matrix, Cllr and EER
  return(list(likelihood_ratio_matrix = likelihood_ratio_matrix, cllr = c, eer = eer))
}

# MVKD_loop <- function(data, test_speakers = NULL, bg_speakers = NULL, data_col = NULL){
#   # Discard unused data columns if specified
#   if(!is.null(data_col)) data <- data[, c(1, data_col)]
#
#   # Rename speaker column to ensure uniformity
#   colnames(data)[1] <- "speaker"
#
#   # Use all speakers as test speakers if not specified
#   if(is.null(test_speakers)) test_speakers <- as.vector(unique(data$speaker))
#
#   # Use test speakers as background if not specified
#   if(is.null(bg_speakers)) bg_speakers <- test_speakers
#
#   # Count number of test speakers
#   num_test_speakers <- length(test_speakers)
#
#
#   # Set up LR matrix, SS matrix and DS matrix
#   likelihood_ratio_matrix <- data.frame(matrix(0, nrow = num_test_speakers, ncol = num_test_speakers))
#   colnames(likelihood_ratio_matrix) <- as.character(test_speakers)
#   rownames(likelihood_ratio_matrix) <- as.character(test_speakers)
#   LR_SS <- numeric()
#   LR_DS <- numeric()
#
#   # Check if test speakers are part of background
#   test_in_bg <- length(setdiff(test_speakers, bg_speakers)) == 0
#
#   # Create full background data and index
#   background_data_full <- data[which(data$speaker %in% bg_speakers), ]
#
#   # Calculate speaker means and covariance matrices in advance, so each multivar_kernel_LR only has to extract rows
#   background_mean <- unlist(summarise_all(select(background_data_full, -speaker), mean))
#   background_stat <- select(mutate(tidyr::nest(group_by(background_data_full, speaker)),
#                                    n = map_int(data, nrow),
#                                    sse = map(data, sse),
#                                    bg_means = map(data, colMeans),
#                                    bg_within_covar = map2(n, sse, ~ .x * .y),
#                                    bg_between = map2(n, bg_means, ~ .x * tcrossprod(.y - background_mean))),
#                             speaker, bg_means, bg_within_covar, bg_between)
#
#   # If test speakers are not part of background, set up background for LR comparison
#   if(!test_in_bg){
#     background_data <- background_data_full
#     bg_stat <- background_stat
#   }
#
#   # Retrieve range of row indices for each speaker in data
#   by_speaker_data <- mutate(tidyr::nest(group_by(data, speaker)),
#                             ntoken = map_int(data, nrow),
#                             firsthalf = map2(data, ntoken, ~ .x[1:ceiling(.y/2),]),
#                             secondhalf = map2(data, ntoken, ~ .x[(ceiling(.y/2) + 1):.y,]))
#
#   # If test speakers are in background, retrieve where range of row indices for each speaker in background data frame
#   if(test_in_bg) within_bg_index <- summarise(group_by(tibble::rowid_to_column(background_data_full, "id")[1:2], speaker), first = min(id), last = max(id))
#
#   # Iterate through test speaker set to perform same- and different-speaker LR comparison
#   for(sus in test_speakers){
#     # Assign suspect data: first half of data (rounded up)
#     suspect_data <- as.matrix(by_speaker_data$firsthalf[which(by_speaker_data$speaker == sus)][[1]])
#
#     # If suspect is part of background, select from background for removal later
#     if(test_in_bg) isus_in_bg <- within_bg_index[which(within_bg_index$speaker == sus), ]
#
#     for(off in test_speakers){
#       # Assign offender data: second half of data (rounded down)
#       offender_data <- as.matrix(by_speaker_data$secondhalf[which(by_speaker_data$speaker == off)][[1]])
#
#       # Create background data and background index and run LR
#       # If test speakers are part of background, remove suspect and another speaker from background
#       if(test_in_bg){
#         if(sus == off){
#           # Generate random speaker to be removed
#           rand_speaker <- sample(setdiff(bg_speakers, sus), 1)
#           removed_speaker <- within_bg_index[which(within_bg_index$speaker == rand_speaker), ]
#         } else removed_speaker <- within_bg_index[which(within_bg_index$speaker == off), ]
#
#         # Remove suspect and random speaker/offender from background for same/different-speaker comparison
#         background_data <- background_data_full[-c(isus_in_bg$first:isus_in_bg$last, removed_speaker$first:removed_speaker$last), ]
#         bg_stat <- background_stat[-which(background_stat$speaker %in% c(sus, removed_speaker$speaker)), -1]
#       }
#
#       # Run LR on reduced/full background set
#       lr <- multivar_kernel_LR(suspect_data, offender_data, background_data, bg_stat)
#
#       # Store LR in LR matrix
#       likelihood_ratio_matrix[as.character(sus), as.character(off)] <- lr
#
#       # If same-speaker comparison, store LR in SS
#       if(sus == off) LR_SS[length(LR_SS) + 1] <- lr
#       # If different-speaker comparison, store LR in DS
#       else LR_DS[length(LR_DS) + 1] <- lr
#     }
#   }
#
#   # Cllr
#   c <- cllr(LR_SS, LR_DS, log = FALSE)
#
#   # EER
#   eer <- EER_linear(LR_SS, LR_DS)
#   #eer <- EER(LR_SS, LR_DS)
#
#   # Output: LR matrix, Cllr and EER
#   return(list(likelihood_ratio_matrix = likelihood_ratio_matrix, cllr = c, eer = eer))
# }
