#' Calculate SS and DS mean scores from likelihood ratio matrix
#'
#' Not ready for use
#'
#' @param LR.matrix Data frame of likelihood ratio scores. Row and column names must correspond to the speaker names/identifiers. LR scores should not be transformed in any way (e.g. logarithmically).
#' @return Data frame where each row contains the SS and DS mean scores for a speaker.
#' @importFrom dplyr arrange group_by mutate right_join summarise
#' @importFrom tidyr gather
#' @export

zooify <- function(LR.matrix){
  # Convert raw LR matrix to data frame where each row is a log10LR entry
  LRlog <- arrange(gather(tibble::rownames_to_column(log10(LR.matrix), "suspect"), "offender", "LLR", -suspect), suspect)
  SS <- LRlog[LRlog$suspect == LRlog$offender, ]
  DS <- LRlog[LRlog$suspect != LRlog$offender, ]

  LR.zoo <- summarise(group_by(mutate(SS, speaker = as.character(suspect)), speaker), SS = mean(LLR))
  LR.zoo <- right_join(summarise(group_by(gather(DS, "identity", "speaker", -LLR), speaker), DS = mean(LLR)), LR.zoo, by = "speaker")

  return(zoo = LR.zoo)
}

#' Zooplot
#'
#' @param LR.zoo Zoo data frame returned by \code{zooify}.
#' @param quantile Plot line segments marking quantiles of scores if TRUE.
#' @import ggplot2 ggrepel
#' @export

zooplot <- function(LR.zoo, quantile = TRUE){
  zoop <- ggplot(LR.zoo, aes(SS, DS)) + geom_point() + scale_y_reverse() + geom_text_repel(aes(label = speaker))
  if(quantile){
    SS.low <- quantile(LR.zoo$SS, 0.25)
    SS.high <- quantile(LR.zoo$SS, 0.75)
    DS.low <- quantile(LR.zoo$DS, 0.25)
    DS.high <- quantile(LR.zoo$DS, 0.75)
    zoop <- zoop + geom_segment(aes(x = SS.low, y = DS.low, xend = -Inf, yend = DS.low)) +
      geom_segment(aes(x = SS.high, y = DS.low, xend = Inf, yend = DS.low)) +
      geom_segment(aes(x = SS.low, y = DS.high, xend = -Inf, yend = DS.high)) +
      geom_segment(aes(x = SS.high, y = DS.high, xend = Inf, yend = DS.high)) +
      geom_segment(aes(x = SS.low, y = DS.high, xend = SS.low, yend = Inf)) +
      geom_segment(aes(x = SS.low, y = DS.low, xend = SS.low, yend = -Inf)) +
      geom_segment(aes(x = SS.high, y = DS.high, xend = SS.high, yend = Inf)) +
      geom_segment(aes(x = SS.high, y = DS.low, xend = SS.high, yend = -Inf))
  }
  return(zoop)
}
