`%||%` <- function(x, y) if(!is.null(x)) x else y

array_to_list <- function(data_array, data_dim) lapply(split(data_array, slice.index(data_array, ifelse(data_dim == 1, 1, 3))), function(x) matrix(x, nrow = data_dim))

extract_data_half <- function(data_frame, half = c("first", "second")){
  if(half == "first") return(data_frame[1:ceiling(nrow(data_frame)/2), ])
  else return(data_frame[-(1:ceiling(nrow(data_frame)/2)), ])
}

log1pexp_individual <- function(x){
  if(x <= -37) return(exp(x))
  else if(x <= 18) return(log1p(exp(x)))
  else if(x <= 33.3) return(x + exp(-x))
  else return(x)
}

log1pexp <- function(x) sapply(x, log1pexp_individual)

matrix_to_rows <- function(data_matrix) split(data_matrix, slice.index(data_matrix, 1))

mean_distance <- function(source_means, speaker_mean) sweep(-source_means, 2, speaker_mean, "+")

non_diag <- function(m) t(m)[upper.tri(m) | lower.tri(m)]

speaker_filter <- function(data_list, speakers, remove = FALSE){
  if(!remove) return(data_list[names(data_list) %in% speakers])
  else return(data_list[!(names(data_list) %in% speakers)])
}

split_by_speaker <- function(data) split(subset(data, select = -speaker), data$speaker)

# Sum of standard error function
sse <- function(x){
  xm <- colMeans(x)
  xc <- as.matrix(sweep(x, 2, xm, '-'))
  return(crossprod(xc))
}
