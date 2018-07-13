
fillTheBlanks <- function(x, missing=""){
  rle <- rle(as.character(x))
  empty <- which(rle$value == missing)
  rle$values[empty] <- rle$value[empty-1] 
  return(inverse.rle(rle))
}

normRat <- function(dat, mult = 1){
  med <- median(dat[dat != 0])
  dat[dat != 20] <- (dat[dat != 20]/med)*mult
  return(dat)
}

meanAbsDev <- function(dat){
  return(mean(abs(dat - mean(dat))))
}

normalizeRatios <- function(dat){
  dat <- dat %>% dplyr::group_by(Sample) %>%
    dplyr::mutate(norm_ratio = normRat(ratio)) %>%
    dplyr::ungroup()
  return(dat)
}

getMedianRatios <- function(dat){
  medians <- dat %>% dplyr::group_by(Sample) %>%
    dplyr::summarise(median = median(ratio),
                     norm_median = median(norm_ratio)) %>%
    dplyr::ungroup()
  return(medians)
}

# uses standard inter quartile range method to remove outliers from numeric vector
removeOutliers <- function(vec){
  quartiles <- quantile(vec, c(0.25, 0.75))
  iqr = quartiles['75%'] - quartiles['25%']
  lowerBound <- quartiles['25%'] - (iqr * 1.5)
  upperBound <- quartiles['75%'] + (iqr * 1.5)
  
  return(vec[vec >= lowerBound & vec <= upperBound])
}


