##
## Helper functions taken from Yun's package ruok
##


#' Replace the content of vector with different content.
#'
#' @param x a vector.
#' @param replace_vec a vector. The dictionary for replacement.
#'
#' @return a vector after replacement.
#' @export
#' @author Yun Yan
#' @examples
#' replace_vector(c('a', 'b', 'b'), c('a'='A', 'b'=100))
#' ## [1] "A"   "100" "100"
replace_vector <- function(x, replace_vec){
  # replace_vector(c('a', 'b', 'b'), c('a'='A', 'b'=100))
  y <- x
  x_from <- names(replace_vec)
  for (xx in x_from){
    y[x == xx] <- replace_vec[xx]
  }
  return(y)
}


#' Enframe a list
#'
#' @param l
#' @param name
#' @param value
#' @export
#' @author Yun Yan
enframe_list <- function(l, name='name', value='value'){
  names_l <- names(l)
  values_l <- unlist(l)
  if (all(is.null(names_l))){names_l <- 1:length(l)}
  names_l <- rep(names_l, times=sapply(l, length))
  res <- data.frame(name=names_l, value=values_l, row.names=NULL)
  colnames(res) <- c(name, value)
  res
}


#' Deframe to a list
#'
#' @param df 
#'
#' @return
#' @export
#' @author Yun Yan
#' @examples
#' deframe_to_list(iris[, c(5, 1)])
deframe_to_list <- function(df){
  stopifnot(ncol(df)==2)
  df <- as.data.frame(df)
  lv <- gtools::mixedsort(unique(df[, 1]))
  res <- lapply(lv, function(lvx){
    return(df[df[,1]==lvx, 2])
  })
  names(res) <- lv
  return(res)
}


#' Get the mode of a vector of numbers
#'
#' @param v numeric. A vector of numbers.
#'
#' @return numeric. The mode of the input numbers
#' @export
#' @author Yun Yan
#' @examples
#' math_get_mode(c(5, 10, 30, 20, 20))
math_get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



#' Tapply on a matrix
#'
#' @param mat Each row is the sample and each column is a feature.
#' @param INDEX The clusters/groups each row of the mat belongs to.
#' @param FUN Function to run apply, e.g., mean
#' @param ... 
#'
#' @return matrix
#' @export
#' @import gtools
#' @author Yun Yan
#' @examples
#' X <- as.matrix(iris[, 1:4])
#' Y <- iris$Species
#' mat_tapply(X, Y, mean, rm.na=T)
mat_tapply <- function(mat, INDEX, FUN, ...){
  # mat: a row as samples
  # INDEX: length == nrow(mat)
  n <- nrow(mat)
  p <- ncol(mat)
  stopifnot(n==length(INDEX))
  levs <- gtools::mixedsort(unique(INDEX))
  
  res <- sapply(levs, function(l){
    idx <- INDEX %in% l
    mat_l <- mat[idx, , drop=F]
    apply(mat_l, 2, FUN=FUN, ...)
  })
  return(t(res))
  
}


#' Pretty correlation result
#'
#' @param x Numeric value
#' @param y Numeric value
#' @param to_str Default: False.
#' @param ... Other parameters to pass to `cor.test`.
#'
#' @return correlation result
#' @author Yun Yan
#' @export
#'
#' @examples
#' pretty_cor_test(iris$Sepal.Length, iris$Sepal.Width)
pretty_cor_test <- function(x, y, to_str=T, ...){
  obj <- cor.test(x, y, ...)
  res <- c(obj$estimate, `pval`=obj$p.value)
  if (!to_str) {return(res)}
  res <- prettyNum(res, digits=3)
  paste(sprintf('%s=%s', names(res), as.character(res)), 
        collapse = ' ')
}