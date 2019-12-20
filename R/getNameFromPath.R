#' Get name from path
#'
#' Obtains sample names from paths.
#' @param path Character or character vector with the paths to convert to sample names.
#' @param suffix Suffix to remove from the file name.
#' @param prefix Prefix to remove from the file name.
#' @export
getNameFromPath <- function(path,
                            suffix="",
                            prefix="") {
  names <- gsub(prefix, "", gsub(suffix, "", basename(path)))
  return(names)
}
