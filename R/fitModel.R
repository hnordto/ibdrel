#' Fit a relationship classification model
#'
#' Train a non-parametric Naive Bayes classification model on IBD segment lengths.
#'
#' A probability density function (PDF) for each combination of class and feature
#' is estimated by kernel density estimation (KDE). Segments below the specified
#' threshold (`cutoff`) are treated as no IBD.
#'
#' @param data IBD segment data as a list with the following elements:
#'    * `segments`: A nested list where each outer element is a pedigree relationship and each
#'    inner element is a vector of IBD segment lengths in cM for one instance of sharing. Typically returned by [ibdSimulations()] and [lengthIBD()].
#'    * `peds`: A list of the pedigrees in `segments` on `pedtools::ped` format.
#'
#'
#' @param data Nested list of training data.
#'   A hierarchical list structured as:
#'   \code{data[[class_level]][[class]] -> list of training instances},
#'   where:
#'   \itemize{
#'     \item \strong{class_level}: grouping of relationship classes
#'     \item \strong{class}: a specific relationship
#'     \item each element: a training instance (feature vector)
#'   }
#'   Allowed class levels and classes depend on your dataset and should match the expected hierarchy.
#'
#' @param features Character vector specifying which features to use.
#'   Allowed values: \code{"count"}, \code{"total"}, \code{"median"}, \code{"longest"}, \code{"shortest"}.
#'
#' @param cutoff Numeric value indicating the minimum segment length to be considered IBD.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \strong{pdfs}: nested list of PDFs, following the structure of the input \code{data}
#'     \item \strong{cutoff}: the cutoff used for IBD segment filtering
#'     \item \strong{features}: the features used for the model
#'   }
#'
#' @examples
#' ### Example 1: Fit default model ###
#' model <- fitModel()
#' model
#'
#' @export
fitModel <- function(data = NULL,
                     features = "default",
                     cutoff = 7) {

  if (is.null(data)) {
    data = loadData()
  } else {
    data = loadData(data)
  }

  metadata = data$metadata
  peds = data$peds
  segments = data$segments
  names(peds) = sapply(peds, annotatePedigree)


  if (all(features == "default")) {
    features = DEFAULT.FEATURES

    pdfs = lapply(segments, preparePdfs, featureSel = features, cutoff = cutoff)


  } else {
    if (any(!(features %in% VALID.FEATURES))) {
      stop("Invalid feature(s). Must be one or multiple of ",
           paste(VALID.FEATURES, collapse = ", "))
    }

    pdfs = lapply(segments, preparePdfs, featureSel = features, cutoff = cutoff)
  }

  names(pdfs) <- names(segments)

  model = list(
    pdfs = pdfs,
    cutoff = cutoff,
    features = features,
    peds = peds,
    segments = segments,
    metadata = metadata
  )

  class(model) <- "ibdrelModel"

  return (model)

}

#'
#' @export
print.ibdrelModel <- function(x, ...) {
  cat("<ibdrelModel>\n")
  cat("Number of probability density functions (sex-aware pedigree relationships):", length(x$pdfs), "\n")
  cat("Features:", x$features, " (default: ", DEFAULT.FEATURES, ")", "\n")
  cat("Cutoff (cM):", x$cutoff, "\n")
}
