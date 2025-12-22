#'
#' @export
anchor <- function(ped, ids = ibdrel:::identifyLeaves(ped)) {
  # Only support relationships connected through one ancestral path


  verb <- verbalisr::verbalise(ped, ids)
  if (length(verb) > 1) {
    stop("Relationship is connected through more than one ancestral path.")
  }

  comAnc <- pedtools::commonAncestors(ped, ids)

  if (length(comAnc) == 1) {
    return (comAnc)
  }

  chr <- pedtools::children(ped, comAnc)

  anch = chr[!(chr %in% ids)]

  return (anch)

}

#'
#' @export
pedVariants <- function(ped, ids = identifyLeaves(ped)) {

  verb <- verbalisr::verbalise(ped, ids)

  type <- verb[[1]]$type
  degree = verb[[1]]$degree
  cousDeg = if (type == "cousin") min(verb[[1]]$nSteps)-1 else NA
  half = !(verb[[1]]$full)
  removal = verb[[1]]$removal
  nGen = if(type == "lineal" || type == "avuncular") removal-1 else NA


  path = ibdrel:::getPath(ped, ids)
  anchors = anchor(ped, ids)
  free = path[!(path %in% anchors)]

  anchors.idx = which(path %in% anchors)
  free.idx = which(path %in% free)

  anchor.comb <- sexCombinations(length(anchors), FALSE)
  free.comb <- sexCombinations(length(free), FALSE)

  if (length(free.comb) == 0) {
    free.comb = "m" # Arbitrary
  }

  nvariants = length(anchor.comb)*length(free.comb)




  variants = matrix(NA, nrow = nvariants, ncol = 7)

  i = 1
  for (anchor.sexes in anchor.comb) {
    anchor.sexes.split = strsplit(anchor.sexes, "")[[1]]


    for (free.sexes in free.comb) {
      free.sexes.split = strsplit(free.sexes, "")[[1]]

      sexPath = character(length = length(path))
      sexPath[anchors.idx] = anchor.sexes.split
      sexPath[free.idx] = free.sexes.split
      sexPath = paste(sexPath, collapse = "")

      variants[i,] = c(type, degree, removal, nGen, cousDeg, half, sexPath)
      i = i + 1
    }

  }

  colnames(variants) <- c("type", "degree", "removal", "nGen", "cousDeg", "half",
                          "sexPath")
  variants = as.data.frame(variants)
  return (variants)


}
