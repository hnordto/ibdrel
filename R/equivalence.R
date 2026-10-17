#'
#' @export
anchor <- function(ped, ids = ibdrel:::specifyIndividuals(ped)) {
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

equivalenceClass <- function(ped, ids = specifyIndividuals(ped)) {

  verb <- verbalisr::verbalise(ped, ids)

  if (length(verb) > 1) {
    stop("Relationship is not simple unilineal.")
  } else {
    verb <- verb[[1]]
  }

  type = verb$type
  half = !(verb$full)
  sexpath = verb$sexPath

  path = getPath(ped, ids)

  anchors = anchor(ped, ids)
  free = path[!(path %in% anchors)]

  anchor.sex <- pedtools::getSex(ped, anchors)
  free.sex <- pedtools::getSex(ped, free)

  switch(type,
         lineal = {
           classcode = "D"
         },
         avuncular = {
           classcode = if(half) "H" else "A"
         },
         cousin = {
           classcode = if(half) "H" else "C"
         },
         sibling = {
           classcode = if(half) "H" else NA
         }
  )

  sexInt2Chr <- c("m", "f")


  if (classcode == "C") {
    if (length(unique(anchor.sex)) == 2) {
      anchor.code = "_{fm}"
    } else {
      anchor.code = paste0(sexInt2Chr[anchor.sex], collapse = "")
      anchor.code = paste0("_{",anchor.code, "}")
    }
  } else if(classcode == "H" || classcode == "A") {
    anchor.code = paste0(sexInt2Chr[anchor.sex],collapse = "")
    anchor.code = paste0("_{",anchor.code, "}")
  } else {
    anchor.code = "_{}"
  }

  n.sex.free = table(sexInt2Chr[free.sex])

  nf = n.sex.free["f"]
  nm = n.sex.free["m"]

  nf[is.na(nf)] <- 0
  nm[is.na(nm)] <- 0

  eqclass = paste0(classcode,anchor.code,"(",nf,",",nm,")")

  eqclass


}

equivalenceRep <- function(equivalence.class, ignoreSex = FALSE) {
  code.split = strsplit(equivalence.class, "_")

  classcode = code.split[[1]][1]

  anchor.sex = sub(".*\\{([^}]*)\\}.*","\\1", code.split[[1]][2])
  free.sex = strcapture( "\\((-?[0-9]+),\\s*(-?[0-9]+)\\)",
                         code.split[[1]][2],
                         proto = list(nf = integer(), nm = integer()))
  nf = free.sex$nf
  nm = free.sex$nm
  n = nf+nm
  z = nchar(anchor.sex)


  # Standard sex-aware representative: Least number of removals
  # Assign anchors: m = paternal, f = maternal
  # Assign free from left to right: rep("m", nf)+rep("p", nm)
  # See Vigeland, 2022 for computations

  sexChr2Int <- list("m" = 1, "f" = 2)

  if (classcode == "C") {
    d = ((n+1)*(z+1))/3+2

    l <- degree_to_l(d, full = TRUE)
    l <- l[!duplicated(l),]

    l <- l[(l$l1 > 1 & l$l2 > 1),]

    cousDeg <- pmin(l$l1, l$l2) - 1
    removal <- abs(l$l1-l$l2)

    cousDeg <- unique(cousDeg[which(removal == min(removal))])

    ped <- pedtools::cousinPed(degree = cousDeg, removal = min(removal), half = FALSE)



  } else if (classcode == "H") {
    d = ((n+1)*(z+1))/2+1

    l <- degree_to_l(d, full = FALSE)
    l <- l[!duplicated(l),]

    cousDeg <- pmin(l$l1, l$l2) - 1
    removal <- abs(l$l1-l$l2)

    cousDeg <- unique(cousDeg[which(removal == min(removal))])

    # Avuncular pedigrees can be constructed as cousinPeds with cousDeg = 0

    ped <- pedtools::cousinPed(degree = cousDeg, removal = min(removal), half = TRUE)

  } else if (classcode == "A") {
    d = ((n+1)*(z+1))/2+1

    l <- degree_to_l(d, full = TRUE)
    l <- l[!duplicated(l),]

    l <- l[(l$l1 > 1 & l$l2 == 1),]

    removal <- abs(l$l1-l$l2)

    ped <- pedtools::avuncularPed(removal = min(removal), half = FALSE)

  } else {
    d = (n+1)*(z+1)


    ped <- pedtools::linearPed(d)

    annot <- annotatePedigree(ped,)

  }

  path = getPath(ped, ids = specifyIndividuals(ped))

  anchors = anchor(ped, specifyIndividuals(ped))
  free = path[!(path %in% anchors)]

  anchor.sexes <- strsplit(anchor.sex, "")[[1]]

  swap.anchor <- anchors[which(anchor.sexes == "f")]
  swap.free <- free[0:nf]

  ped <- ped |> pedtools::swapSex(c(swap.anchor, swap.free), verbose = FALSE)

  annot = annotatePedigree(ped, ignoreSex = ignoreSex)

  return(annot)

}

#'
#' @export
pedVariants <- function(ped, ids = specifyIndividuals(ped), onlyExtremes = FALSE) {

  verb <- verbalisr::verbalise(ped, ids)

  type <- verb[[1]]$type
  degree = verb[[1]]$degree
  cousDeg = if (type == "cousin" || type == "sibling") min(verb[[1]]$nSteps)-1 else NA
  half = !(verb[[1]]$full)
  removal = verb[[1]]$removal
  nGen = if(type == "lineal" || type == "avuncular") degree+as.integer(isFALSE(half)) else NA


  path = getPath(ped, ids)
  anchors = anchor(ped, ids)
  free = path[!(path %in% anchors)]

  anchors.idx = which(path %in% anchors)
  free.idx = which(path %in% free)

  anchor.comb <- sexCombinations(length(anchors), FALSE)
  free.comb <- sexCombinations(length(free), FALSE)

  if (length(free.comb) == 0) {
    free.comb = "m" # Arbitrary
  }

  if (length(anchor.comb) == 0) {
    anchor.comb = "m"  # "Arbitrary
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

      if (identical(strsplit(sexPath, "")[[1]], character(0))) {
        sexPath = NA
      }


      variants[i,] = c(type, degree, removal, nGen, cousDeg, half, sexPath)
      i = i + 1
    }

  }

  colnames(variants) <- c("type", "degree", "removal", "nGen", "cousDeg", "half",
                          "sexPath")
  variants = as.data.frame(variants)

  if (isTRUE(onlyExtremes)) {
    extreme.idx = which(sapply(variants$sexPath, function(s) length(unique(strsplit(s, "")[[1]])) == 1))
    variants = variants[extreme.idx,]
  }

  return (variants)


}
