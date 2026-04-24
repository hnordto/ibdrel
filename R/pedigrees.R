#'
#' @export
degree_to_l = function(degree, full = T) {


  l <- data.frame(degree = integer(),
                  l1 = integer(),
                  l2 = integer())

  if (isTRUE(full)) {
    gamma = 1
  } else {
    gamma = 0
  }

  # If l1 or l2 = 0, gamma = 0 (special case for lineal rels)



  for (l1 in 0:degree) {

    if (l1 == 0) {
      l2 = degree - l1
    } else {
      l2 = degree - l1 + gamma
    }


    l <- rbind(l, data.frame(degree = degree,
                             l1 = l1,
                             l2 = l2))

  }

  for (l2 in 0:degree) {

    if (l2 == 0) {
      l1 = degree - l2
    } else {
      l1 = degree - l2 + gamma
    }


    l <- rbind(l, data.frame(degree = degree,
                             l1 = l1,
                             l2 = l2))
  }



  return (l)

}

#'
#' @export
identify_path_lengths = function(min_path_degree = 0,
                                 max_path_degree) {
  l <- data.frame(degree = integer(),
                  l1 = integer(),
                  l2 = integer())

  for (degree in min_path_degree:max_path_degree) {
    l <- rbind(l, degree_to_l(degree))
  }

  return (l)
}



#'
#' @export
listRelationships = function(l_df,
                              ignoreSex = F,
                              ignoreSymmetries = F,
                              full = T) {
  rels = data.frame(type = character(),
                    degree = integer(),
                    removal = integer(),
                    nGen = integer(),
                    cousDeg = integer(),
                    half = logical(),
                    sexPath = character())

  for (i in 1:nrow(l_df)) {
    l1 = l_df[i, "l1"]
    l2 = l_df[i, "l2"]

    if (l1 == 0 || l2 == 0) {
      type = "lineal"
    } else if (l1 == 1 || l2 == 1) {
      type = "avuncular"
    } else {
      type = "cousin"
    }

    nSteps = c(l1, l2)


    removal = abs(diff(nSteps)) # Removal


    # Number of generations (for lineal and avuncular)
    if (type == "lineal" || type == "avuncular") {
      if (removal > 0) {
        nGen = removal - 1
      } else {
        nGen = 0
      }
    } else {
      nGen = NA
    }

    # Cousin degree
    if (type == "cousin") {
      cousDeg = min(l1, l2) - 1
    } else {
      cousDeg = NA
    }

    if(!(type == "lineal")) {
      for (half in c(FALSE, TRUE)) {
        degree = defineDegree(nSteps, half)

        if (degree > max(l_df$degree)) next # Skip if degree exceeds the max listed (case for half relationships)

        sexPaths = defineSexPath(nSteps, half)

        if(isFALSE(ignoreSex)) {
          for (sexPath in sexPaths) {
            rels_tmp = data.frame(type = type,
                                  degree = degree,
                                  removal = removal,
                                  nGen = nGen,
                                  cousDeg = cousDeg,
                                  half = half,
                                  sexPath = sexPath)
            rels = rbind(rels, rels_tmp)
          }
        } else {
          rels_tmp = data.frame(type = type,
                                degree = degree,
                                removal = removal,
                                nGen = nGen,
                                cousDeg = cousDeg,
                                half = half,
                                sexPath = NA)
          rels = rbind(rels, rels_tmp)
        }

      }
    } else {
      half = NA
      degree = defineDegree(nSteps, half)
      sexPaths = defineSexPath(nSteps, half)

      if(isFALSE(ignoreSex)) {
        for (sexPath in sexPaths) {
          rels_tmp = data.frame(type = type,
                                degree = degree,
                                removal = removal,
                                nGen = nGen,
                                cousDeg = cousDeg,
                                half = NA,
                                sexPath = sexPath)
          rels = rbind(rels, rels_tmp)
        }
      } else {
        rels_tmp = data.frame(type = type,
                              degree = degree,
                              removal = removal,
                              nGen = nGen,
                              cousDeg = cousDeg,
                              half = NA,
                              sexPath = NA)
        rels = rbind(rels, rels_tmp)
      }




    }

  }

  # Remove duplicates

  rels = rels[!duplicated(rels),]


  return (rels)

}

defineDegree = function(nSteps, half) {
  return (sum(nSteps) - as.integer(isFALSE(half)))
}

defineSexPath = function(nSteps, half) {
  degree = defineDegree(nSteps, half)

  if (degree == 0) {
    sexPathLength = 0
  } else {
    sexPathLength = degree - 1
  }

  sexPaths = sexCombinations(sexPathLength, ordered = T)

  if (!(identical(sexPaths, character(0)))) {
    sexPaths = sexPaths
  } else {
    sexPaths = NA
  }

  return (sexPaths)

}

sexCombinations = function(sexPathLength, ordered = T) {
  sexes = c("p", "m")
  combinations = expand.grid(rep(list(sexes), sexPathLength))
  combinations = apply(combinations, 1, paste, collapse = "")
  if (isFALSE(ordered)) {
    combinations.ordered <- sapply(combinations, function(x) paste0(sort(strsplit(x, "")[[1]]), collapse = ""))
    combinations = unique(combinations.ordered)
  }

  return (combinations)
}

strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="") # Get the reverse of a stirng, e.g. a sexpat


#'
#' @export
removeSymmetries <- function(rels.df) {
  rels.df$sexPathRev = strReverse(rels.df$sexPath)
  rels.df <- rels.df[!duplicated(apply(rels.df, 1, function(x)
    paste(sort(x), collapse = ''))),]
  rels.df <- subset(rels.df, select = -c(sexPathRev))
  return (rels.df)
}


# ---- CREATING PEDIGREES ----

#' @export
identifyLeaves = function(pedigree) {
  leaves = pedtools::leaves(pedigree)

  # "Leaves" are not defined for lineal relationships
  if (length(leaves) < 2) {
    leaves = c(pedigree$ID[1], leaves)
  }

  return (leaves)
}


# The path depends on whether the relationship is full or half
getPath = function(pedigree, leaves = NULL) {

  if (is.null(leaves)) {
    leaves = identifyLeaves(pedigree)
  }

  rel = verbalisr::verbalise(pedigree, ids = leaves)

  # Only supporting unilineal relationships as of now

  if (!(length(rel) == 1)) {
    stop("Complex relationship.")
  }

  for (i in 1:length(rel)) {
    path = rel[[i]]$path
    sexPath = rel[[i]]$sexPath
    full = rel[[i]]$full

    # Split ID to path (keep ancestor if relationship is half)
    path = gsub("\\[([0-9]+,[0-9]+)\\]", "", path)
    path = gsub("\\[([0-9]+)\\]", "-\\1", path)
    path = as.numeric(unlist(strsplit(path, "-")))
    path = path[!is.na(path)]

    # Check if start and end correspond to "leaves"

    if (path[1] %in% leaves) {
      path = path[-1]
    } else {
      stop()
    }

    if (path[length(path)] %in% leaves) {
      path = path[-length(path)]
    }
  }

  return (path)

}

sexstr_to_int = function(sexstr) {
  sexes = list("p" = 1, "m" = 2)

  return (sexes[sexstr][[1]])

}

getSwpSexPath = function(ped, sexPath) {
  pathIds = getPath(ped)
  actual_sexPath_int = pedtools::getSex(ped, ids = pathIds)

  assigned_sexPath_split = strsplit(sexPath, "")[[1]]

  swpIds = c()

  for (i in 1:length(pathIds)) {
    assigned_sex = assigned_sexPath_split[[i]]
    actual_sex = actual_sexPath_int[i]
    id = pathIds[i]

    assigned_sex_int = sexstr_to_int(assigned_sex)

    if(!(actual_sex == assigned_sex_int)) {
      swpIds = c(swpIds, id)
    }
  }

  return (swpIds)

}


#' @export
constructPedigrees = function(pedigrees_df) {

  pedigrees = list()

  for (i in 1:nrow(pedigrees_df)) {

    nGen = as.integer(pedigrees_df[i, "nGen"])
    removal = as.integer(pedigrees_df[i, "removal"])
    cousDeg = as.integer(pedigrees_df[i, "cousDeg"])
    sexPath = as.character(pedigrees_df[i, "sexPath"])
    type = as.character(pedigrees_df[i, "type"])
    half = as.logical(pedigrees_df[i, "half"])

    if (type == "lineal") {
      ped = pedtools::linearPed(removal)


      if(!(is.na(sexPath))) {
        swpIds = getSwpSexPath(ped, sexPath)
        ped = pedtools::swapSex(ped, swpIds)
      }


    } else if (type == "avuncular") {
      ped = pedtools::avuncularPed(rem = removal, half = half)


      if(!(is.na(sexPath))) {
        swpIds = getSwpSexPath(ped, sexPath)
        ped = pedtools::swapSex(ped, swpIds)
      }


    } else {
      ped = pedtools::cousinPed(deg = cousDeg, removal = removal, half = half)

      if (!(is.na(sexPath))) {
        swpIds = getSwpSexPath(ped, sexPath)
        ped = pedtools::swapSex(ped, swpIds)
      }

    }

    # If the pedigree only contain one individual
    if (length(ped$ID) == 1) next

    pedigrees = append(pedigrees, list(ped), after = length(pedigrees))

  }

  return (pedigrees)

}

#' @export
annotatePedigree = function(ped, ids = NULL) {
  if (is.null(ids)) {
    leaves = identifyLeaves(ped)
  }

  verb <- verbalisr::verbalise(ped, ids = leaves)

  # Iterat over unilineal relationships
  for (i in 1:length(verb)) {
    rel <- verb[[i]]
    type = rel$type
    degree = rel$degree
    half = isFALSE(rel$full)
    removal = rel$removal
    ng = if (removal > 1) removal - 1 else 0


    switch(type,
           lineal = {
             annot = paste0("L",degree)
           },
           sibling = {
             annot = if(half) "hS" else "S"
           },
           avuncular = {
             annot = paste0(if(half) "h", if(ng > 0) paste0("G(",ng,")"), "A")
           },
           cousin = {
             cousDeg = min(length(rel$v1),length(rel$v2))-1
             annot = paste0(if(half) "h",cousDeg,"C",removal,"R")
           }
    )

    sexpath = rel$sexPath

    sexpath.annot <- condenseSexpath(sexpath)

    if (sexpath == "") {
      annotation = paste0(annot)
    } else {
      annotation = paste0(annot, "-", sexpath.annot)
    }

  }

  return (annotation)

}

condenseSexpath <- function(sexpath) {
  sexpath.distinct = unique(strsplit(sexpath, "")[[1]])

  if (length(sexpath.distinct) == 1) {
    return (sexpath.distinct[1])
  } else {
    return (sexpath)
  }
}

#' @export
pedsMetadata = function(pedlist) {

  if (is.null(names(pedlist))) {
    names(pedlist) = sapply(pedlist, annotatePedigree)
  }

  metadata = data.frame(rel = names(pedlist))
  metadata$code = sapply(pedlist, pedCode)
  metadata$details = lapply(pedlist, pedDetails)
  metadata$degree = sapply(pedlist, pedDegree)

  metadata$kappa0 = sapply(pedlist, pedKappa, 1)
  metadata$kappa1 = sapply(pedlist, pedKappa, 2)
  metadata$kappa2 = sapply(pedlist, pedKappa, 3)
  metadata$kappa = paste0(metadata$kappa0, "-", metadata$kappa1, "-",
                          metadata$kappa2)

  metadata$kinship = sapply(pedlist, pedKinship)

  donnelly.classes = groupDonnelly(pedlist, names(pedlist))

  metadata = merge(metadata, donnelly.classes,
                   by = "rel", all.x = TRUE, all.y = FALSE)

  return(metadata)
}

pedCode = function(ped) {
  verbalisr::verbalise(ped, ids = identifyLeaves(ped))[[1]]$code
}

pedDetails = function(ped) {
  verbalisr::verbalise(ped, ids = identifyLeaves(ped))[[1]]$rel
}

pedDegree = function(ped) { # Only supporting unilineal relationships as of now
  verbalisr::verbalise(ped, ids = identifyLeaves(ped))[[1]]$degree
}

pedKappa = function(ped, which.kappa) { # Only supporting unilineal relationships as of now
  ribd::kappaIBD(ped, ids = identifyLeaves(ped))[which.kappa]

}

pedKinship= function(ped) { # Only supporting unilineal relationships as of now
  ribd::kinship(ped, ids = identifyLeaves(ped))
}


donnellyRep <- function(degree, half) {

  # Special case: Degree 1 (and 0) have no cousin relationships. These are sibling
  if (degree %in% c(0,0)) {
    return (NA)
  }

  full <- isFALSE(half)

  l <- degree_to_l(degree, full = full)
  l <- l[!duplicated(l),]

  cousDeg <- pmin(l$l1, l$l2)-1
  removal <- abs(l$l1-l$l2)


  cousDeg.rep <- cousDeg[which(removal == min(removal))[1]] # First with lowest removal

  return (list(cousDeg = cousDeg.rep, removal = min(removal)))
}

donnellyAnnot <- function(degree, half) {
  rep <- donnellyRep(degree, half)
  ped <- pedtools::cousinPed(degree = rep$cousDeg, removal = rep$removal, half = half)
  annot <- annotatePedigree(ped)
  annot <- strsplit(annot, "-")[[1]][1] # Remove sexpath, to be appended later
  return (annot)
}

# More efficient approach using verbalisr() and skipping the simulation
# Does not work on sex-specific paths
# Only support unilieal relationships as of now
# Avuncular and cousin are equal in this setting
#' @export
groupDonnelly = function(pedlist, annotation) {
  i = 1
  for (ped in pedlist) {
    verb = verbalisr::verbalise(ped, ids = identifyLeaves(ped))
    rel = annotation[i]
    type = verb[[1]]$type
    full = verb[[1]]$full
    code = verb[[1]]$code
    degree = verb[[1]]$degree
    l1 = verb[[1]]$v1
    l2 = verb[[1]]$v2
    l.l1 = length(l1)
    l.l2 = length(l2)
    sexpath = condenseSexpath(verb[[1]]$sexPath)
    nSteps = sum(verb[[1]]$nSteps)

    if (type == "cousin") {

      if(isTRUE(full)) {
        annot = donnellyAnnot(degree, half = F)

        class.identifier.detailed = paste0(annot,"-",sexpath)
        class.identifier = annot
      } else {
        annot = donnellyAnnot(degree, half = T)

        class.identifier.detailed = paste0(annot,"-",sexpath)
        class.identifier = annot
      }

    } else if (type == "avuncular")  {

      # Half avuncular grouped with half cousins
      # Full avuncular are separate classes

      if(isTRUE(full)) {
        class.identifier.detailed = rel
        class.identifier = strsplit(rel, "-")[[1]][1]
      } else {
        annot = donnellyAnnot(degree, half = TRUE)

        class.identifier.detailed = paste0(annot, "-", sexpath)
        class.identifier = annot
      }


    } else {
      class.identifier.detailed = rel
      class.identifier = strsplit(rel, "-")[[1]][1]
    }


    if (i == 1) {
      grouper = data.frame(eqclass.detailed = class.identifier.detailed,
                           eqclass = class.identifier,
                           rel = rel)
    } else {
      grouper = rbind(grouper, data.frame(eqclass.detailed = class.identifier.detailed,
                                          eqclass = class.identifier,
                                          rel = rel))
    }

    i = i + 1

  }

  return (grouper)

}

numtimes = function(n) {
  if(n < 0) stop2("`n` must be nonnegative")
  if(n == 0) return("")
  if(n == 1) return("once")
  if(n == 2) return("twice")
  paste(n, "times")
}

is_number <- function(x) grepl("^[0-9]+$", x)


#' Translate a relationship class code on standardized ibdrel format to readable relationships
#'
#' @export
classTranslator <- function(class, resolution) {

  dictionary <- list("m" = "maternal",
                     "p" = "paternal")

  if (resolution == "eqclass.detailed" || resolution == "eqclass") {
    # Split on sex path

    class.split.rel <- strsplit(class, "-")[[1]][1]
    class.split.sex <- strsplit(class, "-")[[1]][2]

    # Split on half
    rel.half <- strsplit(class.split.rel, "h")[[1]]

    if (length(rel.half) > 1) {
      rel <- rel.half[2]
      half.write = "half "
    } else {
      rel <- rel.half[1]
      half.write = ""
    }

    rel.split <- strsplit(rel, "")[[1]]

    # Check rel type. If first is numeric (valid), reltype is cousin
    if (is_number(rel.split[1])) {
      cousDeg <- scales::ordinal(as.numeric(rel.split[1]))
      rem <- rel.split[3]
      #rem <- numtimes(as.numeric(rel.split[3]))

      if (rem == 0) {
        rel.write <- paste0(half.write, cousDeg, " cous")

      } else {
        rel.write <- paste0(half.write, cousDeg, " cous ", rem, "R")
      }


    } else if (rel.split[1] == "G" || rel.split[1] == "A") {

      if (grepl("\\(\\d+\\)", rel)) {
        rem <- gsub(".*\\((\\d+)\\).*", "\\1", rel)
      } else {
        rem <- NA
      }

      if (!is.na(rem)) {
        rel.write <- paste0(half.write, " G*", rem, " avunc")
      } else {
        rel.write <- paste(half.write, "avunc")
      }

    } else if (rel.split[1] == "S") {

      rel.write <- paste(half.write, "sib")

    } else if (rel.split[1] == "L") {
      linDeg <- as.integer(strsplit(rel, "L")[[1]][2])

      if (linDeg > 2) {
        rel.write <- paste0("G*", linDeg-1, " parent")
      } else if (linDeg == 2) {
        rel.write <- paste0("G parent")
      } else {
        rel.write <- paste0("parent")
      }
    }

    if (!(is.na(class.split.sex))) {
      sexpath.split = strsplit(class.split.sex, "")[[1]]
      sexpath.expanded = dictionary[sexpath.split]
      sexpath.write = paste0("(",paste(sexpath.split, collapse = "-"),")")
    } else {
      sexpath.write = ""
    }

    rel.write = paste(rel.write, sexpath.write)

  } else if (resolution == "kappa") {
    kappa.split <- as.numeric(strsplit(class, "-")[[1]])
    kappa0 = as.character(MASS::fractions(kappa.split[1]))
    kappa1 = as.character(MASS::fractions(kappa.split[2]))
    kappa2 = as.character(MASS::fractions(kappa.split[3]))

    rel.write <- paste0("(", kappa0, " , ", kappa1, " , ", kappa2, ")")

  } else if (resolution == "kinship") {
    rel.write <- as.character(MASS::fractions(as.numeric(class)))
  } else if (resolution == "degree") {
    rel.write <- class
  } else {
    stop("Unknown classification resolution")
  }

  rel.write <- trimws(rel.write)
  return(rel.write)
}

#'
#'@export
lookupClass = function(class, to, from, metadata) {
  lookup = metadata |>
    dplyr::select(from, !!rlang::sym(to)) |>
    tibble::deframe()

  if (!is.null(class)) {
    return (lookup[class])
  } else {
    return (lookup)
  }

}
