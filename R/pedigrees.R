#' Convert a pedigree path degree to possible path lengths
#'
#' @param degree An integer: The pedigree (path) degree
#' @param full Logical: Is the relationship full (TRUE) or half (FALSE)
#'
#' @return A DataFrame

degree_to_l = function(degree, full = T) {

  l <- data.frame(degree = integer(),
                  l1 = integer(),
                  l2 = integer())

  if (isTRUE(full)) {
    gamma = 1
  } else {
    gamma = 0
  }

  for (l1 in 0:degree) {
    l2 = degree - l1 + gamma

    l <- rbind(l, data.frame(degree = degree,
                             l1 = l1,
                             l2 = l2))

  }

  for (l2 in 0:degree) {
    l1 = degree - l2 + gamma

    l <- rbind(l, data.frame(degree = degree,
                             l1 = l1,
                             l2 = l2))
  }

  return (l)

}

#' Convert multiple path degrees to possible path lengths
#'
#' @param min_path_degree An integer: Smallest path degree
#' @param max_path_degree An integer: Largest path degree
#'
#' @return A DataFrame
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
      }
    } else {
      half = F
      degree = defineDegree(nSteps, half)
      sexPaths = defineSexPath(nSteps, half)

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


}

sexCombinations = function(sexPathLength, ordered = T) {
  if (isTRUE(ordered)) {
    sexes = c("p", "m")
    combinations = expand.grid(rep(list(sexes), sexPathLength))
    combinations = apply(combinations, 1, paste, collapse = "")
  }

  return (combinations)
}

# ---- CREATING PEDIGREES ----

swpSx = function(x, ids)
  swapSex(x, ids, verbose = F)

lin = function(deg, swp = NULL)
  swpSx(linearPed(deg), swp)

av = function(rem = 1, swp = NULL)
  swpSx(avuncularPed(removal = rem), swp)

cous = function(deg, rem = 0, swp = NULL)
  swpSx(cousinPed(deg, removal = rem), swp)

hcous = function(deg, rem = 0, swp = NULL)
  swpSx(halfCousinPed(deg, removal = rem), swp)

# Need a function for converting sexpaths to swp

identifyLeaves = function(pedigree) {
  leaves = leaves(pedigree)

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
  actual_sexPath_int = getSex(ped, ids = pathIds)

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


# Should also create a "metadata"-like data frame
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
      ped = linearPed(nGen)


      if(!(is.na(sexPath))) {
        swpIds = getSwpSexPath(ped, sexPath)
        ped = swapSex(ped, swpIds)
      }


    } else if (type == "avuncular") {
      ped = avuncularPed(rem = removal, half = half)


      if(!(is.na(sexPath))) {
        swpIds = getSwpSexPath(ped, sexPath)
        ped = swapSex(ped, swpIds)
      }


    } else {
      ped = cousinPed(deg = cousDeg, removal = removal, half = half)

      if (!(is.na(sexPath))) {
        swpIds = getSwpSexPath(ped, sexPath)
        ped = swapSex(ped, swpIds)
      }

    }

    # If the pedigree only contain one individual
    if (length(ped$ID) == 1) next

    pedigrees = append(pedigrees, list(ped), after = length(pedigrees))

  }

  return (pedigrees)

}

annotatePedigree = function(ped, ids = NULL) {

  if (is.null(ids)) {
    leaves = identifyLeaves(ped)
  }

  relationships = verbalisr::verbalise(ped, ids = leaves)


  for (i in 1:length(relationships)) {
    relationship = relationships[[i]]
    relationship_str = relationship$rel
    relationship_sexpath <- relationship$sexPath

    annotation = paste0(relationship_str, "-",relationship_sexpath)

  }

  return (annotation)

}

annotatePedigrees = function(pedlist) {

  annotation_all = c()

  for (i in 1:length(pedlist)) {
    ped = pedlist[[i]]

    annotation = annotatePedigree(ped)

    annotation_all = c(annotation_all, annotation)

  }

  return (annotation_all)
}

pedigreesMetadata = function(pedlist) {

  nVar = 2

  metadata = matrix(nrow = length(pedlist), ncol = nVar)

  for (i in 1:length(pedlist)) {
    ped = pedlist[[i]]

    verb = verbalisr::verbalise(ped, ids = identifyLeaves(ped))

    # Only support unilineal relationships as of now
    relationship = verb[[1]]$rel
    degree = verb[[1]]$degree

    metadata[i,] = c(relationship, degree)
  }

  metadata = as.data.frame(metadata)
  colnames(metadata) = c("relationship", "degree")

  return (metadata)
}

pedigreesMetadata = function(pedlist) {
  metadata = data.frame(Relationship = names(pedlist))
  metadata$degree = sapply(pedlist, pedDegree)

  Donnelly.classes = groupDonnelly(pedlist, N = 100, seed = 1234)

  metadata = merge(metadata, Donnelly.classes, by = "Relationship",
                   all.x = T, all.y = F)

  metadata



}

pedDegree = function(ped) { # Only supporting unilineal relationships as of now
  verbalisr::verbalise(ped, ids = identifyLeaves(ped))[[1]]$degree
}

pedKappa = function(ped) { # Only supporting unilieal relationships as of now
  ribd::kappaIBD(ped, ids = identifyLeave(ped))[[1]]$degree
}

pedKinship= function(ped) { # Only supporting unilineal relationships as of now
  ribd::kinship(ped, ids = identifyLeaves(ped))[[1]]$degree
}

# Donnelly equivalences ---------
# A simulation-based approach to identify potential Donnelly-equivalences

groupDonnelly <- function(pedlist, N, seed) {
  donnelly = list()

  i = 1
  for (ped in pedlist) {
    sim = ibdsim2::ibdsim(ped, N = N, seed = seed, ids = identifyLeaves(ped))
    segments = postprocess.simulation(sim, ped)

    if (i == 1) {
      donnelly <- append(donnelly, list(segments$length))
      donnelly.rels <- data.frame(Relationship = verbalisr::verbalise(ped,
                                                                      ids = identifyLeaves(ped))[[1]]$rel,
                                  class = factor(1))

    } else {
      if (list(segments$length) %in% donnelly) {
        index = which(sapply(donnelly, function(x) identical(x, segments$length)))
        donnelly.rels <- rbind(donnelly.rels, data.frame(Relationship = verbalisr::verbalise(ped,
                                                                                             ids = identifyLeaves(ped))[[1]]$rel,
                                                         class = factor(index)))
      } else {
        donnelly <- append(donnelly, list(segments$length))
        donnelly.rels <- rbind(donnelly.rels, data.frame(Relationship = verbalisr::verbalise(ped,
                                                                                             ids = identifyLeaves(ped))[[1]]$rel,
                                                         class = factor(i)))
      }
    }

    i = i + 1
  }

  return (donnelly.rels)
}
