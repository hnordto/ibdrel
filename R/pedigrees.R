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

strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="") # Get the reverse of a stirng, e.g. a sexpat

removeSymmetries <- function(rels.df) {
  rels.df$sexPathRev = strReverse(rels.df$sexPath)
  rels.df <- rels.df[!duplicated(apply(rels.df, 1, function(x)
    paste(sort(x), collapse = ''))),]
  rels.df <- subset(rels.df, select = -c(sexPathRev))
  return (rels.df)
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

# Helper functions from verbalisr

vrb = function(x, ids = leaves(x), paths = FALSE, ...) {
  format(verbalise(x, ids), includePaths = paths, ...)
}

vrbAbbr = function(x, ids = leaves(x)) {
  vrb(x, ids, abbreviate = TRUE, simplify = TRUE, cap = FALSE)
}

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


    switch(type,
           lineal = {
             annot = paste0("L-",degree)
           },
           sibling = {
             annot = if(half) "hS" else "fS"
           },
           avuncular = {
             annot = paste0(if(half) "hA-" else "fA-",degree)
           },
           cousin = {
             cousDeg = min(length(rel$v1),length(rel$v2))-1
             removal = rel$removal
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
  sexpath.distinc = unique(strsplit(sexpath, ""[[1]]))

  if (lengt(sexpath.distinct) == 1) {
    return (sexpath.distinct[1])
  } else {
    return (sexpath)
  }
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

# Donnelly equivalences ---------
# A simulation-based approach to identify potential Donnelly-equivalences

groupDonnelly.dep <- function(pedlist, N, seed) {
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

# More efficient approach using verbalisr() and skipping the simulation
# Does not work on sex-specific paths
# Only support unilieal relationships as of now
# Avuncular and cousin are equal in this setting
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
      if (isTRUE(full)) { # if full
        class.identifier.detailed = paste0("fC-",degree,"-",sexpath)
        class.identifier = paste0("fC-",degree)
      } else {
        class.identifier.detailed = paste0("hC-", degree, "-", sexpath)
        class.identifier = paste0("hC-", degree)
      }
    }

    if (type == "avuncular") {
      if (isTRUE(full)) {
        class.identifier.detailed = paste0("fA-",degree,"-",sexpath)
        class.identifier = paste0("fA-",degree)
      } else {
        class.identifier.detailed = paste0("hA-", degree, "-", sexpath)
        class.identifier = paste0("hA-", degree)
      }
    }

    if (type == "lineal") {

      if (degree == 1) {
        class.identifier.detailed = paste0("L-",degree)
        class.identifier = class.identifier.detailed
      } else {
        class.identifier.detailed = paste0("L-", degree, "-", sexpath)
        class.identifier = paste0("L-", degree)
      }
    }

    if (type == "sibling") {
      if (isTRUE(full)) {
        class.identifier.detailed = paste0("fS-", degree)
        class.identifier = class.identifier.detailed
      } else {
        class.identifier.detailed = paste0("hS-",degree,"-",sexpath)
        class.identifier = paste0("hS-",degree)
      }
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

#' Translate a relationship class code on standardized ibdrel format to readable relationships
#'
classTranslator <- function(class, resolution) {

  if (resolution == "eqclass.detailed" || resolution == "eqclass") {
    dictionary <- list("L" = "Linear",
                       "fA" = "Full avuncular",
                       "hA" = "Half avuncular",
                       "fC" = "Full cousins",
                       "hC" = "Half cousins",
                       "fS" = "Full siblings",
                       "hS" = "Half siblings",
                       "m" = "maternal",
                       "p" = "paternal")

    class.split <- strsplit(class, "-")[[1]]
    reltype = class.split[1]
    degree = class.split[2]

    if (!is.na(class.split[3])) {
      sexpath = class.split[3]

      sexpath.split = strsplit(sexpath, "")[[1]]
      sexpath.expanded = dictionary[sexpath.split]
      sexpath.write = paste0("(",paste(sexpath.expanded, collapse = "-"),")")
    } else {
      sexpath.write = ""
    }

    rel.write <- paste(dictionary[reltype], "of degree", degree, sexpath.write)
  } else if (resolution == "kappa") {
    class.split <- as.numeric(strsplit(class, "-")[[1]])
    kappa0 = as.character(MASS::fractions(class.split[1]))
    kappa1 = as.character(MASS::fractions(class.split[2]))
    kappa2 = as.character(MASS::fractions(class.split[3]))

    rel.write <- paste("\u03BA\u2080 =", kappa0, ",",
                       "\u03BA\u2081 =", kappa1, ",",
                       "\u03BA\u2082 =", kappa2)
  } else if (resolution == "kinship") {
    rel.write <- paste("\u03C6 =", as.character(MASS::fractions(as.numeric(class))))
  } else if (resolution == "degree") {
    rel.write <- paste("Degree", class)
  } else {
    warning("Unknown classification resolution.")
    rel.write <- NULL
  }




  return(rel.write)
}
