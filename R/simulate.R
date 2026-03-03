###########################
# SIMULATION              #
###########################

# Wrapper function for ibdsim2::ibdsim for running multiple simulations
# with the same parameters

#'
#'@export
ibdSimulations = function(pedlist,
                           N = 1000,
                           seed = NULL,
                           ids = NULL,
                           map = "decode19",
                           model = "chi",
                           keep = "everything",
                           cutoff = 0) {

  if (!(length(seed) == length(pedlist) | length(seed) == 1 | is.null(seed))) {
    stop("Seed must either be of length ",length(pedlist), " or 1")
  }

  if (is.null(seed)) {
    seeds = rep(NULL, length(pedlist))
  } else {
    seeds = seed
  }

  simulations = list()

  if (keep == "everything") {

    for (i in 1:length(pedlist)) {

      ped = pedlist[[i]]
      seed = seeds[i]

      if (is.null(ids)) { # Currently only supporting pedigrees with natural "leaves"
        ids = identifyLeaves(ped)
      }

      simulation = ibdsim2::ibdsim(x = ped, N = N, seed = seed,
                                   ids = ids, map = map, model = model)

      simulations = append(simulations, list(simulation), after = length(simulations))

      ids = NULL

    }
  } else if (keep == "nonzero") {

    seedlst = list()

    for (i in 1:length(pedlist)) {

      nonzero = 0

      ped = pedlist[[i]]

      if (is.null(ids)) {
        ids = identifyLeaves(ped)
      }

      simulations_tmp = list()

      while (nonzero < N) {

        trySeed = sample(1:10000000,1)

        simulation = ibdsim2::ibdsim(x = ped, N = 1, seed = trySeed,
                                     ids = ids, map = map, model = model,
                                     verbose = FALSE)
        segments = ibdsim2::findPattern(simulation, pattern = list(carriers = ids),
                                        cutoff = cutoff, unit = "cm")

        if (nrow(segments) > 0) {
          simulation = list(simulation)
          names(simulation) = trySeed
          simulations_tmp = append(simulations_tmp, simulation, after = length(simulations_tmp))
          nonzero = nonzero + 1
        }

      }

      simulations = append(simulations, list(simulations_tmp), after = length(simulations_tmp))

      ids = NULL
    }
  }



  return (simulations)

}


###########################
# POSTPROCESSING          #
###########################


# Updated method for postprocessing

extractLength <- function(mtrx) {
  unname(mtrx[,5]-mtrx[,4])
}

#'
#'@export
ibdfindr2ibdrel <- function(ibdfindr) {
  segments = ibdfindr$segments
  unname(segments[,3]-segments[,2])
}

#'
#'@export
lengthIBD = function(simlist,
                     pedlist,
                     annotationlist,
                     cutoff = 0) {

  segmentLst = list()

  for (i in 1:length(simlist)) {
    sim = simlist[[i]]
    ped = pedlist[[i]]
    annotation = annotationlist[i]

    segments = ibdsim2::findPattern(sim, pattern = list(carriers = identifyLeaves(ped)),
                                    cutoff = cutoff, unit = "cm")
    segments = lapply(segments, extractLength)

    segmentLst = append(segmentLst, list(segments))

  }

  names(segmentLst) = annotationlist

  return (segmentLst)

}

#'
#'@export
aggregateSegments = function(segments, metadata, metadata.agg.column) {

  lookup = lookupClass(NULL, metadata.agg.column, "rel", metadata)

  agg.ids = lookup[names(segments)]

  aggregatedSegments = purrr::map2(segments, agg.ids, ~ tibble::tibble(agg.id = .y, segment = list(.x))) |>
    dplyr::bind_rows() |>
    dplyr::group_by(agg.id) |>
    dplyr::summarise(segments = list(unlist(segment, recursive = FALSE)),
                     .groups = "drop") |>
    tibble::deframe()

  return (aggregatedSegments)
}

