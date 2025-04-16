###########################
# SIMULATION              #
###########################

# Wrapper function for ibdsim2::ibdsim for running multiple simulations
# with the same parameters

ibdSimulations = function(pedlist,
                           N = 1000,
                           seed = NULL,
                           ids = NULL,
                           map = "decode19",
                           model = "chi") {

  if (!(length(seed) == length(pedlist) | length(seed) == 1 | is.null(seed))) {
    stop("Seed must either be of length ",length(pedlist), " or 1")
  }

  if (is.null(seed)) {
    seeds = seq(1:length(pedlist))
  }

  simulations = list()
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

  return (simulations)

}


###########################
# POSTPROCESSING          #
###########################

postprocess.simulation = function(sim,
                                  ped,
                                  ids = NULL,
                                  cutoff = 7,
                                  unit = "cm") {
  if (is.null(ids)) {
    ids = identifyLeaves(ped)
  }

  ibd_pattern = findPattern(sim,
                            pattern = list(carriers = ids),
                            cutoff = cutoff,
                            unit = "cm")

  ibd = bind_rows(lapply(ibd_pattern, data.frame), .id = "sim")

  colnames(ibd) <- c("sim", "chrom", "startMB", "endMB", "startCM",
                        "endCM", "x.P", "x.M", "y.P", "y.M", "IBD", "sigma")

  if (unit == "cm") {
    ibd$length = ibd$endCM - ibd$startCM
  } else if (unit == "mb") {
    ibd$length = ibd$endMB - ibd$startMB
  } else {
    stop("Invalid unit. Must be one of 'cm', 'mb'.")
  }

  return (ibd)
}

postprocessSimulations = function(simlist,
                                   pedlist,
                                   annotationlist,
                                   cutoff = 7,
                                   unit = "cm") {

  if(!(length(simlist) == length(pedlist))) {
    stop("Number of pedigrees and number of simulation runs must be equal.")
  }

  if(!(length(pedlist) == length(annotationlist))) {
    stop("All pedigrees need annotation.")
  }

  for (i in 1:length(simlist)) {
    sim = simlist[[i]]
    ped = pedlist[[i]]
    annotation = annotationlist[i]

    ibd_segments = postprocess.simulation(sim = sim, ped = ped,
                                          cutoff = cutoff,
                                          unit = unit)

    ibd_segments$kinship = annotation
    ibd_segments$degree = verbalisr::verbalise(ped, ids = identifyLeaves(ped))[[1]]$degree

    if (i == 1) {
      allsegments = ibd_segments
    } else {
      allsegments = rbind(allsegments, ibd_segments)
    }
  }

  return (allsegments)

}


lengthIBD = function(segments) {
  relationships = unique(segments$kinship)

  segmentsLst = list()

  for (relationship in relationships) {
    segmentLengths = segments |>
      filter(kinship == relationship) |>
      group_by(sim) |> # !
      summarise(length = list(length)) |>
      ungroup()

    for (i in 1:nrow(segmentLengths)) {
      if (i == 1) {
        segmentLst = segmentLengths$length[i]
      } else {
        segmentLst = append(segmentLst, segmentLengths$length[i])
      }
    }

    segmentsLst = append(segmentsLst, list(segmentLst))

  }

  names(segmentsLst) = relationships
  return (segmentsLst)

}
