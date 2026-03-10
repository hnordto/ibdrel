# Requires ibdrel::ibdrel_unilineal

library(ibdrel)

peds <- ibdrel::ibdrel_unilineal$peds

annotation <- sapply(peds, annotatePedigree)

 # Change seed
seeds <- 1000:(length(peds)-1+1000)

sims_0 <- ibdSimulations(peds, N = 100, seed = NULL, keep = "nonzero", cutoff = 0)
sims_2 <- ibdSimulations(peds, N = 100, seed = NULL, keep = "nonzero", cutoff = 2)
sims_5 <- ibdSimulations(peds, N = 100, seed = NULL, keep = "nonzero", cutoff = 5)
sims_7 <- ibdSimulations(peds, N = 100, seed = NULL, keep = "nonzero", cutoff = 7)
sims_9 <- ibdSimulations(peds, N = 100, seed = NULL, keep = "nonzero", cutoff = 9)
sims_12 <- ibdSimulations(peds, N = 100, seed = NULL, keep = "nonzero", cutoff = 12)
sims_15 <- ibdSimulations(peds, N = 100, seed = NULL, keep = "nonzero", cutoff = 15)

# Note: Seeds contained i sims/segment objects

segments_0 <- lengthIBD(sims_0, peds, names(peds))
segments_2 <- lengthIBD(sims_2, peds, names(peds))
segments_5 <- lengthIBD(sims_5, peds, names(peds))
segments_7 <- lengthIBD(sims_7, peds, names(peds))
segments_9 <- lengthIBD(sims_9, peds, names(peds))
segments_12 <- lengthIBD(sims_12, peds, names(peds))
segments_15 <- lengthIBD(sims_15, peds, names(peds))

ibdrel_unilineal_test <- list(
  segments_0 = segments_0,
  segments_2 = segments_2,
  segments_5 = segments_5,
  segments_7 = segments_7,
  segments_9 = segments_9,
  segments_12 = segments_12,
  segments_15 = segments_15,
  peds = peds
)

usethis::use_data(ibdrel_unilineal_test, overwrite = T)
