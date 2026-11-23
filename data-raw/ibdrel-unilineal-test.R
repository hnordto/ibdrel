# Requires ibdrel::ibdrel_unilineal

library(ibdrel)

peds <- ibdrel::ibdrel_unilineal$peds

annotation <- sapply(peds, annotatePedigree)

 # Change seed
seeds <- 1000:(length(peds)-1+1000)

sims_0 <- ibdSimulations(peds, N = 100, seed = NULL, keep = "nonzero", cutoff = 0)
sims_7 <- ibdSimulations(peds, N = 100, seed = NULL, keep = "nonzero", cutoff = 7)


# Note: Seeds contained i sims/segment objects

segments_0 <- lengthIBD(sims_0, peds, names(peds))
segments_7 <- lengthIBD(sims_7, peds, names(peds))

ibdrel_unilineal_test <- list(
  segments_0 = segments_0,
  segments_7 = segments_7,
  peds = peds
)

usethis::use_data(ibdrel_unilineal_test, overwrite = T)
