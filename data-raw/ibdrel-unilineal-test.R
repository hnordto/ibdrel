# Requires ibdrel::ibdrel_unilineal

peds <- ibdrel::ibdrel_unilineal$peds

annotation <- sapply(peds, annotatePedigree)

 # Change seed
seeds <- 1000:(length(peds)-1+1000)

sims <- ibdSimulations(peds, N = 100, seed = seeds)
segments <- lengthIBD(sims, peds, names(peds))

ibdrel_unilineal_test <- list(
  segments = segments,
  peds = peds
)

usethis::use_data(ibdrel_unilineal_test, overwrite = T)
