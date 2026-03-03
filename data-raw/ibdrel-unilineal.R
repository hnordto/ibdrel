library(pedtools)
library(ibdrel)
library(dplyr)

l = identify_path_lengths(0,11)
rels = listRelationships(l, ignoreSex=T)
#rels <- removeSymmetries(rels)

# Define "close" relationships: Where relationships with all combinations of sex paths should be present
close.rels <- c(1,2,3)

rels_close <- rels |> filter(degree %in% close.rels)
rels_distant <- rels |> filter(!(degree %in% close.rels))

# Only keep full maternal/paternal lineages for "distant" relationships
#rels_paternal <- rels_distant |> group_by(type, degree, removal, nGen, cousDeg, half) |>
#  slice_head(n = 1)
#rels_maternal <- rels_distant |> group_by(type, degree, removal, nGen, cousDeg, half) |>
#  slice_tail(n = 1)

peds_close <- constructPedigrees(rels_close)
peds_distant <- constructPedigrees(rels_distant)

rels_close_variants = do.call(rbind.data.frame, lapply(peds_close, pedVariants))
rels_distant_variants = do.call(rbind.data.frame, lapply(peds_distant, pedVariants, onlyExtremes = TRUE))

rels <- rbind(rels_close_variants, rels_distant_variants)
peds <- constructPedigrees(rels)

annotation = sapply(peds, annotatePedigree)
names(peds) <- annotation


metadata = pedsMetadata(peds)


# SIMULATE TRAINING DATA

seed = seq(1,length(peds))+20260303

sims_training <- ibdSimulations(peds, N = 1000, seed = seed)
segments <- lengthIBD(sims_training, peds, annotation)

ibdrel_unilineal <- list(
  segments = segments,
  peds = peds
)

usethis::use_data(ibdrel_unilineal, overwrite = T)
