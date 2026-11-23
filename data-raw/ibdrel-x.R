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

# Remove siblings
sib.idx <- which(names(peds) == "S")
peds = peds[-sib.idx]
annotation = annotation[-sib.idx]

# All configurations of leaves

peds.males <- peds
peds.females <- lapply(peds, pedtools::swapSex, ids = specifyIndividuals)

peds.both.right <- lapply(peds, function(x) {
  id = specifyIndividuals(x)[1]
  pedtools::swapSex(x, id)
})

peds.both.left <- lapply(peds, function(x) {
  id = specifyIndividuals(x)[2]
  pedtools::swapSex(x,id)
})

peds.all <- c(peds, peds.females, peds.both.right, peds.both.left)
annotation.all <- sapply(peds.all, annotatePedigree, includeInds = TRUE)
names(peds.all) <- annotation.all




# SIMULATE TRAINING DATA

seed = seq(1,length(peds.all))+20260911

sims_training <- ibdSimulations(peds.all, N = 1000, seed = seed, chrom = "X")
segments <- lengthIBD(sims_training, peds.all, annotation.all)

segments.filtered <- Filter(
  function(x) !all(lengths(x) == 0),
  segments
)

peds.filtered <- peds.all[names(segments.filtered)]

metadata = pedsMetadata(peds.filtered, includeInds = T)

ibdrel_x <- list(
  segments = segments.filtered,
  peds = peds.filtered,
  metadata = metadata
)

usethis::use_data(ibdrel_x, overwrite = T)
