library(ibdrel)



rels <- readRDS("~/ibdclassifier/inst/data/unilineal_relationships_degree11.rds")

rels_paternal <- rels |> group_by(type, degree, removal, nGen, cousDeg, half) |>
  slice_head(n = 1)
rels_maternal <- rels |> group_by(type, degree, removal, nGen, cousDeg, half) |>
  slice_tail(n = 1)
rels <- rbind(rels_paternal, rels_maternal)
peds <- constructPedigrees(rels)
annotation <- annotatePedigrees(peds)

peds_training <- peds[1:77] # Only paternal, for simplicity
names(peds_training) <- names(segments_training_rel)


aggregate.segments <- function(segments, metadata, metadata.agg.column) {
  aggregatedSegments = list()

  for (i in 1:length(segments)) {
    rel.id = names(segments)[i]
    agg.id = metadata |>
      dplyr::filter(rel == rel.id) |>
      dplyr::select(metadata.agg.column) |>
      as.character()

    segments.rel <- segments[[i]]

    for (segment in segments.rel) {
      if (agg.id %in% names(aggregatedSegments)) {
        aggregatedSegments[[agg.id]][[length(aggregatedSegments[[agg.id]])+1]] <- segment
      } else {
        aggregatedSegments[[agg.id]][[1]] <- segment
      }
    }
  }

  return (aggregatedSegments)

}


##########################
# GENERATE TRAINING DATA #
##########################

# Load rels from file


sims_training <- ibdSimulations(peds, N = 1000, seed = NULL)

segments_df_training <- postprocessSimulations(sims_training, peds,
                                               annotation,
                                               cutoff = 0)
segments_df_training$sim <- paste0(segments_df_training$sim, "-", segments_df_training$kinship) # Unique simulation identifier, keep all segments intact
segments_df_training$kinship <- sub("-([pm]+)?($|\\s.*)", "", segments_df_training$kinship) # Ignore sex paths

segments_training_rel <- lengthIBD(segments_df_training)

# Aggregation

metadata = pedsMetadata(peds_training)

segments_training_donnelly <- aggregate.segments(segments_training_rel,
                                                 metadata,
                                                 "class")

segments_training_kappa <- aggregate.segments(segments_training_rel,
                                              metadata,
                                              "kappa")
segments_training_kinship <- aggregate.segments(segments_training_rel,
                                                metadata,
                                                "kinship")
segments_training_degree <- aggregate.segments(segments_training_rel,
                                               metadata,
                                               "degree")


# Save data

saveRDS(segments_training_rel, file = "inst/data/segments_unilineal_rel.rds")
saveRDS(segments_training_donnelly, file = "inst/data/segments_unilineal_donnelly.rds")
saveRDS(segments_training_kappa, file = "inst/data/segments_unilineal_kappa.rds")
saveRDS(segments_training_kinship, file = "inst/data/segments_unilineal_kinship.rds")
saveRDS(segments_training_degree, file = "inst/data/segments_unilineal_degree.rds")
saveRDS(peds_training, file = "inst/data/peds_unilineal.rds")





##########################
# GENERATE TEST DATA     #
##########################

sims_test <- ibdSimulations(peds, N = 100, seed = NULL)


segments_df_test <- postprocessSimulations(sims_test,
                                           peds,
                                           annotation,
                                           cutoff = 7)
segments_df_test$sim <- paste0(segments_df_test$sim, "-", segments_df_test$kinship)
segments_df_test$kinship <- sub("-([pm]+)?($|\\s.*)", "", segments_df_test$kinship) # Ignore sex paths

segments_test <- lengthIBD(segments_df_test)

saveRDS(segments_test, file = "inst/data/segments_unilineal_rel_test.rds")




