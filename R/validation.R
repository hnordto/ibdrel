
topPred = function(pred, n = NULL) {
  if (is.null(n)) {
    n = nrow(pred[[1]])
  }

  tp = lapply(pred, function(mat) {
    (apply(mat, 1, function(row) {
      colnames(mat)[order(row, decreasing = TRUE)[1:n]]
    }))
  })

  return (tp)
}

#'
#'@export
sensitivity <- function(pred, true, k=1) {
  classes <- colnames(pred)
  true <- factor(true, levels = classes)
  true_chr <- as.character(true)

  if (any(is.na(true))) {
    stop()
  }

  res <- lapply(k, function(kk) {

    topk_idx <- t(apply(pred, 1, function(x) {
      order(x, decreasing = TRUE)[1:kk]
    }))

    topk_class <- matrix(
      classes[topk_idx],
      nrow = nrow(pred),
      ncol = kk
    )

    res.k <- lapply(classes, function(c) {
      is_true_c <- (true_chr == c)

      is_pred_in_top_k <- apply(topk_class, 1, function(x) {
        c %in% x
      })

      TP <- sum(is_true_c & is_pred_in_top_k)
      FN <- sum(is_true_c & !is_pred_in_top_k)

      sens <- if ((TP+FN) > 0) TP / (TP+FN) else NA_real_

      data.frame(
        class = c,
        sensitivity = sens,
        k = kk
      )

    })

    do.call(rbind, res.k)

  })

  do.call(rbind, res)

}

ppv <- function(pred, true) {
  classes <- colnames(pred)
  true <- factor(true, levels = classes)
  true_chr <- as.character(true)

  top_idx <- t(apply(pred, 1, function(x) {
    order(x, decreasing = TRUE)[1]
  }))

  top_class = classes[top_idx]

  res <- lapply(classes, function(c) {
    is_true_c <- (true_chr == c)
    is_pred_c <- (top_class == c)

    TP <- sum(is_true_c & is_pred_c)
    FP <- sum(!is_true_c & is_pred_c)

    ppv <- TP / (TP+FP)

    data.frame(
      class = c,
      ppv = ppv
    )

  })

  do.call(rbind, res)


}

reliability <- function(pred, true, nbins = 10) {
  classes <- colnames(pred)
  true <- factor(true, levels = classes)

  res <- lapply(classes, function(c) {
    y_bin <- as.integer(true == c)
    p_hat = pred[,c]

    bins <- cut(
      p_hat,
      breaks = seq(0, 1, length.out = nbins + 1),
      include.lowest = TRUE,
      right = TRUE
    )

    tibble(
      y_bin = y_bin,
      p_hat = p_hat,
      bin = bins
    ) |>
      group_by(bin) |>
      summarise(
        bin_center = mean(p_hat), # Empirisk posisjon: Vil variere mellom klasser,
        mean_pred = mean(p_hat),
        obs_rate = mean(y_bin),
        n = n(),
        .groups = "drop"
      ) |>
      mutate(class = c)

  })

  do.call(rbind, res)
}

model.logloss <- function(segments.test, features, cutoff) {
  metadata = pedsMetadata(ibdrel_unilineal$peds)

  true.eqdetailed <- factor(trueClasses(segments.test, metadata, "eqclass.detailed"))
  true.eq <- factor(trueClasses(segments.test, metadata, "eqclass"))
  true.deg <- factor(trueClasses(segments.test, metadata, "degree"))

  res <- data.frame(features = character(),
                    logloss.eqdetailed = numeric(),
                    logloss.eq = numeric(),
                    logloss.eeg = numeric())

  for (i in 1:length(features)) {
    feature <- features[[i]]

    model = fitModel(features = feature, cutoff = cutoff, collapseDistant = FALSE)

    pred = predict(segments.test, model, sort = FALSE)

    pred.eqdetailed <- pred$eqclass.detailed
    pred.eq <- pred$eqclass
    pred.deg <-pred$deg

    true.eqdetailed <- factor(true.eqdetailed, levels = colnames(pred.eqdetailed))
    true.eq <- factor(true.eq, levels = colnames(pred.eq))
    true.deg <- factor(true.deg, levels = colnames(pred.deg))

    w.eqdetailed <- 1 / table(true.eqdetailed)[as.character(true.eqdetailed)]
    w.eq <- 1 / table(true.eq)[as.character(true.eq)]
    w.deg <- 1 / table(true.deg)[as.character(true.deg)]

    logloss.eqdetailed <- SLmetrics::weighted.logloss(actual = true.eqdetailed, response = pred.eqdetailed, w = w.eqdetailed)
    logloss.eq <- SLmetrics::weighted.logloss(actual = true.eq, response = pred.eq, w = w.eq)
    logloss.deg <- SLmetrics::weighted.logloss(actual = true.deg, response = pred.deg, w = w.deg)

    res <- rbind(res, data.frame(features = toString(feature),
                                 logloss.eqdetailed = logloss.eqdetailed,
                                 logloss.eq = logloss.eq,
                                 logloss.deg = logloss.deg))
  }

  return (res)
}

smith.prepareTrain <- function(cutoff, path) {
  metadata = pedsMetadata(ibdrel_unilineal$peds)

  data = aggregateSegments(ibdrel_unilineal$segments, metadata, "degree", T)

  # Should be trained with 'count' + 'total'
  features <- lapply(data, prepareFeatures, featureSel = c("count", "total"), cutoff = cutoff)

  genomeLength <- 3391.354

  counts.train <- lapply(features, "[[", "count")
  total.train <- lapply(features, "[[", "total")

  # Degree 1 contains only PO. Smith Bayes classifier needs variation in training data to construct PDFs
  # Add small random variation
  counts.train[[1]] <- counts.train[[1]] + runif(length(counts.train[[1]]), -0.25,0.25)
  total.train[[1]] <- total.train[[1]] + runif(length(total.train[[1]]), -0.25, 0.25)



  if (length(counts.train) != length(total.train)) {
    stop()
  }

  for (i in 1:length(counts.train)) {
    counts <- counts.train[[i]]
    total <- total.train[[i]]

    df.tmp <- data.frame("ibdprop" = total /genomeLength,
                         "ibdcount" = counts,
                         "degree" = i,
                         "name" = paste0("degree-", i, "-", "pair-", seq(1,length(counts))))

    if (i == 1) {
      df.train <- df.tmp
    } else {
      df.train <- rbind(df.train, df.tmp)
    }
  }

  for (deg in seq(1,7)) {
    write.table(df.train[df.train$degree == deg,],
                file = paste0(path, "train_", cutoff, "/train_",deg,".csv"),
                sep = ",",
                col.names = F, quote = F, row.names = F)
  }

  cat(nrow(df.train[df.train$degree %in% c(1,2,3,4,5,6,7),]))

}

smith.prepareTest <- function(cutoff, path) {
  metadata = pedsMetadata(ibdrel_unilineal_test$peds)

  testdata = ibdrel_unilineal_test[[paste0("segments_",cutoff)]]
  data = aggregateSegments(testdata, metadata, "degree", T)

  features <- lapply(data, prepareFeatures, featureSel = c("count", "total"), cutoff = cutoff)

  genomeLength <- 3391.354

  counts.test <- lapply(features, "[[", "count")
  total.test <- lapply(features, "[[", "total")

  counts.test[[1]] <- counts.test[[1]] + runif(length(counts.test[[1]]), -0.25,0.25)
  total.test[[1]] <- total.test[[1]] + runif(length(total.test[[1]]), -0.25, 0.25)

  if (length(counts.test) != length(total.test)) {
    stop()
  }

  for (i in 1:length(counts.test)) {
    counts <- counts.test[[i]]
    total <- total.test[[i]]

    df.tmp <- data.frame("ibdprop" = total / genomeLength,
                         "ibdcount" = counts,
                         "degree" = i,
                         "name" = paste0("degree-",i,"-","pair-",seq(1,length(counts))))
    if (i == 1) {
      df.test <- df.tmp
    } else {
      df.test <- rbind(df.test, df.tmp)
    }
  }

  for (deg in seq(1,7)) {
    write.table(df.test[df.test$degree==deg,],
                file = paste0(path,"test_",cutoff,"/test_",deg,".csv"),
                sep = ",", col.names = F, quote = F, row.names = F)
  }

  cat(nrow(df.test[df.test$degree %in% c(1,2,3,4,5,6,7),]))

}

smith.readResult <- function(path) {
  input <- readLines(path)

  m <- matrix(NA, nrow = 7, ncol = 7,
              dimnames = list(Prediction = c(seq(1,7)),
                              Reference = c(seq(1,7))))
  for (deg in seq(1,7)) {
    corr = sub(".*\\(([^()]*)\\).*", "\\1", input[2*deg])
    x <- input[2*deg+1]
    matches <- regmatches(x, gregexpr(":\\s*\\d+", x))[[1]]
    numbers <- as.numeric(sub(":\\s*", "", matches))


    numbers <- append(numbers, corr, after = deg - 1)

    m[,deg] = numbers
  }

  m <- matrix(as.numeric(m), nrow = nrow(m), ncol = ncol(m))
  rownames(m) <- colnames(m) <- seq(1,7)
  m <- caret::confusionMatrix(m)

  m
}

ibdkinship.prepareTest <- function(cutoff, path, sims = NULL) {
  metadata = pedsMetadata(ibdrel_unilineal$peds)
  peds.to.include <- metadata |>
    filter(eqclass %in% c("L1", "hS", "A", "L2", "1C0R",
                          "1C1R", "2C0R", "2C1R", "3C")) |>
    select(rel) |> as.vector()
  distant.peds.to.include <- metadata |>
    filter(degree > 7) |>
    select(rel) |> as.vector()
  peds.to.include <- c(unique(peds.to.include$rel), unique(distant.peds.to.include$rel))

  #testdata = ibdrel_unilineal_test[[paste0("segments_",cutoff)]]
  #testdata <- testdata[peds.to.include]
  #data = aggregateSegments(testdata, metadata, "eqclass", T)

  peds <- ibdrel_unilineal$peds[peds.to.include]
  annotation <- sapply(peds, annotatePedigree)

  if (is.null(sims)) {
    sims <- ibdSimulations(peds, N = 100, keep = "nonzero", cutoff = 5)
  }




  ibdkinship_names <- c("pach" = "L1",
                        "hsibling" = "hS",
                        "avuncular" = "A",
                        "grand" = "L2",
                        "1cousin" = "1C0R",
                        "1cor" = "1C1R",
                        "2cousin" = "2C0R",
                        "2cor" = "2C1R",
                        "3cousin" = "3C0R",
                        "unrelated" = "distant")

  segments <- lapply(seq_along(sims), function(sidx) {
    s = sims[[sidx]]
    ped = peds[[sidx]]
    ibdsim2::findPattern(s, pattern = list(carriers = identifyLeaves(ped)),
                         cutoff = cutoff, unit = "cm")
  })
  names(segments) <- names(sims)

  lookup = lookupClass(NULL, "eqclass", "rel", metadata, T)
  group = lookup[names(segments)]
  merged_segments <- lapply(split(seq_along(segments), group), function(i) unlist(segments[i], recursive = F))


  for (i in 1:length(merged_segments)) {
    sim_rel = merged_segments[[i]]
    rel = names(merged_segments)[i]
    ibdkinship_name = names(ibdkinship_names)[which(ibdkinship_names == rel)]

    for (chr in seq(1,22)) {
      rel_chr <- lapply(seq_along(sim_rel), function(id) {
        x = sim_rel[[id]]
        d <- as.data.frame(x)
        d <- d[d$chr==chr,]

        if (nrow(d) > 0) {
          data.frame(sample = id, ibd = d$endCM-d$startCM)
        } else {
          NULL
        }


      })
      rel_chr <- do.call(rbind, rel_chr)

      filename = paste0("chr",chr,"-",ibdkinship_name,".txt")
      writepath = paste0(path, filename)

      write.table(rel_chr, file = writepath, col.names = F,
                  quote = F, row.names = F)

    }

  }


}
