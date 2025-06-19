# Utils


# Plots ------------------------

outlierPlot <- function(features, obs_features, ped) {
  count = features[[ped]][["countPdf"]]
  total = features[[ped]][["totalPdf"]]

  count.obs = obs_features[[1]][["countPdf"]]
  total.obs = obs_features[[1]][["totalPdf"]]

  data = data.frame(count, total)
  data.obs = data.frame(count = count.obs, total = total.obs)

  ggplot() +
    geom_point(data = data, mapping = aes(x = count, y = total)) +
    geom_point(data = data.obs, mapping = aes(x = count, y = total),
               shape = 4, stroke = 2, size = 5, colour = "darkred") +
    labs(x = "Number of segments",
         y = "Total IBD segment length (cM)") +
    theme_minimal()
}

varScatterplot <- function(features, obs_features, ped, var1_selection, var2_selection) {
  var1 = features[[ped]][[var1_selection]]
  var2 = features[[ped]][[var2_selection]]

  var1.obs = obs_features[[1]][[var1_selection]]
  var2.obs = obs_features[[1]][[var2_selection]]

  data = data.frame(var1, var2)
  data.obs = data.frame(var1 = var1.obs, var2 = var2.obs)

  ggplot() +
    geom_point(data = data, mapping = aes(x = var1, y = var2)) +
    geom_point(data = data.obs, mapping = aes(x = var1, y = var2),
               shape = 4, stroke = 2, size = 5, colour = "darkred") +
    labs(x = var1_selection,
         y = var2_selection) +
    theme_minimal()
}



filterMetadata <- function(metadata,
                           targetCol) {
  metadata |>
    group_by(targetCol) |>
    slice(1)

}

relsAtLevel <- function(metadata,
                        targetCol,
                        targetRow) {

  metadata |>
    dplyr::filter(.data[[targetCol]] == targetRow) |>
    dplyr::select(rel) |>
    as.vector() -> res
  res[[1]]
}


resultTable <- function(metadata,
                        posteriors,
                        outliers,
                        lofs,
                        df,
                        aggLevel) {

  if (aggLevel == "rel") {
    mergeCol = "rel"
  } else if (aggLevel == "class") {
    mergeCol = "class"
  } else if (aggLevel == "kappa") {
    mergeCol = "kappa"
  } else if (aggLevel == "kinship") {
    mergeCol = "kinship"
  } else if (aggLevel == "degree") {
    mergeCol = "degree"
  }

  selectCol = c("Relationship", "Posterior", "Outlier", "LOF")


  results <- data.frame(Relationship = names(posteriors),
                        Posterior = round(as.numeric(posteriors), 4),
                        Outlier = outliers,
                        LOF = lofs)

  results <- merge(results, metadata, sort = FALSE,
                   by.x = "Relationship", by.y = mergeCol,
                   all.x = TRUE, all.y = FALSE)

  if (mergeCol == "class") {
    results$class <- results$Relationship
  }

  hideCol = setdiff(colnames(results), selectCol)

  #results$Outlier = ifelse(isTRUE(results$Outlier), "Yes", "No")

  n_classes <- length(unique(results$class))
  palette <- RColorBrewer::brewer.pal(min(max(n_classes, 3), 9), "Blues")

  if (n_classes > length(palette)) {
    palette <- colorRampPalette(palette)(n_classes)
  }

  results |>
    gt() |>
#    data_color(
#      columns = Relationship,
#      colors = col_factor(
#        palette = palette,
#        domain = unique(results$class)
#      )
#    ) |>
#    data_color(
#      columns = Outlier,
#      colors = scales::col_factor(
#        palette = c("darkgreen", "darkred"),
#        domain = c(TRUE, FALSE)
#      )
#    ) |>
    text_transform(
      locations = cells_body(columns = Relationship),
      fn = function(x) {
        lapply(x, function(val) {
          val_js <- gsub("'", "&#39", jsonlite::toJSON(val, auto_unbox = TRUE))
          val_html <- htmltools::htmlEscape(val)
          gt::html(
            paste0(
              "<a href='#' onclick='Shiny.setInputValue(\"name_clicked\", ", val_js, ", {priority: \"event\"})'>",
              val_html,
              "</a>"
            )
          )
        })
      }
    ) |>
    gt_theme_538() -> results


  results <- results |> cols_hide(columns = hideCol)
  results
}


resultTable.dep <- function(metadata, posteriors, outliers, mdists, df, slice) {


  if (!is.null(metadata)) {
    df |> # Ordering
      group_by(class) |>
      mutate(max_prob = max(Posterior)) |>
      arrange(desc(max_prob), desc(Posterior)) |>
      select(-max_prob) |>
      ungroup() -> df

  } else {
    df = results
    df$class = df$Relationship

    if (!is.null(slice)) {
      top_group <- unique(df$class)[1:slice]

      df <- df[df$class %in% top_group,]
    }

  }

  palette <- hue_pal()(length(unique(df$class)))


  df |>
    gt() |>
    data_color(
      columns = class,
      target_columns = Relationship,
      palette = palette
    ) |>
    data_color(
      columns = Outlier,
      colors = scales::col_factor(
        palette = c("green", "red"),
        domain = c("No", "Yes")
      )
    ) |>
    text_transform(
      locations = cells_body(columns = Relationship),
      fn = function(x) {
        lapply(x, function(val) {
          gt::html(
            paste0(
              "<a href='#' onclick=\"Shiny.setInputValue('name_clicked', '", val, "', {priority: 'event'})\">",
              val,
              "</a>"
            )
          )
        })
      }
    ) |>
    gt_theme_538() -> df

  # Formatting
  df$Relationship = str_to_sentence(df$Relationship, locale = "en")
  df$Outlier = ifelse(df$Outlier, "Yes", "No")

  df
}

# GT hide everything except
