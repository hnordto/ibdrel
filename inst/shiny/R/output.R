# Plots ------------------------

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
               shape = 4, size = 5, colour = "red")
}


formatTable <- function(table) {
  table |>
    gt() |>
    data_color(columns = vars(class),
               colors = "relationship",
               colors = PALETTE)
}


resultTable <- function(metadata, posteriors, outliers, mdists, df, slice) {

  results = data.frame(Relationship = names(posteriors),
                       Posterior = round(as.numeric(posteriors),4),
                       Outlier = outliers,
                       Distance = mdists,
                       Distance_p = pchisq(mdists, df, lower.tail = F))

  if (!is.null(metadata)) {
    df = merge(results, metadata, sort = FALSE)
    df$Group = factor(df$class, ordered = TRUE) # Already ordered

    if (!is.null(slice)) {
      top_group <- unique(df$Group)[1:slice]

      df <- df[df$Group %in% top_group,]
    }


    df |> # Ordering
      group_by(Group) |>
      mutate(max_prob = max(Posterior)) |>
      arrange(desc(max_prob), desc(Posterior)) |>
      select(-max_prob) |>
      ungroup() -> df

  } else {
    df = results
    df$Group = df$Relationship

    if (!is.null(slice)) {
      top_group <- unique(df$Group)[1:slice]

      df <- df[df$Group %in% top_group,]
    }

  }

  palette <- hue_pal()(length(unique(df$Group)))

  # Formatting
  df$Relationship = str_to_sentence(df$Relationship, locale = "en")
  df$Outlier = ifelse(df$Outlier, "Yes", "No")

  df |>
    gt() |>
    data_color(
      columns = Group,
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
    gt_theme_538() -> df
  df
}
