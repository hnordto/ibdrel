helpWidget = function(id) {
  tooltip(
    actionBttn(
      inputId = NS(id, "help"),
      label = NULL,
      icon = icon("question")
    ),
    title = "See help on this panel."
  )
}

checkLikelihood = function(likelihood) {
  if (all(likelihood == -Inf)) {
    validate(need(FALSE, "Waiting for input.\nNeed at least one segment with length >= cutoff."))
  }
}

likelihoodBoxUI = function(id) {

  tagList(
    tabBox(
      id = NS(id, "resTabs"),
      width = NULL,
      selected = "eqclass.detailed",
      collapsible = FALSE,

      title = div(
        div(
          "Result table",
        ),
        div(
          helpWidget(id)
        )
      ),

      tabPanel(
        title = "Sex-specific pedigree",
        value = "eqclass.detailed",
        div(
          class = "table-box",
          gt::gt_output(NS(id, "resTableDetailedEq"))
        )
      ),

      tabPanel(
        title = "Pedigree",
        value = "eqclass",
        div(
          class = "table-box",
          gt::gt_output(NS(id, "resTableEq"))
        )
      ),

      tabPanel(
        title = "Kappa",
        value = "kappa",
        div(
          class = "table-box",
          gt::gt_output(NS(id, "resTableKappa"))
        )
      ),

      tabPanel(
        title = "Kinship",
        value = "kinship",
        div(
          class = "table-box",
          gt::gt_output(NS(id, "resTableKinship"))
        )
      ),

      tabPanel(
        title = "Degree",
        value = "degree",
        div(
          class = "table-box",
          gt::gt_output(NS(id, "resTableDegree"))
        )
      )
    ),

    uiOutput(NS(id,"multiplePedBox"))
  )
}

likelihoodBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {

    # Set active tab
    observe({
      req(input$resTabs)
      data[["selectedTab"]] <- input$resTabs
      data[["selectedClass"]] <- input$selectClass
    })


    output$resTableDetailedEq = gt::render_gt({
      res <- "eqclass.detailed"
      checkLikelihood(data[["likelihoods"]][[res]])
      likelihoodTable(data, res, session)
    })

    output$resTableEq = gt::render_gt({
      res <- "eqclass"
      checkLikelihood(data[["likelihoods"]][[res]])
      likelihoodTable(data, res, session)
    })

    output$resTableKappa = gt::render_gt({
      res <- "kappa"
      checkLikelihood(data[["likelihoods"]][[res]])
      likelihoodTable(data, res, session)
    })

    output$resTableKinship = gt::render_gt({
      res <- "kinship"
      checkLikelihood(data[["likelihoods"]][[res]])
      likelihoodTable(data, res, session)
    })

    output$resTableDegree = gt::render_gt({
      res <- "degree"
      checkLikelihood(data[["likelihoods"]][[res]])
      likelihoodTable(data, res, session)
    })


    # Help
    observeEvent(input$help, {
      shinyalert(
        className = "helpbox",
        html = TRUE,
        text = read_file("help/likelihoodBox.html"),
        showConfirmButton = FALSE,
        closeOnClickOutside = TRUE
      )
    })

  })
}

# Format table

likelihoodTable <- function(data, resolution, session) {

  likelihoods = data[["likelihoods"]][[resolution]] # XX FIX XX
  mdists = data[["mdists"]][[resolution]]
  filter = data[["filter"]]
  normalize = data[["normalize"]]
  threshold = data[["mdistthreshold"]]

  metadata = data[["metadata"]]


  if (filter) {
    inlier.classes <- names(which(mdists < threshold))
    likelihoods <- likelihoods[inlier.classes]
  }

  if (normalize) {
    likelihoods = ibdrel::normalizeLikelihoods(likelihoods)
  }

  likelihoods <- sort(likelihoods, decreasing = TRUE)
  #classlabels <- sapply(names(likelihoods), ibdrel::classTranslator, resolution)

  mdists <- mdists[names(likelihoods)] # Sorted mdists
  isoutlier <- mdists < threshold

  metadata.subset <- metadata[metadata$class %in% names(likelihoods),]
  metadata.subset <- metadata.subset[!duplicated(metadata.subset$class),]
  metadata.subset <- metadata.subset[match(names(likelihoods), metadata.subset$class),] # Ensure correct sort

  # Check if class contains multiple pedigrees

  classlabels = sapply(names(likelihoods), function(x) {
    metadata.class <- metadata[metadata$class == x,]$class
    lab = ibdrel::classTranslator(x, resolution)
    if (length(metadata.class) > 1) {
      paste(lab, "[+]")
    } else {
      lab
    }
  })
  classlabels = firstup(classlabels) # First letter to uppercase

  classlabels <- setNames(classlabels, names(likelihoods))

  kappa0 = as.character(MASS::fractions(metadata.subset$kappa0))
  kappa1 = as.character(MASS::fractions(metadata.subset$kappa1))
  kappa2 = as.character(MASS::fractions(metadata.subset$kappa2))
  kinship = as.character(MASS::fractions(metadata.subset$kinship))
  degree = as.character(metadata.subset$degree)

  # Base dataframe, independent of resolution
  df <- data.frame(Rank = seq_along(likelihoods),
                   Outlier = isoutlier,
                   Class = names(likelihoods),
                   Likelihood = round(likelihoods, 2))

  if (resolution == "eqclass.detailed" || resolution == "eqclass") {
    df <- cbind(df,
                Kappa0 = kappa0,
                Kappa1 = kappa1,
                Kappa2 = kappa2,
                Kinship = kinship,
                Degree = degree)
    colnames(df)[5:8] <- c("\u03BA\u2080", "\u03BA\u2081", "\u03BA\u2082", "\u03C6")
  }

  if (resolution == "kappa") {
    df <- cbind(df,
                Kinship = kinship,
                Degree = degree)
    colnames(df)[5] <- c("\u03C6")
  }


  if (resolution == "kinship") {
    df <- cbind(df,
                Degree = metadata.subset$degree)
  }



  res <- df |>

    gt() |>

    text_transform(
      locations = cells_body(columns = c(Outlier)),
      fn = function(x) {
        ifelse(
          x,
          "<span>&#9989;</span>",
          "<span>&#10060;</span>"
        )
      }
    ) |>

    text_transform(
      locations = cells_body(columns = c(Class)),
      fn = function(x) {
        sprintf(
          "<div class='clickable-cell' id='rel_%s'>%s</div>",
          x, classlabels[x]
        )
      }
    ) |>
    tab_style(
      style = cell_fill(color = "forestgreen", alpha = .25),
      locations = cells_body(columns = c(Class, Likelihood))
    ) |>
    gt_theme_538(quiet = TRUE)

  session$onFlushed(function() {
    session$sendCustomMessage(
      "initClickableCells",
      list(ids = df$Class, ns = session$ns("selectClass"))
    )
  }, once = TRUE)

  return (res)
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
