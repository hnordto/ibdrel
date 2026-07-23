helpWidget = function(id) {
  tooltip(
    actionBttn(
      inputId = NS(id, "help"),
      label = NULL,
      style = "jelly",
      size = "sm",
      icon = icon("question"),
      color = "royal"
    ),
    title = "See help on this panel"
  )
}

w_resolution = function(id) {
  tooltip(
    radioGroupButtons(
      inputId = NS(id,"resTabs"),
      label = NULL,
      choices = c("Sex" = "eqclass.detailed",
                  "Rel" = "eqclass",
                  "Kap" = "kappa",
                  "Kin" = "kinship",
                  "Deg" = "degree"),
      justified = TRUE,
      individual = TRUE,
      status = "olive"
    ),
    title = "Select resolution"
  )
}

checkLikelihood = function(likelihood) {
  if (all(likelihood == -Inf, na.rm = TRUE)) {
    validate(need(FALSE, "Waiting for input.\nNeed at least one segment with length above cutoff value."))
  }
}

likelihoodBoxUI = function(id) {
  bs4Card(
    id = NS(id, "likelihoodBox"),
    title = div(class = "box-title-flex",
      div(class = "leftcolumn", "Relatedness"),
      div(class = "rightcolumn", helpWidget(id))
    ),
    width = NULL,
    collapsible = FALSE,
    status = "olive",

    div(w_resolution(id)),

    DT::DTOutput(NS(id, "resTable"),
                 width = "fit-content"),
    footer = div(
      checkboxGroupInput(NS(id, "likelihoodTableSettings"), "",
                         c("Normalize" = "normalize",
                           "Filter" = "filter"),
                         selected = "normalize",
                         inline = TRUE)
    )
  )
}

likelihoodBoxServer = function(id, data) {
  moduleServer(id, function(input, output, session) {

    observe({
      data[["selectedTab"]] <- input$resTabs
    })


    table_data <- reactive({
      req(data[["likelihoods"]])
      res <- input$resTabs
      checkLikelihood(data[["likelihoods"]][[res]])
      computeLikelihoodTable(data, res)
    })

    output$resTable <- DT::renderDataTable({
      tab <- table_data()
      likelihoodTable(tab$df)
    })

    observeEvent(input$resTable_rows_selected, {
      idx <- input$resTable_rows_selected
      if (length(idx) == 1) {
        tab <- table_data()
        class_codes <- tab$class_codes
        selected_class <- class_codes[idx]
        data[["selectedClass"]] <- selected_class
      }
    })

    observe({
      data[["normalize"]] = "normalize" %in% input$likelihoodTableSettings
      data[["filter"]] = "filter" %in% input$likelihoodTableSettings
    })

    # Help
    observeEvent(input$help, {
      shinyalert(
        className = "helpbox",
        html = TRUE,
        text = read_file("help/likelihoodBox.html"),
        showConfirmButton = FALSE,
        closeOnClickOutside = TRUE,
        size = "m"
      )
    })

  })
}

# Table maps

computeLikelihoodTable <- function(data, resolution) {
  likelihoods = data[["likelihoods"]][[resolution]]
  mdists = data[["mdists"]][[resolution]]
  filter = data[["filter"]]
  normalize = data[["normalize"]]
  threshold = data[["mdistthreshold"]][[resolution]]
  metadata = data[["metadata"]]

  if (normalize) {
    likelihoods = ibdrel::normalizeLikelihoods(likelihoods)
  }


  if (filter) {
    inlier.classes <- names(which(mdists < threshold))
    likelihoods <- likelihoods[inlier.classes]
  }


  likelihoods <- sort(likelihoods, decreasing = TRUE)
  #classlabels <- sapply(names(likelihoods), ibdrel::classTranslator, resolution)

  mdists <- mdists[names(likelihoods)] # Sorted mdists
  threshold <- threshold[names(likelihoods)] # Sorter outlier thresholds
  isoutlier <- mdists < threshold

  metadata.subset <- metadata[metadata$class %in% names(likelihoods),]
  metadata.subset <- metadata.subset[!duplicated(metadata.subset$class),]
  metadata.subset <- metadata.subset[match(names(likelihoods), metadata.subset$class),] # Ensure correct sort

  # Check if class contains multiple pedigrees

  classlabels = sapply(names(likelihoods), function(x) {
    metadata.class <- metadata[metadata$class == x,]$class
    ibdrel::classTranslator(x, resolution)
  })
  classlabels = firstup(classlabels) # First letter to uppercase

  classlabels <- setNames(classlabels, names(likelihoods))

  if (resolution == "eqclass.detailed") {
    eqclass = as.character(metadata.subset$eqclass)
    eqclass = sapply(eqclass, function(x) {
      classTranslator(x, "eqclass")
    })
  }

  kappa0 = as.character(MASS::fractions(metadata.subset$kappa0))
  kappa1 = as.character(MASS::fractions(metadata.subset$kappa1))
  kappa2 = as.character(MASS::fractions(metadata.subset$kappa2))

  kappa = as.character(paste0("(",kappa0, ",", kappa1, ",", kappa2, ")"))
  kinship = as.character(MASS::fractions(metadata.subset$kinship))
  degree = as.character(metadata.subset$degree)

  # Base dataframe, independent of resolution
  df <- data.frame(Rank = seq_along(likelihoods),
                   Outlier = isoutlier,
                   Class = names(likelihoods),
                   Conf = round(likelihoods, 2))

  if (resolution == "eqclass.detailed") {
    df <- cbind(df,
                Rel = eqclass,
                Kappa = kappa,
                Kinship = kinship,
                Degree = degree)
    colnames(df)[6:7] <- c("\u03BA", "\u03C6")
    class_rename = "Sex"
  }

  if (resolution == "eqclass") {
    df <- cbind(df,
                Kappa = kappa,
                Kinship = kinship,
                Degree = degree)
    colnames(df)[6:7] <- c("\u03BA", "\u03C6")
    class_rename = "Rel"
  }

  if (resolution == "kappa") {
    df <- cbind(df,
                Kinship = kinship,
                Degree = degree)
    colnames(df)[6] <- c("\u03C6")
    class_rename = "\u03BA"
  }


  if (resolution == "kinship") {
    df <- cbind(df,
                Degree = metadata.subset$degree)
    class_rename = "\u03C6"
  }

  if (resolution == "degree") {
    class_rename = "Degree"
  }

  df$Outlier <- ifelse(
    df$Outlier,
    "&#9989;",
    "&#10060;"
  )

  df$Class <- classlabels[df$Class]
  colnames(df)[1:2] <- c(" ", " ")
  colnames(df)[3] <- class_rename

  #colnames(df)[colnames(df) == "Rank"] <- ""
  #colnames(df)[colnames(df) == "Outlier"] <- ""
  #colnames(df)[colnames(df) == "Class"] <- class_rename
  #colnames(df)[colnames(df) == "Likelihood"] = "Lik"

  list(
    df = df,
    class_codes = names(likelihoods),
    class_label = class_rename
  )
}

# Format table

likelihoodTable <- function(df) {


  res <- DT::datatable(
    df,
    rownames = FALSE,
    escape = FALSE,
    selection = "single",
    options = list(
      dom = "t",
      autoWidth = FALSE,
      ordering = FALSE,
      scrollX = TRUE,
      scrollY = "500px",
      scrollCollapse = TRUE,
      paging = FALSE,
      columnDefs = list(
        list(className = "bold-col", targets = 3),
        list(className = "dt-center", targets = "_all")
      )
    ),
    class = "compact stripe hover"
  )

  return (res)
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
