library(ggplot2)

source("R/simulate.R")


segmentDensity = function(segments, var = "length") {
  
  if (is.null(nrow(segments))) {
    segments = bind_rows(lapply(segments, data.frame), .id="id")
  } else {
    segments$id = factor(1)
  }
  
  if (var == "length") {
    p <- ggplot(segments, aes(x = length)) +
      geom_density(aes(fill = id, colour = id), alpha = .75, linewidth = 1.25) +
      facet_grid(rows = vars(id)) +
      theme_classic() +
      labs(x = "Individual IBD segment length (cM)",
           y = "Density") +
      theme(strip.background = element_blank(),
            strip.text.y = element_blank(),
            legend.position = "none")
  } else if (var == "number") {
    p <- segments |> 
      group_by(id,sim) |> 
      summarise(count = n()) |> 
      ggplot(aes(x = count)) +
        geom_density(aes(fill = id), alpha = .75) +
        facet_grid(rows = vars(id)) +
        theme_classic() +
        labs(x = "Number of IBD segments",
             y = "Density") +
        theme(strip.background = element_blank(),
              strip.text.y = element_blank(),
              legend.position = "none")
      
  } else if (var == "total") {
    p <- segments |> 
      group_by(id,sim) |> 
      summarise(total = sum(length)) |> 
      ggplot(aes(x = total)) +
        geom_density(aes(fill = id), alpha = .75) +
        facet_grid(rows = vars(id)) +
        theme_classic() +
        labs(x = "Total length IBD (cM)",
             y = "Density") +
        theme(strip.background = element_blank(),
              strip.text.y = element_blank(),
              legend.position = "none")
  }
  
  return (p)
  
}


