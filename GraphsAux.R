suppressPackageStartupMessages(library(ggplot2))
library(grid)

# Width (in pt) of the caption on the right; 0 if there is none

legend_width_pt <- function(p) {
  g <- ggplotGrob(p + theme(legend.position = "right"))
  i <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(i) == 0) return(0)
  convertWidth(sum(g$grobs[[i]]$widths), "pt", valueOnly = TRUE)
}

# Equalize PANEL width by adding padding to the right margin of the one with the narrowest caption

match_panel_width_right <- function(p1, p2,
                                    base_margin1 = margin(8, 10, 8, 8, unit = "pt"),
                                    base_margin2 = margin(8, 10, 8, 8, unit = "pt"),
                                    extra_pad_pt = 0) {
  w1 <- legend_width_pt(p1)
  w2 <- legend_width_pt(p2)
  target <- max(w1, w2)
  
  add1 <- max(0, target - w1 + extra_pad_pt)
  add2 <- max(0, target - w2 + extra_pad_pt)
  
  base_margin1[2] <- base_margin1[2] + unit(add1, "pt")
  base_margin2[2] <- base_margin2[2] + unit(add2, "pt")
  
  list(
    p1 + theme(plot.margin = base_margin1),
    p2 + theme(plot.margin = base_margin2)
  )
}



