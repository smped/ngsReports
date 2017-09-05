# A commonly used function for setting up dendrograms for interactive plots
ggdend <- function(df) {
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = df, ggplot2::aes_string("x","y", xend = "xend", yend = "yend")) +
    ggdendro::theme_dendro()
}
