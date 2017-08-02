library(dplyr)
library(scales)
library(ggplot2)

df <- ngsReports::Per_base_sequence_content(fdl)
df <- dplyr::mutate(df,
                    Start = gsub("([0-9]*)-[0-9]*", "\\1", Base),
                    Start = as.integer(Start))
df <- mutate(df, colour = rgb(rescale(floor(A)+0.5*floor(G)), rescale(floor(T)+0.5*floor(G)), rescale(floor(C))))

basicStat <- Basic_Statistics(x) %>% dplyr::select(Filename, Longest_sequence)

df1 <- df %>% dplyr::right_join(basicStat, by = "Filename") %>%
  dplyr::select(Filename, Start, colour, Longest_sequence) %>%
  tidyr::spread(Start, colour)
unique(df$colour)

splitLengths <- df1 %>% split(f = .['Longest_sequence'])

dfInner <- splitLengths %>% lapply(function(k){
  k <- k %>% dplyr::select(Filename, dplyr::one_of(as.character(1:k$Longest_sequence[1]))) %>%
    t() %>%
    zoo::na.locf() %>%
    t() %>%
    as.data.frame()
}) %>% do.call(plyr::rbind.fill, .) %>% magrittr::set_rownames(.$Filename)

  dfLong <- dfInner %>% gather("Start", "Colour", 2:ncol(.))
  dfLong$Colour <- dfLong$Colour
  dfLong$Start <- as.integer(dfLong$Start)

  t <- fastqcReports::getSummary(x) %>% dplyr::filter(Category == "Per base sequence quality")
  t <- dplyr::full_join(dfInner["Filename"], t, by = "Filename")
  t$Filename <- with(t, factor(Filename, levels=Filename))
  key <- t["Filename"]

  cols <- scales::col_numeric(df$colour, domain = NULL)
  p <- ggplot(dfLong, aes(x = Start, y = Filename, fill = Colour)) + geom_tile() + theme(legend.position = "none") +
    scale_fill_manual(values = dfLong$Colour) +
    ggplot2::theme(panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())

  d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
    ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                              panel.background = element_blank(),
                                                              legend.position="none",
                                                              axis.title=element_blank(),
                                                              axis.text=element_blank(),
                                                              axis.ticks=element_blank())
  d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))

 BQheatmap <- plotly::subplot(d, p, widths = c( 0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "Sequencing Cycle"))


}

#do not cluster Filenames
if(!clusterNames){
  dfLong <- dfInner %>% tidyr::gather("Start", "Mean", 2:ncol(.))
  dfLong$Mean <- as.integer(dfLong$Mean)
  dfLong$Start <- as.integer(dfLong$Start)

  t <- getSummary(x) %>% dplyr::filter(Category == "Per base sequence quality")
  t <- dplyr::full_join(dfInner["Filename"], t, by = "Filename")
  t$Filename <- with(t, factor(Filename, levels=Filename))
  key <- t["Filename"]



  p <- ggplot2::ggplot(dfLong, aes(x = Start, y = Filename)) + ggplot2::geom_tile(aes(fill = Mean), color = "white", size = 5) +
    ggplot2::scale_fill_gradientn(colours = c(col["FAIL"], col["FAIL"],col["WARN"], col["WARN"], col["PASS"], col["PASS"]),
                                  values = rescale(c(0,20,20,30,30,40)),
                                  guide = "colorbar", limits=c(0, 40)) +
    ggplot2::theme(panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())

  d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
    ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                              panel.background = element_blank(),
                                                              legend.position="none",
                                                              axis.title=element_blank(),
                                                              axis.text=element_blank(),
                                                              axis.ticks=element_blank())
  d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


  BQheatmap <- plotly::subplot(d, p, widths = c(0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "Sequencing Cycle"))
}
}

