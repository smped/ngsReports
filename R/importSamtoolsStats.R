x <- c(
    list.files("inst/extdata/", pattern = "samt", full.names = TRUE)
)



## Try figuring out how to get a boxplot of quality scores
getBxp <- function(df){
    df <- mutate(df, Q = cumsum(Count) / sum(Count))
    Q25 <- filter(df, Q > 0.25)$Score[1]
    Q50 <- filter(df, Q > 0.5)$Score[1]
    Q75 <- filter(df, Q > 0.75)$Score[1]
    iqr <- Q75 - Q25
    lwr <- max(min(df$Score), Q25 - 1.5*iqr)
    upr <- min(max(df$Score), Q75 + 1.5*iqr)
    mn <- with(df, sum(Score * Count) / sum(Count))
    df[1:2][1,] %>%
        mutate(lwr = lwr, Q25 = Q25, Q50 = Q50, Q75 = Q75, upr = upr, mn = mn)
}
