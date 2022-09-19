library(tidyverse)
library(interp)
library(ggrepel)
library(plotly)
library(threejs)
#########################
# plot landscapes

setwd("/Users/eamenho/Qiu_Lab/sars_epistatis/")

# recipe: https://rviews.rstudio.com/2020/12/14/plotting-surfaces-with-r/
plotly.land <- function(proj_file) {
  fn <- read_csv(proj_file)
  origdata <- as.data.frame(list(x = fn$Dim_1, y = fn$Dim_2, z = fn$fit))
  grid <- with(origdata, interp::interp(x, y, z))
  plot_ly(x = grid$x, y = grid$y, z = grid$z, name = paste("landscape", proj_file)) %>% 
    add_surface()
}

plotly.land("cov-pca.csv")
plotly.land("cov-wuhan-pca.csv")
plotly.land("cov-beta-pca.csv")
plotly.land("cov-e484k-pca.csv")
plotly.land("cov-delta-pca.csv")
plotly.land("cov-n501y-pca.csv")


plot.land <- function(proj_file) {
  fn <- read_csv(proj_file)
  origdata <- as.data.frame(list(x = fn$Dim_1, y = fn$Dim_2, z = fn$fit))
  grid <- with(origdata, interp::interp(x, y, z))
  griddf <- subset(data.frame(x = rep(grid$x, nrow(grid$z)),
                              y = rep(grid$y, each = ncol(grid$z)),
                              z = as.numeric(grid$z)), !is.na(z))
  ggplot(griddf, aes(x, y, z = z)) +
    geom_contour_filled(bins = 10, alpha = 0.8) +
    labs(title = paste("Fitness landscape: ", proj_file)) +
    theme_bw()
}

plot.land("cov-pca.csv")
plot.land("cov-wuhan-pca.csv")
plot.land("cov-beta-pca.csv")
plot.land("cov-e484k-pca.csv")
plot.land("cov-delta-pca.csv")
plot.land("cov-n501y-pca.csv")
