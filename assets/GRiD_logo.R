R_Package_lib <- "/nas/longleaf/home/catererz/R/x86_64-pc-linux-gnu-library/4.5"
install.packages("hexSticker", lib = R_Package_lib)

library(hexSticker, lib.loc = R_Package_lib)
library(showtext, lib.loc = R_Package_lib)
library(ggplot2, lib.loc = R_Package_lib)
library(rstudioapi, lib.loc = R_Package_lib)

## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Gochi Hand", "gochi")
## Automatically use showtext to render text for future devices
showtext_auto()

outdir <- dirname(rstudioapi::getSourceEditorContext()$path)

counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)

# p <- ggplot(aes(x = mpg, y = wt), data = mtcars) +
#     geom_point()
p <- p + theme_void() + theme_transparent()

## use the ggplot2 example
sticker(
  p,
  package = "GRiD", p_size = 22, 
  s_x = 1, s_y = .75, s_width = 1.3, s_height = 1,
  h_fill = "#D91E36", h_color = "black",
  p_family = "gochi", 
  filename = file.path("grid_logo.png"),
  url = "Genomic Repeat inference from Depth",
  u_size = 4, # smaller font size
  u_color = "black",
  u_family = "gochi"
)