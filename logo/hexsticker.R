install.packages("hexSticker")

library(hexSticker)
library(showtext)
library(ggplot2)
library(rstudioapi)

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
p <- ggplot()
p <- p + theme_void() + theme_transparent() # + p

## use the ggplot2 example

sticker(
    p,
    package = "Unnamed", p_size = 18,
    s_x = 1, s_y = .75, s_width = 1.3, s_height = 1,
    h_fill = "#3E8A4A", h_color = "black",
    p_family = "gochi",
    filename = file.path(outdir,"unnamed.png"),
    # url = "doi.org/10.1093/bioinformatics/btaf540",
    u_size = 4, # smaller font size
    u_color = "black",
    u_family = "gochi"
)
