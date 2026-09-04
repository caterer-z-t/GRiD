# R_Package_lib <- "/nas/longleaf/home/catererz/R/x86_64-pc-linux-gnu-library/4.5"
# install.packages("hexSticker", lib = R_Package_lib)

# library(hexSticker, lib.loc = R_Package_lib)
# library(showtext, lib.loc = R_Package_lib)
# library(ggplot2, lib.loc = R_Package_lib)
# library(rstudioapi, lib.loc = R_Package_lib)

library(hexSticker)
library(showtext)
library(ggplot2)
library(rstudioapi)
library(magick)

## Loading Google fonts
font_add_google("Gochi Hand", "gochi")
showtext_auto()

# Get the script directory
outdir <- dirname(rstudioapi::getSourceEditorContext()$path)
bkg_img <- file.path(outdir, "grid_bkg.png")
print_file <- file.path(outdir, "grid_1200x1200.png")
web_file <- file.path(outdir, "grid_600x600.png")

# Common sticker parameters
common_params <- list(
  subplot = bkg_img,
  package = "GRiD",
  p_family = "gochi",
  p_color = "#FFFFFF",
  p_y = 0.4,
  s_x = 1.0,
  s_y = 1.2,
  s_width = 0.85,
  s_height = 1,
  h_fill = "#F47B7B",
  h_color = "#F47B7B",
  h_size = 0,
  asp = 1 # Ensures exact 1:1 square aspect ratio
)

# 1. Printable image (1200 x 1200 px)
do.call(sticker, c(common_params, list(
  p_size = 60,
  filename = print_file,
  out_width = 2,
  out_height = 2,
  dpi = 600
)))

# Make absolutely sure the printable image is 1200 x 1200
img_print <- image_read(print_file)

img_print <- image_extent(
  img_print,
  geometry = "1200x1200",
  gravity = "center",
  color = "none"
)

image_write(img_print, path = print_file)

# 2. Website image (600 x 600 px)
img_web <- image_read(print_file)

img_web <- image_resize(img_web, "600x600")

img_web <- image_extent(
  img_web,
  geometry = "600x600",
  gravity = "center",
  color = "none"
)

image_write(img_web, path = web_file)
