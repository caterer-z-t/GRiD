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

## Loading Google fonts
font_add_google("Gochi Hand", "gochi")
showtext_auto()

# Get the script directory
outdir <- dirname(rstudioapi::getSourceEditorContext()$path)

# Path to your background image in the same directory
bkg_img <- file.path(outdir, "grid_bkg.png")

## Generate sticker using grid_bkg.png as the image
sticker(
  subplot = bkg_img,
  package = "GRiD",
  p_size = 60,
  p_family = "gochi",
  p_color = "#FFFFFF",
  p_y = .4, # Text positioning above background image
  s_x = 1.0,
  s_y = 1.2,
  s_width = 0.85, # Adjusted width/height to fit cleanly inside hex
  s_height = 1,
  h_fill = "#F47B7B",
  h_color = "#F47B7B",
  h_size = 0,
  # url = "Genomic Repeat inference from Depth",
  # u_size = 9.3,
  # u_color = "#FFFFFF",
  # u_family = "gochi",
  # u_y = 0.1, # Footer text positioning
  filename = file.path(outdir, "grid_logo.png"),
  dpi = 600
)
