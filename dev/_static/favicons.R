library(magick)

# Source logo path
logo_path <- "/Users/zc/Library/CloudStorage/OneDrive-UCB-O365/Stanislabski/thesis_chapters/GRiD/assets/grid_logo.png"

# Destination directory
static_dir <- "/Users/zc/Library/CloudStorage/OneDrive-UCB-O365/Stanislabski/thesis_chapters/GRiD/docs/source/_static"

# Load source image
img <- image_read(logo_path)

# Helper function to crop square, resize, and save
generate_favicon <- function(image, size, filename) {
    # Trim transparent padding, fit square, and scale
    resized <- image %>%
        image_trim() %>%
        image_resize(paste0(size, "x", size)) %>%
        image_extent(paste0(size, "x", size), gravity = "Center", color = "none")

    image_write(resized, path = file.path(static_dir, filename))
}

# 1. Standard Web Favicons
generate_favicon(img, 16, "favicon-16x16.png")
generate_favicon(img, 32, "favicon-32x32.png")
generate_favicon(img, 32, "favicon.ico") # 32x32 PNG saved as .ico for modern browsers

# 2. Apple Touch Icon
generate_favicon(img, 180, "apple-touch-icon.png")

# 3. Android Chrome Icons
generate_favicon(img, 192, "android-chrome-192x192.png")
generate_favicon(img, 512, "android-chrome-512x512.png")

cat("All favicons successfully updated in _static/\n")
