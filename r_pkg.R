#Sys.setenv("_R_CHECK_SYSTEM_CLOCK_" = 0)

# Create r package folder
usethis::create_package("D:/My Drive/GitKraken/NCC")

# Copy in R folder the functions of the r package
setwd("D:/My Drive/GitKraken/NCC")
devtools::document()
devtools::load_all()

# Build & check the package
devtools::build(pkg = "D:/My Drive/GitKraken/NCC", path = NULL, binary = FALSE, manual = TRUE, vignettes = TRUE)
devtools::check_built(path = "D:/My Drive/GitKraken/NCC", cran = TRUE, manual = TRUE, incoming = TRUE)
devtools::build_manual(pkg = "D:/My Drive/GitKraken/NCC", path = NULL)

#create vignette
usethis::use_vignette("my-vignette")

pkgdown::build_site(pkg = "D:/My Drive/GitKraken/NCC")

# https://www.r-bloggers.com/2017/08/building-a-website-with-pkgdown-a-short-guide/
# https://r-pkgs.org/vignettes.html
