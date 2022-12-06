#Sys.setenv("_R_CHECK_SYSTEM_CLOCK_" = 0)

# Create r package folder
usethis::create_package("C:/Users/pavla/Nextcloud/GitKraken/NCC")

# Copy in R folder the functions of the r package
setwd("C:/Users/pavla/Nextcloud/GitKraken/NCC")
devtools::document()
devtools::load_all()

# Build & check the package
devtools::build(pkg = "C:/Users/pavla/Nextcloud/GitKraken/NCC", path = NULL, binary = FALSE, manual = TRUE, vignettes = TRUE)
devtools::check_built(path = "C:/Users/pavla/Nextcloud/GitKraken/NCC", cran = TRUE, manual = TRUE)
devtools::build_manual(pkg = "C:/Users/pavla/Nextcloud/GitKraken/NCC", path = NULL)

#create vignette
usethis::use_vignette("my-vignette")

pkgdown::build_site(pkg = "C:/Users/pavla/Nextcloud/GitKraken/NCC")

# https://www.r-bloggers.com/2017/08/building-a-website-with-pkgdown-a-short-guide/
# https://r-pkgs.org/vignettes.html
