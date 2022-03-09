# Install development version from GitHub
devtools::install_github("pavlakrotka/NCC", build = TRUE, force=T)
# Run once to configure your package to use pkgdown
usethis::use_pkgdown()
# Run to build the website
pkgdown::build_site() 
