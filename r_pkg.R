

# Create r package folder
usethis::create_package("C:/Users/pavla/Desktop/NCC")

# Copy in R folder the functions of the r package 
setwd("X:/EU-PEARL/NCC")
devtools::document()
devtools::load_all()

# Build & check the package
devtools::build(pkg = "X:/EU-PEARL/NCC", path = NULL, binary = FALSE, manual = TRUE)
devtools::check_built(path = "X:/EU-PEARL/NCC", cran=TRUE, manual = TRUE)
devtools::build_manual(pkg = "X:/EU-PEARL/NCC", path = NULL)
