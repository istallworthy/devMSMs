
#Building devMSMs website
#https://pkgdown.r-lib.org/articles/pkgdown.html


library(pkgdown)

usethis::use_pkgdown_github_pages()

use_pkgdown_github_pages(config_file = "_pkgdown.yml", destdir = "docs")
