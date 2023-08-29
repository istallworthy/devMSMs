
#Building devMSMs website
#https://pkgdown.r-lib.org/articles/pkgdown.html
#https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax#headings

library(pkgdown)
usethis::use_git(message = "Initial commit")
usethis::use_github(private = FALSE)
usethis::use_github_action("pkgdown")
usethis::use_pkgdown_github_pages()

# locally building site
pkgdown::build_site()

usethis::use_pkgdown_github_pages()

# use_pkgdown_github_pages(config_file = "_pkgdown.yml", destdir = "docs")

pkgdown::build_site()
devtools::document()

usethis::use_pkgdown()
usethis::use_github_action(url = "https://raw.githubusercontent.com/r-lib/actions/master/examples/pkgdown.yaml")

build_site_github_pages(
  pkg = ".",
  ...,
  dest_dir = "docs",
  clean = TRUE,
  install = FALSE,
  new_process = FALSE
)


