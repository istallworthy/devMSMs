# ~/.Rprofile
if (interactive()) {
  View <- function(x, title) {
    if (is.data.frame(x) || is.matrix(x) || is.list(x)) {
      utils::View(x, title)  # This should open in RStudio's viewer
    } else {
      print("Object is not viewable in RStudio Viewer")
    }
  }
}



