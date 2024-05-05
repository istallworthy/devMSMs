# Adaptation of gtools::permutations(2, repeats.allowed = TRUE)
perm2 <- function(r, v) {
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r %% 1) != 0) {
    stop("bad value of r")
  }
  if (!is.atomic(v)) {
    stop("v is non-atomic")
  }

  v <- unique(sort(v))
  n <- length(v)

  sub <- function(r, v) {
    if (r == 1) {
      return(matrix(v, ncol = 1))
    }

    inner <- Recall(r - 1, v)
    cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner),
      ncol = ncol(inner), nrow = nrow(inner) * n,
      byrow = TRUE
    ))
  }

  sub(r, v)
}

# Base R rewrite of stringr::str_count(); pattern can only be length 1
string_count <- function(string, pattern = "") {
  lengths(regmatches(string, gregexpr(pattern, string)))
}
