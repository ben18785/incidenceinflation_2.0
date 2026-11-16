# Shared settings and output readers for the medium tier (M*).

m_serial <- list(mean = 15.3, sd = 9.3)
# (m_omega, the serial weights, is set in each M test file: helper-validation-spec.R loads after this file)

# Values of column 'col' at onset day t, one per iteration, ordered by iteration,
# from the iterations after the first half.
m_series <- function(df, t, col, niter) {
  x <- df[col_of(df, "time_onset") == t & col_of(df, "iteration") > niter / 2, , drop = FALSE]
  col_of(x, col)[order(col_of(x, "iteration"))]
}
m_kept <- function(df, col, niter) {
  x <- df[col_of(df, "iteration") > niter / 2, , drop = FALSE]
  col_of(x, col)[order(col_of(x, "iteration"))]
}

