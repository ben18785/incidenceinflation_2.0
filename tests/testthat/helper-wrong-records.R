# the probability of the reports of one onset time, and of every path of the
# unobserved state, worked out in R with wrong records allowed for
wrong_paths_reference <- function(cases, deaths, fatal, survivors, probabilities,
                                  wrong) {
  n_reports <- length(cases)
  p_fatal <- probabilities$p_fatal
  if(is.null(p_fatal))
    p_fatal <- probabilities$p
  states <- expand.grid(r=0:max(cases), alive=0:wrong$max_erroneous,
                        dead=0:wrong$max_erroneous)
  states <- states[states$alive + states$dead <= wrong$max_erroneous, ]
  cases_all <- c(0, cases)
  deaths_all <- c(0, deaths)
  f_transition <- function(j, from, to) {
    a <- to$r - from$r
    survivors_from <- cases_all[j] - from$r - from$alive - from$dead
    survivors_to <- cases_all[j + 1] - to$r - to$alive - to$dead
    b <- survivors_to - survivors_from
    deaths_from <- deaths_all[j] - from$dead
    deaths_to <- deaths_all[j + 1] - to$dead
    deaths_new <- deaths_to - deaths_from
    if(a < 0 | b < 0 | survivors_from < 0 | to$r > fatal | survivors_to > survivors |
       deaths_from < 0 | deaths_new < 0 | deaths_to > to$r | deaths_from > from$r)
      return(0)
    p_wrong <- 0
    for(removed_alive in 0:from$alive)
      for(removed_dead in 0:from$dead) {
        entered_alive <- to$alive - from$alive + removed_alive
        entered_dead <- to$dead - from$dead + removed_dead
        if(entered_alive < 0 | entered_dead < 0)
          next
        p_wrong <- p_wrong +
          stats::dbinom(removed_alive, from$alive, wrong$w[j]) *
          stats::dbinom(removed_dead, from$dead, wrong$w[j]) *
          stats::dpois(entered_alive, wrong$lambda[j] * (1 - wrong$delta)) *
          stats::dpois(entered_dead, wrong$lambda[j] * wrong$delta)
      }
    alive <- from$r - deaths_from
    x <- 0:deaths_new
    p_wrong *
      stats::dbinom(a, fatal - from$r, p_fatal[j]) *
      stats::dbinom(b, survivors - survivors_from, probabilities$p[j]) *
      sum(stats::dbinom(x, a, probabilities$g[j]) *
            stats::dbinom(deaths_new - x, alive, probabilities$q[j]))
  }
  n_states <- nrow(states)
  start <- which(states$r == 0 & states$alive == 0 & states$dead == 0)
  l_transitions <- lapply(seq_len(n_reports), function(j) {
    m <- matrix(0, n_states, n_states)
    for(i in seq_len(n_states))
      for(k in seq_len(n_states))
        m[i, k] <- f_transition(j, states[i, ], states[k, ])
    m
  })
  alpha <- rep(0, n_states)
  alpha[start] <- 1
  for(j in seq_len(n_reports))
    alpha <- as.vector(alpha %*% l_transitions[[j]])
  list(logp=log(sum(alpha)), states=states, start=start,
       transitions=l_transitions)
}
