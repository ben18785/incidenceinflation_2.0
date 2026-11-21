# The share of cases in the heavier tail of the reporting delay, varying by
# onset week (S8 of the sensitivity analyses, sensitivity_extensions.tex of
# the paper). A share q_(w) of the cases of onset week w are reported after an
# exponential delay with mean 'mean', shared by every week, rather than the
# gamma delay:
#
#   S(t | t') = (1 - q_t) {1 - G(t' - t + 1/2)} + q_t exp{-(t' - t + 1/2) / mean},
#
# with q_t the share of the week of onset day t. The logit of the share is a
# random walk over onset weeks, with steps N(0, rw_sd^2); the logit of the
# share of the first week has a normal prior and rw_sd a half-normal one.
# With every share equal this is the single share of the main model.
#
# The walk is carried in the element 'tail' of the revision parameters, a
# named list of:
# - 'weight': the share for each onset week;
# - 'mean': the mean of the exponential delay;
# - 'day_first': the first onset time, from which the weeks are counted (see
#   tail_of_onset);
# - 'rw_sd': the standard deviation of the steps of the walk.
#
# It is allowed for only in the joint model of cases and deaths, where it is
# switched on by priors$reporting_effects$tail_weight_walk (see
# tail_walk_settings).

#' Settings of the random walk of the share of cases in the tail of the
#' reporting delay
#'
#' @param supplied TRUE or a named list of settings to change
#'
#' @return a named list of: 'first_mean' and 'first_sd' (the normal prior of
#' the logit of the share of the first onset week; by default the logit of
#' 0.036 and 1.2, whose quantiles are close to those of a beta(1, 19) prior
#' of a single share), 'step_sd_scale' (the scale of the half-normal prior of
#' the standard deviation of the weekly steps; 0.15), 'initial_rw_sd' (its
#' starting value; 0.1) and 'rw_sd_step' (the standard deviation of the
#' changes proposed in its log; 0.3, and 0 holds it at its starting value)
tail_walk_settings <- function(supplied) {
  walk <- list(first_mean=stats::qlogis(0.036), first_sd=1.2, step_sd_scale=0.15,
               initial_rw_sd=0.1, rw_sd_step=0.3)
  if(isTRUE(supplied))
    supplied <- list()
  if(!is.list(supplied))
    stop("tail_weight_walk must be TRUE, FALSE or a named list of settings.")
  unknown <- setdiff(names(supplied), names(walk))
  if(length(unknown) > 0)
    stop(paste0("Unknown settings for tail_weight_walk: ",
                paste(unknown, collapse=", "), "."))
  walk[names(supplied)] <- supplied
  if(any(vapply(walk, function(x) length(x) != 1 || !is.numeric(x) || !is.finite(x),
                logical(1))))
    stop("The settings of tail_weight_walk must each be a single number.")
  if(any(unlist(walk[c("first_sd", "step_sd_scale", "initial_rw_sd")]) <= 0) |
     walk$rw_sd_step < 0)
    stop("The scales and starting value of tail_weight_walk must be positive, and its step zero or more.")
  walk
}

#' Starting values of the random walk of the share of cases in the tail
#'
#' @param tail the element 'tail' of the initial revision parameters: a named
#' list of 'weight' (a single share, or one for each onset week) and 'mean'
#' @param onset_times the onset times
#' @param walk as returned by \code{tail_walk_settings}
#'
#' @return the element 'tail' of the revision parameters, with a share for
#' each onset week (see the top of deaths-tail-walk.R)
initial_tail_walk <- function(tail, onset_times, walk) {
  n_weeks <- floor((max(onset_times) - min(onset_times)) / 7) + 1
  weight <- as.numeric(tail$weight)
  if(length(weight) != 1 & length(weight) != n_weeks)
    stop(paste0("The tail of the reporting delay must start from a single share or one for each of the ",
                n_weeks, " onset weeks."))
  if(any(!is.finite(weight)) || any(weight <= 0 | weight >= 1))
    stop("The shares of the tail of the reporting delay must start strictly between 0 and 1.")
  if(length(tail$mean) != 1 || !is.finite(tail$mean) || tail$mean <= 0)
    stop("The mean of the tail of the reporting delay must start from a positive number.")
  list(weight=rep(weight, length.out=n_weeks), mean=tail$mean,
       day_first=min(onset_times), rw_sd=walk$initial_rw_sd)
}

#' Log-prior of the random walk of the logit of the share in the tail
#'
#' @param x the logit of the share of each onset week
#' @param rw_sd the standard deviation of the steps
#' @param walk as returned by \code{tail_walk_settings}
#'
#' @return a log-density
tail_walk_log_prior <- function(x, rw_sd, walk) {
  stats::dnorm(x[1], walk$first_mean, walk$first_sd, log=TRUE) +
    sum(stats::dnorm(diff(x), 0, rw_sd, log=TRUE))
}

#' Log of the acceptance ratio of a change in the mean of the tail, proposed
#' on the log scale
#'
#' @param proposed,current the proposed and current means
#' @param logl_proposed,logl their log-likelihoods
#' @param prior named list of 'tail_mean_mu' and 'tail_mean_sigma', the mean
#' and sd of the gamma prior of the mean
#'
#' @return a log-ratio
tail_mean_log_ratio <- function(proposed, current, logl_proposed, logl, prior) {
  logl_proposed - logl +
    dgamma_mean_sd(proposed, prior$tail_mean_mu, prior$tail_mean_sigma, log=TRUE) -
    dgamma_mean_sd(current, prior$tail_mean_mu, prior$tail_mean_sigma, log=TRUE) +
    # Jacobian of the transformation to the log scale
    log(proposed) - log(current)
}

#' The updates of the random walk of the share in the tail, its standard
#' deviation and the mean of the tail, given a log-likelihood
#'
#' In each step: the walk of the logit of the shares is drawn whole by
#' elliptical slice sampling given the standard deviation of its steps (see
#' \code{elliptical_slice}); the standard deviation is then changed with the
#' walk held and with the steps scaled with it, as that of the share of the
#' backlog is (see \code{backlog_slice_steps}), unless walk$rw_sd_step is
#' zero; and the mean of the tail is
#' changed by a random-walk Metropolis step on the log scale.
#'
#' @param tail the element 'tail' of the revision parameters (see the top of
#' deaths-tail-walk.R)
#' @param walk as returned by \code{tail_walk_settings}
#' @param prior named list of 'tail_mean_mu' and 'tail_mean_sigma'
#' @param mean_step the standard deviation of the changes proposed in the log
#' of the mean
#' @param f_logl the log-likelihood of a tail
#' @param n_steps the number of steps
#'
#' @return the tail updated
tail_walk_steps <- function(tail, walk, prior, mean_step, f_logl, n_steps=1) {
  f_accept <- function(log_ratio)
    !is.na(log_ratio) && log(stats::runif(1)) < log_ratio
  f_half_normal <- function(x) stats::dnorm(x, 0, walk$step_sd_scale, log=TRUE)
  f_tail <- function(x, mean=tail$mean) {
    out <- tail
    out$weight <- stats::plogis(x)
    out$mean <- mean
    out
  }
  n <- length(tail$weight)
  x <- stats::qlogis(tail$weight)
  logl <- f_logl(tail)
  for(step in seq_len(n_steps)) {
    # the whole walk, given the standard deviation of its steps
    sd_steps <- tail$rw_sd
    out <- elliptical_slice(
      x, rep(walk$first_mean, n),
      function() cumsum(c(stats::rnorm(1, 0, walk$first_sd),
                          stats::rnorm(n - 1, 0, sd_steps))),
      function(x) f_logl(f_tail(x)), logl)
    x <- out$x
    logl <- out$logl
    tail <- f_tail(x)
    if(n > 1 & walk$rw_sd_step > 0) {
      # its standard deviation with the walk held: no likelihood is needed
      current <- tail$rw_sd
      proposed <- current * exp(stats::rnorm(1, 0, walk$rw_sd_step))
      log_ratio <- tail_walk_log_prior(x, proposed, walk) -
        tail_walk_log_prior(x, current, walk) +
        f_half_normal(proposed) - f_half_normal(current) +
        log(proposed) - log(current)
      if(f_accept(log_ratio))
        tail$rw_sd <- proposed
      # and with its steps scaled with it, the standardised steps held
      current <- tail$rw_sd
      proposed <- current * exp(stats::rnorm(1, 0, walk$rw_sd_step))
      x_proposed <- c(x[1], x[1] + cumsum(diff(x) * proposed / current))
      logl_proposed <- f_logl(f_tail(x_proposed))
      log_ratio <- logl_proposed - logl +
        f_half_normal(proposed) - f_half_normal(current) +
        log(proposed) - log(current)
      if(f_accept(log_ratio)) {
        x <- x_proposed
        logl <- logl_proposed
        tail <- f_tail(x)
        tail$rw_sd <- proposed
      }
    }
    # the mean of the tail, shared by every onset week
    proposed <- tail$mean * exp(stats::rnorm(1, 0, mean_step))
    logl_proposed <- f_logl(f_tail(x, proposed))
    if(f_accept(tail_mean_log_ratio(proposed, tail$mean, logl_proposed, logl,
                                    prior))) {
      tail$mean <- proposed
      logl <- logl_proposed
    }
  }
  tail
}

#' Updates the share of cases in the tail of the reporting delay for each
#' onset week, the standard deviation of its random walk and the mean of the
#' tail, in the joint model of cases and deaths
#'
#' This takes the place of the joint update of the share and mean in
#' \code{deaths_tail_step}. The log-likelihood is that of every term which
#' involves the reporting delay of cases, over all onset times and the
#' uncorrected series (see \code{deaths_reporting_logp}): the cases, fatal
#' and not, newly reported in each report, the deaths newly shown where these
#' depend on it, and the wrong records entered and removed. The updates are
#' those of \code{tail_walk_steps}.
#'
#' @inheritParams deaths_tail_step
#'
#' @return the revision parameters with element 'tail' updated
deaths_tail_walk_step <- function(state, cases_true, reporting_by_onset,
                                  revision_parameters, settings) {
  prior <- settings$tail_prior
  # the reporting delay is held, so its part is calculated once for each
  # value of the tail, for all the onset days of a week
  cache <- new_lag_cache()
  f_logl <- function(tail) {
    proposed <- revision_parameters
    proposed$tail <- tail
    proposed$lag_cache <- cache
    deaths_reporting_logp(state, cases_true, reporting_by_onset, proposed)
  }
  revision_parameters$tail <- tail_walk_steps(
    revision_parameters$tail, prior$tail_weight_walk, prior,
    settings$tail_steps$tail_mean_step, f_logl, settings$n_parameter_steps)
  revision_parameters
}

#' The shares of the tail for each onset week, for storage
#'
#' @param tail the element 'tail' of the revision parameters (see the top of
#' deaths-tail-walk.R)
#'
#' @return a named list of 'tail_weight_week_<w>' for each onset week w and
#' 'tail_weight_rw_sd'
tail_walk_draws <- function(tail) {
  out <- stats::setNames(as.list(tail$weight),
                         paste0("tail_weight_week_", seq_along(tail$weight)))
  out$tail_weight_rw_sd <- tail$rw_sd
  out
}
