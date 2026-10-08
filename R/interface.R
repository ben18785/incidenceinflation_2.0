# A common interface to the ways of fitting the model: the exact sampler
# (mcmc_single) and the relaxed model in Stan (fit_relaxed), fitted by
# Hamiltonian Monte Carlo or by an approximation. The same data, options and
# priors are given whichever is used, and the result has the same form.

#' Options of the model, common to all ways of fitting it
#'
#' @param is_negative_binomial whether the renewal model is negative binomial
#' (the default) or Poisson
#' @param Rt_model how Rt varies over onset times: "random_walk" (its log
#' follows a random walk between periods; relaxed model only), "periods" (it
#' is independent between periods, with a gamma prior) or "changepoints" (it
#' changes at unknown onset times)
#' @param Rt_period_width number of onset times in each period, for
#' "random_walk" and "periods"; periods are aligned to end on the last onset
#' time
#' @param Rt_index the period of each onset time, in place of periods of equal
#' width
#' @param Rt_max_changes for "changepoints" in the relaxed model, the largest
#' number of changes allowed
#' @param Rt_min_length for "changepoints" in the exact sampler, the smallest
#' number of onset times between changes
#' @param delay_model how the reporting delay varies over onset times:
#' "single" (the same throughout), "periods" (its own mean and sd in each of
#' the periods given by reporting_piece_index) or "changepoints" (it changes
#' at unknown onset times)
#' @param reporting_piece_index the period of each onset time, for a delay
#' model of "periods"
#' @param delay_max_changes,delay_min_length as Rt_max_changes and
#' Rt_min_length, for the delay
#' @param weekday_first the day of the week (from 1 to 7) of onset time 1, if
#' there is to be an effect of the day of the week of each report; NULL for
#' none
#' @param tail whether a proportion of cases have a delay which is
#' exponential rather than gamma
#' @param report_effects whether each report has its own random effect on the
#' rate of reporting
#' @param misdating NULL if cases are never misdated; or a named list of
#' 'shifts' and 'shift_probs' (see \code{unrecorded_prob})
#' @param batches NULL for no batches of erroneous records; or a named list
#' as for the 'batches' element of the priors of \code{mcmc_single} (exact
#' sampler only)
#' @param max_erroneous the largest number of erroneous records held for an
#' onset time at once; by default one more than the largest fall in a
#' reported count, and at least 3
#' @param days_sampled in the relaxed model, the number of the latest onset
#' times whose true cases are unknown
#'
#' @return a named list of options
#' @export
model_options <- function(is_negative_binomial=TRUE,
                          Rt_model=c("random_walk", "periods", "changepoints"),
                          Rt_period_width=7, Rt_index=NULL,
                          Rt_max_changes=3, Rt_min_length=7,
                          delay_model=c("single", "periods", "changepoints"),
                          reporting_piece_index=NULL,
                          delay_max_changes=2, delay_min_length=14,
                          weekday_first=NULL, tail=FALSE, report_effects=FALSE,
                          misdating=NULL, batches=NULL, max_erroneous=NULL,
                          days_sampled=35) {
  Rt_model <- match.arg(Rt_model)
  delay_model <- match.arg(delay_model)
  if(delay_model == "periods" & is.null(reporting_piece_index))
    stop("reporting_piece_index must be given if the delay differs between periods.")
  if(!is.null(misdating)) {
    if(is.null(misdating$shifts) | is.null(misdating$shift_probs))
      stop("misdating must be a named list of 'shifts' and 'shift_probs'.")
  }
  list(is_negative_binomial=is_negative_binomial,
       Rt_model=Rt_model, Rt_period_width=Rt_period_width, Rt_index=Rt_index,
       Rt_max_changes=Rt_max_changes, Rt_min_length=Rt_min_length,
       delay_model=delay_model, reporting_piece_index=reporting_piece_index,
       delay_max_changes=delay_max_changes, delay_min_length=delay_min_length,
       weekday_first=weekday_first, tail=tail, report_effects=report_effects,
       misdating=misdating, batches=batches, max_erroneous=max_erroneous,
       days_sampled=days_sampled)
}

#' Priors of the model, common to all ways of fitting it
#'
#' @param ... named priors to change from their defaults, which are those of
#' \code{relaxed_priors}
#'
#' @return a named list of priors
#' @export
model_priors <- function(...) {
  priors <- relaxed_priors()
  changes <- list(...)
  unknown <- setdiff(names(changes), names(priors))
  if(length(unknown) > 0)
    stop(paste0("Unknown priors: ", paste(unknown, collapse=", "), "."))
  priors[names(changes)] <- changes
  priors
}

#' Checks that the options are available for a way of fitting the model
#'
#' @param options as returned by \code{model_options}
#' @param engine as for \code{fit_incidence}
#'
#' @return nothing; stops with a message if an option is not available
check_options_engine <- function(options, engine) {
  if(engine == "mcmc") {
    if(options$Rt_model == "random_walk")
      stop(paste("A random walk for Rt is not available in the exact sampler:",
                 "use an Rt_model of 'periods' or 'changepoints', or another engine."))
  } else {
    if(!is.null(options$batches))
      stop("Batches of erroneous records are only available in the exact sampler (engine 'mcmc').")
    if(options$delay_model == "changepoints" & !is.null(options$misdating))
      stop(paste("In the relaxed model, changes in the delay at unknown times",
                 "cannot be combined with misdating."))
  }
  invisible(NULL)
}

#' The period of each onset time
#'
#' @param n_days number of onset times
#' @param width number of onset times in each period
#' @param index the period of each onset time, if given in place of a width
#'
#' @return a vector of periods, aligned to end on the last onset time
period_index <- function(n_days, width, index=NULL) {
  if(!is.null(index)) {
    if(length(index) != n_days)
      stop("There must be a period for each onset time.")
    return(as.integer(index))
  }
  offset <- (width - n_days %% width) %% width
  as.integer(ceiling((seq_len(n_days) + offset) / width))
}

#' Fits the model by the exact sampler, with arguments and result in the
#' common form
#'
#' @inheritParams fit_incidence
#' @param max_erroneous the largest number of erroneous records
#'
#' @return a named list of 'cases', 'Rt', 'delay' and 'parameters' in the
#' form described in \code{fit_incidence}
#' @importFrom rlang .data
fit_incidence_mcmc <- function(snapshot_df, serial_parameters, options, priors,
                               cases_before, serial_max, max_erroneous, control) {
  defaults <- list(niterations=1500, chains=3, draws=300, cores=1, seed=1,
                   n_parameter_steps=5,
                   reporting_metropolis_parameters=list(mean_step=0.3, sd_step=0.25),
                   revision_metropolis_parameters=list(eta_step=0.3, rho_step=0.3,
                                                       pi_step=0.3),
                   overdispersion_metropolis_sd=8, n_partition_moves=50)
  unknown <- setdiff(names(control), names(defaults))
  if(length(unknown) > 0)
    stop(paste0("Unknown settings for the exact sampler: ",
                paste(unknown, collapse=", "), "."))
  defaults[names(control)] <- control
  control <- defaults

  onset_times <- sort(unique(snapshot_df$time_onset))
  n_days <- length(onset_times)
  Rt_index <- rep(1L, n_days)
  if(options$Rt_model == "periods")
    Rt_index <- period_index(n_days, options$Rt_period_width, options$Rt_index)
  piece_index <- rep(1L, n_days)
  if(options$delay_model == "periods")
    piece_index <- period_index(n_days, 1, options$reporting_piece_index)
  data <- snapshot_df %>%
    dplyr::select("time_onset", "time_reported", "cases_reported") %>%
    dplyr::mutate(Rt_index=Rt_index[match(.data$time_onset, onset_times)],
                  reporting_piece_index=piece_index[match(.data$time_onset, onset_times)])

  priors_mcmc <- list(
    Rt=list(shape=priors$Rt_shape, rate=priors$Rt_rate),
    reporting=list(mean_mu=priors$delay_mean_mu, mean_sigma=priors$delay_mean_sigma,
                   sd_mu=priors$delay_sd_mu, sd_sigma=priors$delay_sd_sigma),
    overdispersion=list(mean=priors$kappa_mu, sd=priors$kappa_sigma),
    revisions=list(eta_mu=priors$eta_mu, eta_sigma=priors$eta_sigma,
                   rho_alpha=priors$rho_alpha, rho_beta=priors$rho_beta,
                   pi_alpha=priors$pi_alpha, pi_beta=priors$pi_beta),
    reporting_effects=list(
      tail_weight_alpha=priors$tail_weight_alpha, tail_weight_beta=priors$tail_weight_beta,
      tail_mean_mu=priors$tail_mean_mu, tail_mean_sigma=priors$tail_mean_sigma,
      weekday_sd=priors$weekday_sd, report_sd_shape=priors$report_sd_shape,
      report_sd_scale=priors$report_sd_scale),
    max_erroneous=max_erroneous)
  change_beta <- ifelse(is.na(priors$change_beta), n_days, priors$change_beta)
  if(options$Rt_model == "changepoints")
    priors_mcmc$Rt_changepoints <- list(alpha=priors$change_alpha, beta=change_beta,
                                        min_length=options$Rt_min_length)
  if(options$delay_model == "changepoints")
    priors_mcmc$reporting_changepoints <- list(
      alpha=priors$change_alpha, beta=change_beta,
      min_length=options$delay_min_length,
      n_partition_moves=control$n_partition_moves)
  if(!is.null(options$batches))
    priors_mcmc$batches <- options$batches

  initial_revisions <- list(eta=priors$eta_mu,
                            rho=priors$rho_alpha / (priors$rho_alpha + priors$rho_beta))
  if(!is.null(options$misdating)) {
    initial_revisions$pi <- priors$pi_alpha / (priors$pi_alpha + priors$pi_beta)
    initial_revisions$shifts <- options$misdating$shifts
    initial_revisions$shift_probs <- options$misdating$shift_probs
  }
  if(options$tail)
    initial_revisions$tail <- list(
      weight=priors$tail_weight_alpha / (priors$tail_weight_alpha + priors$tail_weight_beta),
      mean=priors$tail_mean_mu)
  if(!is.null(options$weekday_first)) {
    initial_revisions$weekday <- TRUE
    initial_revisions$weekday_first <- options$weekday_first
  }
  if(options$report_effects) {
    initial_revisions$report_effects <- TRUE
    if(is.null(initial_revisions$weekday_first))
      initial_revisions$weekday_first <- 1
  }

  largest_reported <- snapshot_df %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::summarise(cases_true=max(.data$cases_reported))
  days_before_end <- max(onset_times) - largest_reported$time_onset
  n_pieces <- max(piece_index)

  f_chain <- function(chain) {
    set.seed(control$seed + 1000 * chain)
    # starting values differ between chains
    inflation <- c(2, 1, 3.5)[(chain - 1) %% 3 + 1]
    initial_cases_true <- largest_reported %>%
      dplyr::mutate(cases_true=ceiling(
        .data$cases_true * (1 + inflation * exp(-days_before_end / 4)) +
          5 * (days_before_end < 5)))
    scale <- c(1, 0.6, 1.5)[(chain - 1) %% 3 + 1]
    initial_reporting <- dplyr::tibble(
      reporting_piece_index=seq_len(n_pieces),
      mean=priors$delay_mean_mu * scale, sd=priors$delay_sd_mu * scale)
    res <- mcmc_single(
      niterations=control$niterations,
      snapshot_with_Rt_index_df=data,
      priors=c(priors_mcmc,
               list(max_cases=ceiling(4 * max(initial_cases_true$cases_true) + 50))),
      serial_parameters=serial_parameters,
      initial_cases_true=initial_cases_true,
      initial_reporting_parameters=initial_reporting,
      initial_Rt=dplyr::tibble(Rt_index=seq_len(max(Rt_index)), Rt=1),
      reporting_metropolis_parameters=control$reporting_metropolis_parameters,
      serial_max=serial_max,
      print_to_screen=FALSE,
      initial_overdispersion=priors$kappa_mu * scale,
      is_negative_binomial=options$is_negative_binomial,
      overdispersion_metropolis_sd=control$overdispersion_metropolis_sd,
      initial_revision_parameters=initial_revisions,
      revision_metropolis_parameters=control$revision_metropolis_parameters,
      cases_before=cases_before,
      n_parameter_steps=control$n_parameter_steps)
    mcmc_common_draws(res, chain, control, onset_times, Rt_index, piece_index)
  }

  chains <- seq_len(control$chains)
  if(control$cores > 1)
    results <- parallel::mclapply(chains, f_chain, mc.cores=control$cores)
  else
    results <- lapply(chains, f_chain)
  failed <- purrr::map_lgl(results, ~inherits(., "try-error"))
  if(any(failed))
    stop(paste("A chain failed:", as.character(results[[which(failed)[1]]]))) # nocov
  stats::setNames(purrr::map(names(results[[1]]), function(name) {
    purrr::map_dfr(results, name)
  }), names(results[[1]]))
}

#' Draws from one chain of the exact sampler, in the common form
#'
#' The first half of the chain is discarded and the rest thinned to at most
#' control$draws draws.
#'
#' @param res as returned by \code{mcmc_single}
#' @param chain index of the chain
#' @param control settings of the sampler
#' @param onset_times the onset times
#' @param Rt_index,piece_index the periods of Rt and of the delay for each
#' onset time
#'
#' @return a named list of 'cases', 'Rt', 'delay' and 'parameters'
#' @importFrom rlang .data
mcmc_common_draws <- function(res, chain, control, onset_times, Rt_index,
                              piece_index) {
  warmup <- floor(control$niterations / 2)
  kept <- unique(round(seq(warmup + 1, control$niterations,
                           length.out=min(control$draws, control$niterations - warmup))))
  f_keep <- function(df) {
    df %>%
      dplyr::filter(.data$iteration %in% kept) %>%
      dplyr::mutate(chain=chain, draw=match(.data$iteration, kept)) %>%
      dplyr::select(-"iteration")
  }
  cases <- f_keep(res$cases) %>%
    dplyr::select("chain", "draw", "time_onset", "cases_true")
  if("time_onset" %in% colnames(res$Rt)) {
    Rt <- f_keep(res$Rt)
  } else {
    Rt <- f_keep(res$Rt) %>%
      dplyr::left_join(dplyr::tibble(time_onset=onset_times, Rt_index=Rt_index),
                       by="Rt_index", relationship="many-to-many")
  }
  Rt <- dplyr::select(Rt, "chain", "draw", "time_onset", "Rt")
  if(!is.null(res$reporting_days)) {
    delay <- f_keep(res$reporting_days)
  } else {
    delay <- f_keep(res$reporting) %>%
      dplyr::left_join(dplyr::tibble(time_onset=onset_times,
                                     reporting_piece_index=piece_index),
                       by="reporting_piece_index", relationship="many-to-many")
  }
  delay <- dplyr::select(delay, "chain", "draw", "time_onset", "mean", "sd")

  parameters <- f_keep(res$revisions) %>%
    dplyr::select("chain", "draw", "eta", "rho", pi_misdated="pi")
  if(all(parameters$pi_misdated == 0))
    parameters$pi_misdated <- NULL
  f_add <- function(df) {
    dplyr::left_join(parameters, dplyr::select(f_keep(df), -dplyr::any_of("iteration")),
                     by=c("chain", "draw"))
  }
  if(!is.null(res$overdispersion))
    parameters <- f_add(dplyr::rename(res$overdispersion, kappa="overdispersion"))
  if(!is.null(res$reporting_effects)) {
    effects <- res$reporting_effects
    parameters <- f_add(effects[, colSums(!is.na(effects)) > 0, drop=FALSE])
  }
  if(!is.null(res$Rt_changepoints))
    parameters <- f_add(dplyr::select(res$Rt_changepoints, "iteration",
                                      Rt_n_changes="n_changes", Rt_p_change="p_change"))
  if(!is.null(res$reporting_changepoints))
    parameters <- f_add(dplyr::select(res$reporting_changepoints, "iteration",
                                      delay_n_changes="n_changes"))
  if(!is.null(res$batch_parameters))
    parameters <- f_add(dplyr::select(res$batch_parameters, "iteration", "n_batches"))
  list(cases=cases, Rt=Rt, delay=delay, parameters=parameters)
}

#' Fits the relaxed model, with arguments and result in the common form
#'
#' @inheritParams fit_incidence_mcmc
#' @param engine one of the engines of \code{fit_incidence} other than "mcmc"
#'
#' @return a named list of 'cases', 'Rt', 'delay' and 'parameters', with the
#' Stan fit as 'raw'
#' @importFrom rlang .data
fit_incidence_relaxed <- function(snapshot_df, serial_parameters, options, priors,
                                  cases_before, serial_max, max_erroneous, engine,
                                  control) {
  data <- relaxed_stan_data(
    snapshot_df, serial_parameters, cases_before=cases_before,
    max_erroneous=max_erroneous, weekday_first=options$weekday_first,
    days_sampled=options$days_sampled, serial_max=serial_max,
    Rt_period_width=options$Rt_period_width, Rt_index=options$Rt_index,
    Rt_random_walk=options$Rt_model == "random_walk",
    reporting_piece_index=if(options$delay_model == "periods")
      options$reporting_piece_index else NULL,
    is_negative_binomial=options$is_negative_binomial, tail=options$tail,
    report_effects=options$report_effects, misdating=options$misdating,
    Rt_max_changes=ifelse(options$Rt_model == "changepoints", options$Rt_max_changes, 0),
    delay_max_changes=ifelse(options$delay_model == "changepoints",
                             options$delay_max_changes, 0),
    priors=priors)
  model <- control$model
  control$model <- NULL
  if(is.null(model))
    model <- relaxed_model()
  defaults <- if(engine == "nuts")
    list(chains=4, parallel_chains=4, iter_warmup=500, iter_sampling=500,
         adapt_delta=0.9, refresh=0)
  else
    list(draws=1000, refresh=0)
  defaults[names(control)] <- control
  fit <- do.call(fit_relaxed, c(list(data=data, method=engine, model=model), defaults))

  draws <- relaxed_draws(fit)
  chain_of_draw <- rep(1L, max(draws$cases$iteration))
  if(engine == "nuts") {
    chain_of_draw <- as.integer(posterior::as_draws_df(fit$draws("eta"))$.chain)
  }
  f_index <- function(df) {
    df %>%
      dplyr::mutate(chain=chain_of_draw[.data$iteration]) %>%
      dplyr::group_by(.data$chain) %>%
      dplyr::mutate(draw=match(.data$iteration, sort(unique(.data$iteration)))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-"iteration") %>%
      dplyr::relocate("chain", "draw")
  }
  parameters <- f_index(draws$parameters)
  # names in common with the exact sampler
  colnames(parameters) <- sub("^weekday_effect\\[(\\d)\\]$", "weekday_\\1",
                              colnames(parameters))
  if(!data$is_weekday)
    parameters <- dplyr::select(parameters, -dplyr::starts_with("weekday_"))
  if(!data$is_Rt_changepoints)
    parameters <- dplyr::select(parameters, -dplyr::any_of("Rt_n_changes"))
  if(!data$is_delay_changepoints)
    parameters <- dplyr::select(parameters, -dplyr::any_of("delay_n_changes"))
  parameters <- dplyr::select(parameters, -dplyr::starts_with("report_effect"),
                              -dplyr::matches("^delay_(mean|sd)"))
  list(cases=f_index(draws$cases), Rt=f_index(draws$Rt), delay=f_index(draws$delay),
       parameters=parameters, raw=fit)
}

#' Potential scale reduction factor (split R-hat) of a quantity
#'
#' Each chain is split in half, so that a drift within chains is detected as
#' well as a difference between them. Values near 1 indicate agreement.
#'
#' @param values draws of the quantity, in order within each chain
#' @param chain the chain of each draw
#'
#' @return the statistic, or NA if it cannot be calculated
#' @export
split_rhat <- function(values, chain) {
  halves <- unlist(lapply(split(values, chain), function(x) {
    n <- floor(length(x) / 2)
    list(x[seq_len(n)], x[n + seq_len(n)])
  }), recursive=FALSE)
  n <- min(lengths(halves))
  if(n < 2)
    return(NA_real_)
  halves <- lapply(halves, function(x) x[seq_len(n)])
  within <- mean(vapply(halves, stats::var, numeric(1)))
  between <- n * stats::var(vapply(halves, mean, numeric(1)))
  if(within == 0)
    return(ifelse(between == 0, 1, Inf))
  sqrt(((n - 1) / n * within + between / n) / within)
}

#' Diagnostics of a fit
#'
#' For engines which sample (the exact sampler and Hamiltonian Monte Carlo):
#' the split R-hat of each parameter, of Rt and the delay for the latest onset
#' time and of the total of the latest 14 onset times' cases; and the largest
#' over all onset times' cases. For approximations: the Pareto k statistic of
#' the importance ratios between the model and the approximation, where the
#' engine provides them (values above about 0.7 indicate that the
#' approximation is poor).
#'
#' @param result the draws in the common form
#' @param engine as for \code{fit_incidence}
#'
#' @return a named list of 'rhat' (a tibble, or NULL), 'rhat_cases_max',
#' 'pareto_k' and, for Hamiltonian Monte Carlo, 'divergences'
#' @importFrom rlang .data
fit_diagnostics <- function(result, engine) {
  out <- list(rhat=NULL, rhat_cases_max=NA_real_, pareto_k=NA_real_)
  if(engine %in% c("mcmc", "nuts")) {
    day_last <- max(result$Rt$time_onset)
    f_order <- function(df) dplyr::arrange(df, .data$chain, .data$draw)
    recent <- result$cases %>%
      dplyr::filter(.data$time_onset > day_last - 14) %>%
      dplyr::group_by(.data$chain, .data$draw) %>%
      dplyr::summarise(value=sum(.data$cases_true), .groups="drop") %>%
      f_order()
    Rt <- f_order(dplyr::filter(result$Rt, .data$time_onset == day_last))
    delay <- f_order(dplyr::filter(result$delay, .data$time_onset == day_last))
    parameters <- f_order(result$parameters)
    quantities <- c(
      list("Rt, latest onset time"=Rt$Rt,
           "Cases, latest 14 onset times"=recent$value,
           "Mean delay, latest onset time"=delay$mean,
           "Sd of delay, latest onset time"=delay$sd),
      as.list(dplyr::select(parameters, -"chain", -"draw")))
    chains <- c(list(Rt$chain, recent$chain, delay$chain, delay$chain),
                rep(list(parameters$chain), ncol(parameters) - 2))
    out$rhat <- dplyr::tibble(
      quantity=names(quantities),
      rhat=purrr::map2_dbl(quantities, chains, split_rhat))
    by_day <- result$cases %>%
      f_order() %>%
      dplyr::group_by(.data$time_onset) %>%
      dplyr::summarise(rhat=split_rhat(.data$cases_true, .data$chain))
    out$rhat_cases_max <- max(by_day$rhat[is.finite(by_day$rhat)], 1)
  }
  if(!is.null(result$raw)) {
    if(engine == "nuts")
      out$divergences <- sum(result$raw$diagnostic_summary(quiet=TRUE)$num_divergent)
    variables <- posterior::variables(result$raw$draws())
    if(all(c("lp__", "lp_approx__") %in% variables)) {
      lp <- posterior::as_draws_df(result$raw$draws(c("lp__", "lp_approx__")))
      ratios <- lp$lp__ - lp$lp_approx__
      out$pareto_k <- tryCatch(
        as.numeric(posterior::pareto_khat(ratios[is.finite(ratios)])$khat),
        error=function(e) NA_real_)
    }
  }
  out
}

#' Fits the model
#'
#' A single way of calling every method of fitting the model. The data, the
#' options of the model and the priors are given in the same form whichever
#' is used, and the result has the same form.
#'
#' The engines are:
#' \itemize{
#' \item "mcmc": the exact sampler (\code{mcmc_single}), in which true case
#' counts are whole numbers. Rt and the delay can change at unknown times
#' without a limit on the number of changes, and batches of erroneous records
#' are available. Rt cannot follow a random walk.
#' \item "nuts": Hamiltonian Monte Carlo on the relaxed model (see
#' \code{fit_relaxed}), in which true case counts are continuous. Changes at
#' unknown times are limited to a maximum number.
#' \item "laplace", "advi", "advi_fullrank", "pathfinder": approximations to
#' the posterior of the relaxed model, which take seconds but can be poor,
#' particularly where Rt or the delay changes at unknown times.
#' }
#' Options which an engine does not have give an error.
#'
#' @param snapshot_df a tibble with columns: time_onset, time_reported,
#' cases_reported. Onset times must be consecutive and start at 1.
#' @inheritParams state_process_logp
#' @param engine the method of fitting (see Details)
#' @param options as returned by \code{model_options}
#' @param priors as returned by \code{model_priors}
#' @inheritParams sample_cases_history
#' @param control a named list of settings for the engine. For "mcmc":
#' 'niterations', 'chains', 'draws' (kept from each chain, after its first
#' half), 'cores' (chains run at once), 'seed', 'n_parameter_steps',
#' 'n_partition_moves' and step sizes (see \code{mcmc_single}). For the
#' others: arguments of the corresponding cmdstanr method (such as 'chains',
#' 'iter_warmup', 'iter_sampling', 'draws' and 'seed'), and optionally
#' 'model', a model already compiled by \code{relaxed_model}.
#'
#' @return an object of class 'incidenceinflation_fit': a named list of
#' \itemize{
#' \item 'cases': a tibble of draws of true cases (chain, draw, time_onset,
#' cases_true)
#' \item 'Rt': a tibble of draws of Rt for each onset time (chain, draw,
#' time_onset, Rt)
#' \item 'delay': a tibble of draws of the mean and sd of the reporting delay
#' for each onset time (chain, draw, time_onset, mean, sd)
#' \item 'parameters': a tibble with a row for each draw and a column for each
#' other parameter (among kappa, eta, rho, pi_misdated, tail_weight,
#' tail_mean, report_sd, weekday_1 to weekday_7, Rt_n_changes and
#' delay_n_changes, as the options require)
#' \item 'diagnostics': see \code{fit_diagnostics}
#' \item 'engine', 'options', 'priors', 'seconds' and, for engines other than
#' "mcmc", 'raw' (the cmdstanr fit)
#' }
#' @export
#' @importFrom rlang .data
fit_incidence <- function(snapshot_df, serial_parameters,
                          engine=c("mcmc", "nuts", "laplace", "advi",
                                   "advi_fullrank", "pathfinder"),
                          options=model_options(), priors=model_priors(),
                          cases_before=NULL, serial_max=40, control=list()) {
  engine <- match.arg(engine)
  check_options_engine(options, engine)
  required <- c("time_onset", "time_reported", "cases_reported")
  if(!all(required %in% colnames(snapshot_df)))
    stop("snapshot_df must have columns: time_onset, time_reported, cases_reported.")
  snapshot_df <- snapshot_df[, required]
  onset_times <- sort(unique(snapshot_df$time_onset))
  if(!all(onset_times == seq_along(onset_times)))
    stop("Onset times must be consecutive and start at 1.")
  max_erroneous <- options$max_erroneous
  if(is.null(max_erroneous))
    max_erroneous <- max(3, max_downward_revision(snapshot_df) + 1)

  time_start <- Sys.time()
  if(engine == "mcmc")
    result <- fit_incidence_mcmc(snapshot_df, serial_parameters, options, priors,
                                 cases_before, serial_max, max_erroneous, control)
  else
    result <- fit_incidence_relaxed(snapshot_df, serial_parameters, options, priors,
                                    cases_before, serial_max, max_erroneous, engine,
                                    control)
  result$diagnostics <- fit_diagnostics(result, engine)
  result$engine <- engine
  result$options <- options
  result$priors <- priors
  result$seconds <- as.numeric(Sys.time() - time_start, units="secs")
  class(result) <- "incidenceinflation_fit"
  result
}

#' Summarises a fit
#'
#' @param object as returned by \code{fit_incidence}
#' @param probs the quantiles to give
#' @param ... not used
#'
#' @return a named list of tibbles of quantiles: 'cases', 'Rt' and 'delay'
#' (the mean of the delay) by onset time, and 'parameters'
#' @export
#' @importFrom rlang .data
summary.incidenceinflation_fit <- function(object, probs=c(0.05, 0.5, 0.95), ...) {
  f_quantiles <- function(values) {
    q <- stats::quantile(values, probs, na.rm=TRUE)
    stats::setNames(as.list(q), paste0("q", formatC(100 * probs, width=2, flag="0")))
  }
  f_by_day <- function(df, column) {
    df %>%
      dplyr::group_by(.data$time_onset) %>%
      dplyr::summarise(dplyr::as_tibble(f_quantiles(.data[[column]])), .groups="drop")
  }
  parameters <- object$parameters %>%
    dplyr::select(-"chain", -"draw") %>%
    tidyr::pivot_longer(dplyr::everything(), names_to="parameter") %>%
    dplyr::group_by(.data$parameter) %>%
    dplyr::summarise(dplyr::as_tibble(f_quantiles(.data$value)),
                     mean=mean(.data$value, na.rm=TRUE), .groups="drop")
  list(cases=f_by_day(object$cases, "cases_true"),
       Rt=f_by_day(object$Rt, "Rt"),
       delay=f_by_day(object$delay, "mean"),
       parameters=parameters)
}

#' Prints a fit
#'
#' @param x as returned by \code{fit_incidence}
#' @param ... not used
#'
#' @return x, invisibly
#' @export
print.incidenceinflation_fit <- function(x, ...) {
  n_draws <- nrow(x$parameters)
  cat("Fit by engine '", x$engine, "': ", n_draws, " draws in ",
      length(unique(x$parameters$chain)), " chain(s), ",
      round(x$seconds), " seconds\n", sep="")
  s <- summary(x)
  day_last <- max(s$Rt$time_onset)
  f_row <- function(row) sprintf("%.2f (%.2f, %.2f)", row$q50, row$q05, row$q95)
  cat("Rt, latest onset time: ", f_row(s$Rt[s$Rt$time_onset == day_last, ]), "\n", sep="")
  cat("Mean delay, latest onset time: ",
      f_row(s$delay[s$delay$time_onset == day_last, ]), "\n", sep="")
  d <- x$diagnostics
  if(!is.null(d$rhat))
    cat("Largest split R-hat: ", sprintf("%.2f", max(d$rhat$rhat, d$rhat_cases_max,
                                                    na.rm=TRUE)), "\n", sep="")
  if(!is.null(d$divergences))
    cat("Divergent transitions: ", d$divergences, "\n", sep="")
  if(!is.na(d$pareto_k))
    cat("Pareto k of the approximation: ", sprintf("%.2f", d$pareto_k), "\n", sep="")
  invisible(x)
}
