// A relaxed version of the model in which nothing is discrete, so that it can
// be fitted by Hamiltonian Monte Carlo or by variational methods:
// - true case counts of recent onset days are continuous (those of older onset
//   days are fixed at their latest reported counts);
// - Rt is constant within periods of onset days given in advance, and either
//   follows a random walk between them or is independent between them;
// - the reporting delay is the same within periods of onset days given in
//   advance.
// Optionally, Rt or the delay can instead change at unknown onset days, up to
// a maximum number of changes given in advance. The segment to which each
// onset day belongs is then a hidden state which can only stay as it is or
// move on by one, and the days of the changes are summed over exactly by a
// forward recursion. This has the prior on changes of the package (each day
// starts a new segment with probability p, which has a beta prior), limited
// to the maximum number, and without a minimum length of segment.
// The model of reporting is otherwise that of the package: a gamma delay,
// optionally with an exponential tail; optional effects of the day of the week
// of each report and of each report; erroneous records which are entered and
// later removed (summed over by a forward recursion); and optionally the
// misdating of cases to neighbouring onset days.
// Optionally the deaths among the reported cases are modelled too (see
// deaths_logp): the case fatality ratio and the delay from onset to a death
// being recorded are each constant within segments of onset days, with
// changes at unknown onset days up to a maximum number, summed over together.
functions {
  // binomial probability mass at k successes of n trials, for real n and k
  real binomial_real_logp(real k, real n, real p) {
    if (k < 0 || k > n) {
      return negative_infinity();
    }
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
           + lmultiply(k, p) + lmultiply(n - k, 1 - p);
  }

  // negative binomial density with mean mu and overdispersion kappa, for real y
  real neg_binomial_real_logp(real y, real mu, real kappa) {
    return lgamma(y + kappa) - lgamma(y + 1) - lgamma(kappa)
           + kappa * (log(kappa) - log(kappa + mu))
           + lmultiply(y, mu / (kappa + mu));
  }

  // Poisson density with mean mu, for real y
  real poisson_real_logp(real y, real mu) {
    return lmultiply(y, mu) - mu - lgamma(y + 1);
  }

  // log of the probability that a gamma delay exceeds d - 0.5 days, for
  // d = 1, ..., n_days. The built-in functions for this have gradients which
  // fail for some values, so the density is integrated numerically instead
  // (Simpson's rule with four steps a day), working down from a point well
  // beyond the last day, where the remaining tail is given by its asymptotic
  // expansion. Only the density at delays of half a day or more is used, which
  // is smooth even where the shape is less than one.
  vector gamma_log_survival(int n_days, real shape, real rate) {
    int n_extra = 60;                       // days integrated beyond the last
    int n_total = n_days - 1 + n_extra;     // whole days integrated over
    int n_points = 4 * n_total + 1;
    vector[n_points] x = linspaced_vector(n_points, 0.5, 0.5 + n_total);
    vector[n_points] density = exp(shape * log(rate) - lgamma(shape)
                                   + (shape - 1) * log(x) - rate * x);
    real x_end = x[n_points];
    vector[n_total + 1] survival;
    survival[n_total + 1] = density[n_points] / rate
                            * (1 + (shape - 1) / (rate * x_end)
                               + (shape - 1) * (shape - 2) / square(rate * x_end));
    for (i in 1:n_total) {
      int day = n_total + 1 - i;            // interval from x[4 day - 3] to x[4 day + 1]
      int k = 4 * day - 3;
      survival[day] = survival[day + 1]
                      + (density[k] + 4 * density[k + 1] + 2 * density[k + 2]
                         + 4 * density[k + 3] + density[k + 4]) / 12;
    }
    return log(survival[1:n_days]);
  }

  // log of the probability that the delay exceeds d - 0.5 days, d = 1, ...,
  // n_days, where a proportion tail_weight of cases have a delay which is
  // exponential with mean tail_mean rather than gamma
  vector delay_log_survival(int n_days, real delay_mean, real delay_sd,
                            real tail_weight, real tail_mean) {
    vector[n_days] out = gamma_log_survival(n_days, square(delay_mean / delay_sd),
                                            delay_mean / square(delay_sd));
    if (tail_weight > 0) {
      for (d in 1:n_days) {
        out[d] = log_sum_exp(log1m(tail_weight) + out[d],
                             log(tail_weight) - (d - 0.5) / tail_mean);
      }
    }
    return out;
  }

  // log-probability that a case remains unreported at the end of each of the
  // days from its onset day on (element k is day onset + k - 1), with the
  // daily hazard of reporting multiplied by that day's multiplier
  vector log_undetected_row(int day_onset, vector log_base, vector hazard_multiplier) {
    int n_days = rows(hazard_multiplier);
    vector[n_days] out = rep_vector(0, n_days);
    out[1] = log_base[1];
    for (k in (day_onset + 1):n_days) {
      out[k - day_onset + 1] = out[k - day_onset]
        - hazard_multiplier[k] * (log_base[k - day_onset] - log_base[k - day_onset + 1]);
    }
    return out;
  }

  // probability that a case has not yet been recorded under its own onset day
  // at the end of each day from its onset day on: it is either unreported, or
  // reported but misdated and not yet corrected
  vector unrecorded_row(vector log_undetected, int n_after, real rho,
                        real p_misdated, array[] int shifts, vector shift_probs) {
    vector[n_after + 1] out;
    real awaiting = 0;
    out[1] = exp(log_undetected[1]);
    for (k in 1:n_after) {
      real p_possible = 0;
      for (s in 1:size(shifts)) {
        if (shifts[s] <= k - 1) {
          p_possible += shift_probs[s];
        }
      }
      awaiting = awaiting * (1 - rho)
                 + p_possible * (exp(log_undetected[k]) - exp(log_undetected[k + 1]));
      out[k + 1] = exp(log_undetected[k + 1]) + p_misdated * awaiting;
    }
    return out;
  }

  // probability that a case of onset day day_true remains unreported at the
  // end of a day (one if the day is before its onset)
  real undetected_at(matrix log_undetected, int day_true, int day) {
    if (day < day_true) {
      return 1;
    }
    return exp(log_undetected[day_true, day - day_true + 1]);
  }

  // expected numbers of records of cases of neighbouring onset days which are
  // misdated to an onset day: element 1 is the number present at its first
  // report, and element j + 1 the number entered between reports j and j + 1
  vector misdated_entries(int day_onset, array[] int days, vector cases_true,
                          matrix log_undetected, real rho, real p_misdated,
                          array[] int shifts, vector shift_probs) {
    int n = size(days);
    int n_days = rows(cases_true);
    vector[n] out = rep_vector(0, n);
    for (s in 1:size(shifts)) {
      int day_true = day_onset - shifts[s];
      if (day_true >= 1 && day_true <= n_days) {
        real rate = p_misdated * shift_probs[s] * cases_true[day_true];
        for (k in (day_onset + 1):days[1]) {
          out[1] += rate
                    * (undetected_at(log_undetected, day_true, k - 1)
                       - undetected_at(log_undetected, day_true, k))
                    * pow(1 - rho, days[1] - k);
        }
        for (j in 2:n) {
          out[j] += rate * (undetected_at(log_undetected, day_true, days[j - 1])
                            - undetected_at(log_undetected, day_true, days[j]));
        }
      }
    }
    return out;
  }

  // log-probability of the reported counts of one onset day given its true
  // cases. log_undetected[k] is the log-probability that a case remains
  // unreported at the end of day (day_onset + k - 1), and unrecorded[k] the
  // probability that it has not been recorded under its onset day by then
  // (the same thing unless cases can be misdated). entries_extra holds the
  // expected numbers of misdated records from neighbouring onset days (see
  // misdated_entries), which are removed in the same way as erroneous records.
  real observation_logp(real cases_true, array[] int days, array[] int counts,
                        int day_onset, vector log_undetected, vector unrecorded,
                        vector entries_extra, real eta, real rho,
                        int max_erroneous) {
    int n = size(days);
    vector[max_erroneous + 1] log_alpha = rep_vector(negative_infinity(),
                                                     max_erroneous + 1);
    if (n <= 1) {
      return 0;
    }
    // erroneous records present at the first report, given its count
    if (days[1] == day_onset || (eta == 0 && entries_extra[1] == 0)) {
      log_alpha[1] = 0;
    } else {
      real mean_erroneous = 0;
      real p_detect = 1 - unrecorded[days[1] - day_onset + 1] / unrecorded[1];
      for (k in (day_onset + 1):days[1]) {
        mean_erroneous += (exp(log_undetected[k - day_onset])
                           - exp(log_undetected[k - day_onset + 1]))
                          * pow(1 - rho, days[1] - k);
      }
      mean_erroneous = eta * cases_true * mean_erroneous + entries_extra[1];
      for (e in 0:min(max_erroneous, counts[1])) {
        if (counts[1] - e <= cases_true) {
          log_alpha[e + 1] = poisson_lpmf(e | mean_erroneous)
                             + binomial_real_logp(counts[1] - e, cases_true, p_detect);
        }
      }
      log_alpha -= log_sum_exp(log_alpha);
    }

    for (j in 2:n) {
      int k_before = days[j - 1] - day_onset + 1;
      int k_after = days[j] - day_onset + 1;
      real p_reported = 1 - unrecorded[k_after] / unrecorded[k_before];
      real mean_entered = eta * cases_true * (exp(log_undetected[k_before])
                                              - exp(log_undetected[k_after]))
                          + entries_extra[j];
      real p_removed = 1 - pow(1 - rho, days[j] - days[j - 1]);
      int change = counts[j] - counts[j - 1];
      vector[max_erroneous + 1] log_alpha_new = rep_vector(negative_infinity(),
                                                           max_erroneous + 1);
      for (e in 0:max_erroneous) {
        if (log_alpha[e + 1] > negative_infinity() && counts[j - 1] - e >= 0) {
          real not_yet_reported = cases_true - (counts[j - 1] - e);
          for (removed in 0:e) {
            real lp_removed = binomial_lpmf(removed | e, p_removed);
            for (e_new in 0:max_erroneous) {
              int entered = e_new - e + removed;
              int reported = change - entered + removed;
              // impossible combinations are left out, since adding a
              // log-probability of minus infinity to another gives a
              // gradient which is not a number
              if (entered >= 0 && reported >= 0 && reported <= not_yet_reported) {
                log_alpha_new[e_new + 1] = log_sum_exp(
                  log_alpha_new[e_new + 1],
                  log_alpha[e + 1] + lp_removed
                  + poisson_lpmf(entered | mean_entered)
                  + binomial_real_logp(reported, not_yet_reported, p_reported));
              }
            }
          }
        }
      }
      log_alpha = log_alpha_new;
    }
    return log_sum_exp(log_alpha);
  }

  // Forward recursion over segments. ll[t, k] is the log-likelihood of day t
  // if it belongs to segment k. Each day after the first starts the next
  // segment with probability p_change. Every path with j changes has weight
  // p_change^j (1 - p_change)^(n - 1 - j), whichever segment it ends in.
  // Returns the matrix whose element (t, k) is the log of the summed weight of
  // all paths over days 1 to t which end in segment k.
  matrix changepoint_forward(matrix ll, real p_change) {
    int n = rows(ll);
    int n_segments = cols(ll);
    matrix[n, n_segments] log_alpha = rep_matrix(negative_infinity(), n, n_segments);
    log_alpha[1, 1] = ll[1, 1];
    for (t in 2:n) {
      for (k in 1:(t < n_segments ? t : n_segments)) {
        // terms of minus infinity are left out of sums, since their
        // gradients are not numbers
        int can_stay = log_alpha[t - 1, k] > negative_infinity();
        int can_move = k > 1 ? log_alpha[t - 1, k - 1] > negative_infinity() : 0;
        if (ll[t, k] > negative_infinity()) {
          if (can_stay && can_move) {
            log_alpha[t, k] = ll[t, k]
                              + log_sum_exp(log_alpha[t - 1, k] + log1m(p_change),
                                            log_alpha[t - 1, k - 1] + log(p_change));
          } else if (can_stay) {
            log_alpha[t, k] = ll[t, k] + log_alpha[t - 1, k] + log1m(p_change);
          } else if (can_move) {
            log_alpha[t, k] = ll[t, k] + log_alpha[t - 1, k - 1] + log(p_change);
          }
        }
      }
    }
    return log_alpha;
  }

  // log of the summed weight of all paths
  real changepoint_marginal(matrix ll, real p_change) {
    int n = rows(ll);
    int n_segments = cols(ll);
    int n_reachable = n < n_segments ? n : n_segments;
    row_vector[n_segments] last = changepoint_forward(ll, p_change)[n];
    real total = negative_infinity();
    for (k in 1:n_reachable) {
      if (last[k] > negative_infinity()) {
        total = total > negative_infinity() ? log_sum_exp(total, last[k]) : last[k];
      }
    }
    return total;
  }

  // a draw of the segment of each day given everything else
  array[] int changepoint_segments_rng(matrix ll, real p_change) {
    int n = rows(ll);
    int n_segments = cols(ll);
    matrix[n, n_segments] log_alpha = changepoint_forward(ll, p_change);
    array[n] int segment;
    int n_reachable = n < n_segments ? n : n_segments;
    segment[n] = categorical_logit_rng(log_alpha[n, 1:n_reachable]');
    for (i in 1:(n - 1)) {
      int t = n - i;
      int k = segment[t + 1];
      if (k == 1) {
        segment[t] = 1;
      } else if (log_alpha[t, k] == negative_infinity()) {
        segment[t] = k - 1;                  // segment k cannot yet have started
      } else {
        real stay = log_alpha[t, k] + log1m(p_change);
        real move = log_alpha[t, k - 1] + log(p_change);
        segment[t] = bernoulli_logit_rng(stay - move) ? k : k - 1;
      }
    }
    return segment;
  }

  // log-likelihood of the true cases of each onset day under the renewal
  // model, for each of a set of values of log Rt (columns)
  matrix renewal_ll(vector cases_true, vector cases_before, vector serial_weights,
                    vector log_Rt, int is_negative_binomial, real kappa) {
    int n_days = rows(cases_true);
    int serial_max = rows(serial_weights);
    int n_values = rows(log_Rt);
    vector[serial_max + n_days] cases_all = append_row(cases_before, cases_true);
    matrix[n_days, n_values] ll = rep_matrix(0, n_days, n_values);
    for (t in 1:n_days) {
      real infectiousness = dot_product(
        serial_weights, reverse(cases_all[t:(t + serial_max - 1)]));
      if (infectiousness > 0) {
        for (k in 1:n_values) {
          real mean_cases = exp(log_Rt[k]) * infectiousness;
          if (is_negative_binomial) {
            ll[t, k] = neg_binomial_real_logp(cases_true[t], mean_cases, kappa);
          } else {
            ll[t, k] = poisson_real_logp(cases_true[t], mean_cases);
          }
        }
      }
    }
    return ll;
  }

  // log-likelihood of the reported counts of each onset day. If
  // is_delay_changepoints there is a column for each segment of the delay;
  // otherwise there is a single column, each onset day having the delay of its
  // period.
  matrix reporting_ll(vector cases_true, array[] int obs_day, array[] int obs_count,
                      array[] int obs_start, array[] int obs_n,
                      vector delay_mean, vector delay_sd, array[] int delay_period,
                      int is_delay_changepoints, real tail_weight, real tail_mean,
                      vector hazard_multiplier, real eta, real rho, int max_erroneous,
                      real p_misdated, array[] int shifts, vector shift_probs) {
    int n_days = rows(cases_true);
    int n_delays = rows(delay_mean);
    int n_columns = is_delay_changepoints ? n_delays : 1;
    matrix[n_days, n_columns] ll = rep_matrix(0, n_days, n_columns);
    array[n_delays] vector[n_days] log_base;
    matrix[n_days, n_days] log_undetected;
    for (p in 1:n_delays) {
      log_base[p] = delay_log_survival(n_days, delay_mean[p], delay_sd[p],
                                       tail_weight, tail_mean);
    }
    if (!is_delay_changepoints) {
      for (t in 1:n_days) {
        log_undetected[t] = log_undetected_row(t, log_base[delay_period[t]],
                                               hazard_multiplier)';
      }
    }
    for (t in 1:n_days) {
      if (obs_n[t] > 1) {
        array[obs_n[t]] int days = obs_day[obs_start[t]:(obs_start[t] + obs_n[t] - 1)];
        array[obs_n[t]] int counts = obs_count[obs_start[t]:(obs_start[t] + obs_n[t] - 1)];
        int n_after = n_days - t;
        for (c in 1:n_columns) {
          vector[n_after + 1] row_undetected;
          vector[n_after + 1] unrecorded;
          vector[obs_n[t]] entries_extra = rep_vector(0, obs_n[t]);
          if (is_delay_changepoints) {
            row_undetected = log_undetected_row(t, log_base[c],
                                                hazard_multiplier)[1:(n_after + 1)];
          } else {
            row_undetected = log_undetected[t, 1:(n_after + 1)]';
          }
          unrecorded = exp(row_undetected);
          if (p_misdated > 0) {
            unrecorded = unrecorded_row(row_undetected, n_after, rho, p_misdated,
                                        shifts, shift_probs);
            entries_extra = misdated_entries(t, days, cases_true, log_undetected, rho,
                                             p_misdated, shifts, shift_probs);
          }
          ll[t, c] = observation_logp(cases_true[t], days, counts, t, row_undetected,
                                      unrecorded, entries_extra, eta, rho,
                                      max_erroneous);
        }
      }
    }
    return ll;
  }

  // Log-probability of the deaths shown in each report of one onset day,
  // given the cases in each report.
  //
  // Each case is in one of three states which the reports show: not yet
  // reported, shown as alive, or shown as dead. If fatal cases are reported
  // as others are, and the delay to a death being recorded is independent of
  // the reporting delay, the chance that a case moves between these states
  // from one report to the next does not depend on its earlier history. The
  // counts in each report then depend only on those in the report before,
  // and neither the number of fatal cases nor how many of the cases shown as
  // alive will die needs to be known:
  // - a case newly reported is shown as dead with probability
  //   cfr * (1 - S(t_j)), where S is the probability a death has not been
  //   recorded;
  // - a case shown as alive at the report before is shown as dead with
  //   probability cfr * (S(t_j-1) - S(t_j)) / (1 - cfr + cfr * S(t_j-1)).
  // This is exact when there are no wrong records.
  //
  // Wrong records shown as dead are entered with the expected numbers given
  // (the first being those present at the first report) and each is removed
  // with probability rho a day. The number present is summed over by a
  // forward recursion, up to max_dead at once. In this the cases newly
  // reported and those shown as alive are taken from the reported counts,
  // which leaves out the few wrong records among them.
  //
  // log_unrecorded[k] is the log-probability that a death has not been
  // recorded by the end of day (day_onset + k - 1).
  real deaths_logp(array[] int days, array[] int cases, array[] int deaths,
                   int day_onset, real cfr, vector log_unrecorded,
                   vector mean_entered, real rho, int max_dead) {
    int n = size(days);
    vector[max_dead + 1] log_alpha = rep_vector(negative_infinity(), max_dead + 1);
    for (j in 1:n) {
      real unrecorded = exp(log_unrecorded[days[j] - day_onset + 1]);
      real unrecorded_before = j > 1 ? exp(log_unrecorded[days[j - 1] - day_onset + 1]) : 1;
      real p_new = cfr * (1 - unrecorded);
      // the difference can be very slightly negative through rounding
      real p_alive = cfr * fmax(unrecorded_before - unrecorded, 0)
                     / (1 - cfr + cfr * unrecorded_before);
      int cases_before = j > 1 ? cases[j - 1] : 0;
      int deaths_before = j > 1 ? deaths[j - 1] : 0;
      int cases_new = cases[j] > cases_before ? cases[j] - cases_before : 0;
      real p_removed = j > 1 ? 1 - pow(1 - rho, days[j] - days[j - 1]) : 0;
      vector[max_dead + 1] log_alpha_new = rep_vector(negative_infinity(), max_dead + 1);
      for (e in 0:max_dead) {
        // before the first report there are no wrong records
        int is_start = j == 1 && e == 0;
        if ((is_start || (j > 1 && log_alpha[e + 1] > negative_infinity()))
            && deaths_before - e >= 0) {
          real lp_from = is_start ? 0 : log_alpha[e + 1];
          int alive = cases_before - (deaths_before - e);
          for (removed in 0:e) {
            real lp_removed = e > 0 ? binomial_lpmf(removed | e, p_removed) : 0;
            for (e_new in 0:max_dead) {
              int entered = e_new - e + removed;
              int deaths_new = (deaths[j] - e_new) - (deaths_before - e);
              if (entered >= 0 && deaths[j] - e_new >= 0 && deaths_new >= 0
                  && deaths_new <= cases_new + alive && alive >= 0
                  && (entered == 0 || mean_entered[j] > 0)) {
                real lp_deaths = negative_infinity();
                for (x in max(0, deaths_new - alive):min(cases_new, deaths_new)) {
                  real term = binomial_lpmf(x | cases_new, p_new)
                              + binomial_lpmf(deaths_new - x | alive, p_alive);
                  lp_deaths = lp_deaths > negative_infinity()
                              ? log_sum_exp(lp_deaths, term) : term;
                }
                if (lp_deaths > negative_infinity()) {
                  real term = lp_from + lp_removed + lp_deaths
                              + (mean_entered[j] > 0
                                 ? poisson_lpmf(entered | mean_entered[j]) : 0);
                  log_alpha_new[e_new + 1] =
                    log_alpha_new[e_new + 1] > negative_infinity()
                    ? log_sum_exp(log_alpha_new[e_new + 1], term) : term;
                }
              }
            }
          }
        }
      }
      log_alpha = log_alpha_new;
    }
    {
      real total = negative_infinity();
      for (e in 0:max_dead) {
        if (log_alpha[e + 1] > negative_infinity()) {
          total = total > negative_infinity() ? log_sum_exp(total, log_alpha[e + 1])
                                              : log_alpha[e + 1];
        }
      }
      return total;
    }
  }

  // probability that each onset day belongs to each segment, given everything
  // else (see changepoint_forward)
  matrix changepoint_smoothed(matrix ll, real p_change) {
    int n = rows(ll);
    int n_segments = cols(ll);
    matrix[n, n_segments] log_alpha = changepoint_forward(ll, p_change);
    matrix[n, n_segments] log_beta = rep_matrix(0, n, n_segments);
    matrix[n, n_segments] out = rep_matrix(0, n, n_segments);
    for (i in 1:(n - 1)) {
      int t = n - i;
      for (k in 1:n_segments) {
        real stay = log1m(p_change) + ll[t + 1, k] + log_beta[t + 1, k];
        log_beta[t, k] = stay;
        if (k < n_segments) {
          real move = log(p_change) + ll[t + 1, k + 1] + log_beta[t + 1, k + 1];
          if (move > negative_infinity()) {
            log_beta[t, k] = stay > negative_infinity() ? log_sum_exp(stay, move) : move;
          }
        }
      }
    }
    for (t in 1:n) {
      real total = negative_infinity();
      for (k in 1:n_segments) {
        real lp = log_alpha[t, k] + log_beta[t, k];
        if (log_alpha[t, k] > negative_infinity() && log_beta[t, k] > negative_infinity()) {
          total = total > negative_infinity() ? log_sum_exp(total, lp) : lp;
        }
      }
      for (k in 1:n_segments) {
        if (log_alpha[t, k] > negative_infinity() && log_beta[t, k] > negative_infinity()) {
          out[t, k] = exp(log_alpha[t, k] + log_beta[t, k] - total);
        }
      }
    }
    return out;
  }

  // Forward recursion over the segments of two quantities which change
  // independently. ll[t, (a - 1) * n_b + b] is the log-likelihood of day t if
  // it belongs to segment a of the first and segment b of the second. Returns
  // the logs of the summed weights of all paths over days 1 to t which end in
  // each pair of segments.
  matrix changepoint2_forward(matrix ll, int n_a, int n_b, real p_a, real p_b) {
    int n = rows(ll);
    matrix[n, n_a * n_b] log_alpha = rep_matrix(negative_infinity(), n, n_a * n_b);
    log_alpha[1, 1] = ll[1, 1];
    for (t in 2:n) {
      for (a in 1:n_a) {
        for (b in 1:n_b) {
          int s = (a - 1) * n_b + b;
          real total = negative_infinity();
          if (ll[t, s] > negative_infinity()) {
            for (move_a in 0:1) {
              for (move_b in 0:1) {
                if (a - move_a >= 1 && b - move_b >= 1
                    && (move_a == 0 || n_a > 1) && (move_b == 0 || n_b > 1)) {
                  real from = log_alpha[t - 1, (a - move_a - 1) * n_b + b - move_b];
                  if (from > negative_infinity()) {
                    real term = from
                                + (n_a > 1 ? (move_a ? log(p_a) : log1m(p_a)) : 0)
                                + (n_b > 1 ? (move_b ? log(p_b) : log1m(p_b)) : 0);
                    total = total > negative_infinity() ? log_sum_exp(total, term) : term;
                  }
                }
              }
            }
            if (total > negative_infinity()) {
              log_alpha[t, s] = ll[t, s] + total;
            }
          }
        }
      }
    }
    return log_alpha;
  }

  // log of the summed weight of all paths of the two quantities
  real changepoint2_marginal(matrix ll, int n_a, int n_b, real p_a, real p_b) {
    int n = rows(ll);
    row_vector[n_a * n_b] last = changepoint2_forward(ll, n_a, n_b, p_a, p_b)[n];
    real total = negative_infinity();
    for (s in 1:(n_a * n_b)) {
      if (last[s] > negative_infinity()) {
        total = total > negative_infinity() ? log_sum_exp(total, last[s]) : last[s];
      }
    }
    return total;
  }

  // logits in which impossible choices have a very large negative value, as
  // the random number functions do not accept minus infinity
  vector finite_logits(vector logits) {
    vector[rows(logits)] out = logits;
    for (i in 1:rows(logits)) {
      if (logits[i] == negative_infinity()) {
        out[i] = -1e30;
      }
    }
    return out;
  }

  // a draw of the pair of segments of each day given everything else: column
  // 1 is the segment of the first quantity and column 2 that of the second
  array[,] int changepoint2_segments_rng(matrix ll, int n_a, int n_b, real p_a,
                                         real p_b) {
    int n = rows(ll);
    matrix[n, n_a * n_b] log_alpha = changepoint2_forward(ll, n_a, n_b, p_a, p_b);
    array[n, 2] int segment;
    int s = categorical_logit_rng(finite_logits(log_alpha[n]'));
    segment[n, 1] = (s - 1) %/% n_b + 1;
    segment[n, 2] = (s - 1) % n_b + 1;
    for (i in 1:(n - 1)) {
      int t = n - i;
      int a = segment[t + 1, 1];
      int b = segment[t + 1, 2];
      vector[4] lp = rep_vector(negative_infinity(), 4);
      for (move_a in 0:1) {
        for (move_b in 0:1) {
          if (a - move_a >= 1 && b - move_b >= 1
              && (move_a == 0 || n_a > 1) && (move_b == 0 || n_b > 1)) {
            lp[1 + 2 * move_a + move_b] =
              log_alpha[t, (a - move_a - 1) * n_b + b - move_b]
              + (n_a > 1 ? (move_a ? log(p_a) : log1m(p_a)) : 0)
              + (n_b > 1 ? (move_b ? log(p_b) : log1m(p_b)) : 0);
          }
        }
      }
      {
        int choice = categorical_logit_rng(finite_logits(lp));
        segment[t, 1] = a - (choice - 1) %/% 2;
        segment[t, 2] = b - (choice - 1) % 2;
      }
    }
    return segment;
  }

  // Log-likelihood of the deaths reported for each onset day, with a column
  // for each pair of a segment of the case fatality ratio and a segment of
  // the delay to a death being recorded: column (a - 1) * n_recording + b.
  // unit_entered[i] is, for row i of the reports, the expected number of
  // wrong records entered since the report before for each case and each
  // unit of eta (see wrong_entry_units).
  matrix deaths_ll(vector cases_true, array[] int obs_day, array[] int obs_count,
                   array[] int obs_deaths, array[] int obs_start, array[] int obs_n,
                   array[] int max_dead, vector cfr, vector recording_mean,
                   vector recording_sd, vector hazard_multiplier,
                   vector unit_entered, real eta, real wrong_dead, real rho) {
    int n_days = rows(cases_true);
    int n_cfr = rows(cfr);
    int n_recording = rows(recording_mean);
    matrix[n_days, n_cfr * n_recording] ll = rep_matrix(0, n_days, n_cfr * n_recording);
    array[n_recording] vector[n_days] log_base;
    for (b in 1:n_recording) {
      log_base[b] = delay_log_survival(n_days, recording_mean[b], recording_sd[b], 0, 1);
    }
    for (t in 1:n_days) {
      if (obs_n[t] > 0) {
        int first = obs_start[t];
        int last = obs_start[t] + obs_n[t] - 1;
        vector[obs_n[t]] mean_entered = eta * wrong_dead * cases_true[t]
                                        * unit_entered[first:last];
        for (b in 1:n_recording) {
          vector[n_days - t + 1] log_unrecorded =
            log_undetected_row(t, log_base[b], hazard_multiplier)[1:(n_days - t + 1)];
          for (a in 1:n_cfr) {
            ll[t, (a - 1) * n_recording + b] = deaths_logp(
              obs_day[first:last], obs_count[first:last], obs_deaths[first:last],
              t, cfr[a], log_unrecorded, mean_entered, rho, max_dead[t]);
          }
        }
      }
    }
    return ll;
  }

  // For each row of the reports, the expected number of wrong records entered
  // since the report before, for each case of the onset day and each unit of
  // eta; for the first report of an onset day, of those entered since the
  // onset day which are still present. Where the delay changes at unknown
  // onset days, the delays of its segments are weighted by the probability
  // that the onset day belongs to each (weights, a row for each onset day).
  vector wrong_entry_units(int n_days, array[] int obs_day, array[] int obs_start,
                           array[] int obs_n, vector delay_mean, vector delay_sd,
                           matrix weights, real tail_weight, real tail_mean,
                           vector hazard_multiplier, real rho) {
    int n_delays = rows(delay_mean);
    vector[size(obs_day)] out = rep_vector(0, size(obs_day));
    array[n_delays] vector[n_days] log_base;
    for (p in 1:n_delays) {
      log_base[p] = delay_log_survival(n_days, delay_mean[p], delay_sd[p],
                                       tail_weight, tail_mean);
    }
    for (t in 1:n_days) {
      if (obs_n[t] > 0) {
        vector[n_days - t + 1] undetected = rep_vector(0, n_days - t + 1);
        for (p in 1:n_delays) {
          if (weights[t, p] > 0) {
            undetected += weights[t, p]
                          * exp(log_undetected_row(t, log_base[p],
                                                   hazard_multiplier)[1:(n_days - t + 1)]);
          }
        }
        for (i in 1:obs_n[t]) {
          int row = obs_start[t] + i - 1;
          if (i == 1) {
            for (k in (t + 1):obs_day[row]) {
              out[row] += (undetected[k - t] - undetected[k - t + 1])
                          * pow(1 - rho, obs_day[row] - k);
            }
          } else {
            out[row] = undetected[obs_day[row - 1] - t + 1]
                       - undetected[obs_day[row] - t + 1];
          }
        }
      }
    }
    return out;
  }

  // multiplier on the hazard of reporting on each day
  vector hazard_multipliers(int n_days, int is_weekday, array[] int weekday_report,
                            vector log_weekday, int is_report_effects,
                            array[] int report_index, vector log_report) {
    vector[n_days] out;
    for (k in 1:n_days) {
      real log_multiplier = 0;
      if (is_weekday && weekday_report[k] > 0) {
        log_multiplier += log_weekday[weekday_report[k]];
      }
      if (is_report_effects && report_index[k] > 0) {
        log_multiplier += log_report[report_index[k]];
      }
      out[k] = exp(log_multiplier);
    }
    return out;
  }
}
data {
  int<lower=1> n_days;                      // onset days
  int<lower=1> n_obs;                       // reported counts, by onset day then report
  array[n_obs] int<lower=1> obs_day;        // day of report
  array[n_obs] int<lower=0> obs_count;
  array[n_days] int<lower=0> obs_start;     // first row of each onset day (0 if none)
  array[n_days] int<lower=0> obs_n;         // number of rows of each onset day
  int<lower=0> max_erroneous;

  int<lower=0> n_sampled;                   // onset days whose true cases are unknown
  array[n_sampled] int<lower=1, upper=n_days> sampled_day;
  vector<lower=0>[n_days] cases_lower;      // lower bound of true cases (sampled days)
  vector<lower=0>[n_days] cases_fixed;      // true cases of the other days

  int<lower=1> serial_max;
  vector<lower=0>[serial_max] serial_weights;
  vector<lower=0>[serial_max] cases_before; // known cases before day 1, earliest first
  int<lower=0, upper=1> is_negative_binomial; // renewal model (Poisson if 0)

  int<lower=1> n_periods;                   // periods within which Rt is constant
  array[n_days] int<lower=1, upper=n_periods> period;
  int<lower=0, upper=1> is_Rt_random_walk;  // else Rt independent between periods
  int<lower=0, upper=1> is_Rt_changepoints; // Rt changes at unknown onset days: the
                                            // periods are then its segments in order
                                            // (n_periods is one more than the maximum
                                            // number of changes) and 'period' is unused

  int<lower=1> n_delay_periods;             // periods within which the delay is the same
  array[n_days] int<lower=1, upper=n_delay_periods> delay_period;
  int<lower=0, upper=1> is_delay_changepoints; // the delay changes at unknown onset
                                            // days: the delay periods are then its
                                            // segments and 'delay_period' is unused
  int<lower=0, upper=1> is_tail;            // exponential tail to the delay

  int<lower=0, upper=1> is_weekday;         // effect of the day of the week of a report
  array[n_days] int<lower=0, upper=7> weekday_report; // day of week of the first report
                                            // on or after each day (0 if none)
  int<lower=0, upper=1> is_report_effects;  // random effect of each report
  int<lower=0> n_reports;
  array[n_days] int<lower=0, upper=n_reports> report_index; // first report on or after
                                            // each day (0 if none)

  int<lower=0, upper=1> is_misdating;
  int<lower=0> n_shifts;
  array[n_shifts] int shifts;               // differences between recorded and true onset
  vector<lower=0>[n_shifts] shift_probs;

  // priors
  real<lower=0> delay_mean_mu;
  real<lower=0> delay_mean_sigma;
  real<lower=0> delay_sd_mu;
  real<lower=0> delay_sd_sigma;
  real<lower=0> eta_mu;
  real<lower=0> eta_sigma;
  real<lower=0> rho_alpha;
  real<lower=0> rho_beta;
  real<lower=0> kappa_mu;
  real<lower=0> kappa_sigma;
  real<lower=0> weekday_sd;
  real<lower=0> log_Rt_first_sd;            // sd of log Rt in the first period (random walk)
  real<lower=0> Rt_step_sd;                 // scale of the prior on the sd of steps in log Rt
  real<lower=0> Rt_shape;                   // gamma prior on Rt (independent periods)
  real<lower=0> Rt_rate;
  real<lower=0> tail_weight_alpha;
  real<lower=0> tail_weight_beta;
  real<lower=0> tail_mean_mu;
  real<lower=0> tail_mean_sigma;
  real<lower=0> report_sd_shape;            // inverse-gamma prior on the variance of the
  real<lower=0> report_sd_scale;            // logs of the effects of reports
  real<lower=0> pi_alpha;
  real<lower=0> pi_beta;
  real<lower=0> change_alpha;               // beta prior on the probability that each
  real<lower=0> change_beta;                // onset day starts a new segment

  // deaths among the reported cases
  int<lower=0, upper=1> is_deaths;
  array[is_deaths ? n_obs : 0] int<lower=0> obs_deaths; // of obs_count, shown as dead
  int<lower=1> n_cfr_segments;              // one more than the maximum number of
  int<lower=1> n_recording_segments;        // changes in each (1 if no deaths)
  array[n_days] int<lower=0> max_dead;      // most wrong records shown as dead at once
  int<lower=0, upper=1> is_recording_weekday; // the delay to a death being recorded has
                                            // its own effects of the day of the week
  real<lower=0> cfr_alpha;                  // beta prior on the case fatality ratio
  real<lower=0> cfr_beta;
  real<lower=0> recording_mean_mu;          // gamma priors on the mean and sd of the
  real<lower=0> recording_mean_sigma;       // delay from onset to a death being recorded
  real<lower=0> recording_sd_mu;
  real<lower=0> recording_sd_sigma;
  real<lower=0> wrong_dead_alpha;           // beta prior on the probability a wrong
  real<lower=0> wrong_dead_beta;            // record is shown as dead
}
parameters {
  vector<lower=0>[n_sampled] cases_extra;   // true cases above their lower bound
  vector<lower=0>[n_delay_periods] delay_mean;
  vector<lower=0>[n_delay_periods] delay_sd;
  sum_to_zero_vector[7] log_weekday;
  real<lower=0> eta;
  real<lower=0, upper=1> rho;
  array[is_negative_binomial] real<lower=0> kappa;
  vector[n_periods] log_Rt_raw;
  array[is_Rt_random_walk] real<lower=0> Rt_step;
  array[is_tail] real<lower=0, upper=1> tail_weight;
  array[is_tail] real<lower=0> tail_mean;
  vector[is_report_effects ? n_reports : 0] report_raw;
  array[is_report_effects] real<lower=0> report_sd;
  array[is_misdating] real<lower=0, upper=1> pi_misdated;
  array[is_Rt_changepoints] real<lower=0, upper=1> Rt_p_change;
  array[is_delay_changepoints] real<lower=0, upper=1> delay_p_change;
  vector<lower=0, upper=1>[is_deaths ? n_cfr_segments : 0] cfr;
  vector<lower=0>[is_deaths ? n_recording_segments : 0] recording_mean;
  vector<lower=0>[is_deaths ? n_recording_segments : 0] recording_sd;
  sum_to_zero_vector[7] log_weekday_recording;
  array[is_deaths] real<lower=0, upper=1> wrong_dead;
  array[is_deaths && n_cfr_segments > 1] real<lower=0, upper=1> cfr_p_change;
  array[is_deaths && n_recording_segments > 1] real<lower=0, upper=1> recording_p_change;
}
transformed parameters {
  vector[n_days] cases_true = cases_fixed;
  vector[n_periods] log_Rt;
  cases_true[sampled_day] = cases_lower[sampled_day] + cases_extra;
  if (is_Rt_random_walk) {
    log_Rt[1] = log_Rt_first_sd * log_Rt_raw[1];
    for (p in 2:n_periods) {
      log_Rt[p] = log_Rt[p - 1] + Rt_step[1] * log_Rt_raw[p];
    }
  } else {
    log_Rt = log_Rt_raw;
  }
}
model {
  vector[n_days] hazard_multiplier = hazard_multipliers(
    n_days, is_weekday, weekday_report, log_weekday, is_report_effects, report_index,
    is_report_effects ? report_sd[1] * report_raw : rep_vector(0, 0));
  real p_misdated = is_misdating ? pi_misdated[1] : 0;

  // priors (gamma distributions given by mean and sd)
  delay_mean ~ gamma(square(delay_mean_mu / delay_mean_sigma),
                     delay_mean_mu / square(delay_mean_sigma));
  delay_sd ~ gamma(square(delay_sd_mu / delay_sd_sigma),
                   delay_sd_mu / square(delay_sd_sigma));
  eta ~ gamma(square(eta_mu / eta_sigma), eta_mu / square(eta_sigma));
  rho ~ beta(rho_alpha, rho_beta);
  log_weekday ~ normal(0, weekday_sd);
  if (is_negative_binomial) {
    kappa ~ gamma(square(kappa_mu / kappa_sigma), kappa_mu / square(kappa_sigma));
  }
  if (is_Rt_random_walk) {
    log_Rt_raw ~ std_normal();
    Rt_step ~ normal(0, Rt_step_sd);
  } else {
    // gamma prior on Rt, with the Jacobian of taking logs
    target += gamma_lpdf(exp(log_Rt_raw) | Rt_shape, Rt_rate) + sum(log_Rt_raw);
  }
  if (is_tail) {
    tail_weight ~ beta(tail_weight_alpha, tail_weight_beta);
    tail_mean ~ gamma(square(tail_mean_mu / tail_mean_sigma),
                      tail_mean_mu / square(tail_mean_sigma));
  }
  if (is_report_effects) {
    report_raw ~ std_normal();
    // inverse-gamma prior on the variance, with the Jacobian of squaring
    target += inv_gamma_lpdf(square(report_sd[1]) | report_sd_shape, report_sd_scale)
              + log(report_sd[1]);
  }
  if (is_misdating) {
    pi_misdated ~ beta(pi_alpha, pi_beta);
  }
  if (is_Rt_changepoints) {
    Rt_p_change ~ beta(change_alpha, change_beta);
  }
  if (is_delay_changepoints) {
    delay_p_change ~ beta(change_alpha, change_beta);
  }
  log_weekday_recording ~ normal(0, weekday_sd);
  if (is_deaths) {
    cfr ~ beta(cfr_alpha, cfr_beta);
    recording_mean ~ gamma(square(recording_mean_mu / recording_mean_sigma),
                           recording_mean_mu / square(recording_mean_sigma));
    recording_sd ~ gamma(square(recording_sd_mu / recording_sd_sigma),
                         recording_sd_mu / square(recording_sd_sigma));
    wrong_dead ~ beta(wrong_dead_alpha, wrong_dead_beta);
    cfr_p_change ~ beta(change_alpha, change_beta);
    recording_p_change ~ beta(change_alpha, change_beta);
  }

  // renewal model
  {
    matrix[n_days, n_periods] ll = renewal_ll(
      cases_true, cases_before, serial_weights, log_Rt, is_negative_binomial,
      is_negative_binomial ? kappa[1] : 1);
    if (is_Rt_changepoints) {
      target += changepoint_marginal(ll, Rt_p_change[1]);
    } else {
      for (t in 1:n_days) {
        target += ll[t, period[t]];
      }
    }
  }

  // reporting
  {
    matrix[n_days, is_delay_changepoints ? n_delay_periods : 1] ll = reporting_ll(
      cases_true, obs_day, obs_count, obs_start, obs_n, delay_mean, delay_sd,
      delay_period, is_delay_changepoints, is_tail ? tail_weight[1] : 0,
      is_tail ? tail_mean[1] : 1, hazard_multiplier, eta, rho, max_erroneous,
      p_misdated, shifts, shift_probs);
    if (is_delay_changepoints) {
      target += changepoint_marginal(ll, delay_p_change[1]);
    } else {
      target += sum(ll);
    }

    // deaths among the reported cases
    if (is_deaths) {
      matrix[n_days, n_delay_periods] weights = rep_matrix(0, n_days, n_delay_periods);
      if (is_delay_changepoints) {
        weights = changepoint_smoothed(ll, delay_p_change[1]);
      } else {
        for (t in 1:n_days) {
          weights[t, delay_period[t]] = 1;
        }
      }
      target += changepoint2_marginal(
        deaths_ll(cases_true, obs_day, obs_count, obs_deaths, obs_start, obs_n,
                  max_dead, cfr, recording_mean, recording_sd,
                  hazard_multipliers(n_days, is_recording_weekday, weekday_report,
                                     log_weekday_recording, 0, report_index,
                                     rep_vector(0, 0)),
                  wrong_entry_units(n_days, obs_day, obs_start, obs_n, delay_mean,
                                    delay_sd, weights, is_tail ? tail_weight[1] : 0,
                                    is_tail ? tail_mean[1] : 1, hazard_multiplier, rho),
                  eta, wrong_dead[1], rho),
        n_cfr_segments, n_recording_segments,
        n_cfr_segments > 1 ? cfr_p_change[1] : 0.5,
        n_recording_segments > 1 ? recording_p_change[1] : 0.5);
    }
  }
}
generated quantities {
  vector[n_days] Rt;
  vector[n_days] delay_mean_day;
  vector[n_days] delay_sd_day;
  int Rt_n_changes = 0;
  int delay_n_changes = 0;
  vector[is_deaths ? n_days : 0] cfr_day;
  vector[is_deaths ? n_days : 0] recording_mean_day;
  vector[is_deaths ? n_days : 0] recording_sd_day;
  vector[is_deaths ? n_days : 0] deaths_true;   // a draw of the fatal cases of each onset day
  int cfr_n_changes = 0;
  int recording_n_changes = 0;
  vector[7] weekday_effect = exp(log_weekday);
  vector[7] weekday_recording_effect = exp(log_weekday_recording);
  vector[is_report_effects ? n_reports : 0] report_effect;
  if (is_report_effects) {
    report_effect = exp(report_sd[1] * report_raw);
  }
  if (is_Rt_changepoints) {
    array[n_days] int segment = changepoint_segments_rng(
      renewal_ll(cases_true, cases_before, serial_weights, log_Rt,
                 is_negative_binomial, is_negative_binomial ? kappa[1] : 1),
      Rt_p_change[1]);
    Rt = exp(log_Rt[segment]);
    Rt_n_changes = segment[n_days] - 1;
  } else {
    Rt = exp(log_Rt[period]);
  }
  if (is_delay_changepoints) {
    array[n_days] int segment = changepoint_segments_rng(
      reporting_ll(cases_true, obs_day, obs_count, obs_start, obs_n, delay_mean,
                   delay_sd, delay_period, 1, is_tail ? tail_weight[1] : 0,
                   is_tail ? tail_mean[1] : 1,
                   hazard_multipliers(n_days, is_weekday, weekday_report, log_weekday,
                                      is_report_effects, report_index,
                                      is_report_effects ? report_sd[1] * report_raw
                                                        : rep_vector(0, 0)),
                   eta, rho, max_erroneous, 0, shifts, shift_probs),
      delay_p_change[1]);
    delay_mean_day = delay_mean[segment];
    delay_sd_day = delay_sd[segment];
    delay_n_changes = segment[n_days] - 1;
  } else {
    delay_mean_day = delay_mean[delay_period];
    delay_sd_day = delay_sd[delay_period];
  }
  if (is_deaths) {
    vector[n_days] hazard_multiplier = hazard_multipliers(
      n_days, is_weekday, weekday_report, log_weekday, is_report_effects, report_index,
      is_report_effects ? report_sd[1] * report_raw : rep_vector(0, 0));
    vector[n_days] hazard_recording = hazard_multipliers(
      n_days, is_recording_weekday, weekday_report, log_weekday_recording, 0,
      report_index, rep_vector(0, 0));
    matrix[n_days, n_delay_periods] weights = rep_matrix(0, n_days, n_delay_periods);
    array[n_days, 2] int segment;
    if (is_delay_changepoints) {
      weights = changepoint_smoothed(
        reporting_ll(cases_true, obs_day, obs_count, obs_start, obs_n, delay_mean,
                     delay_sd, delay_period, 1, is_tail ? tail_weight[1] : 0,
                     is_tail ? tail_mean[1] : 1, hazard_multiplier, eta, rho,
                     max_erroneous, 0, shifts, shift_probs),
        delay_p_change[1]);
    } else {
      for (t in 1:n_days) {
        weights[t, delay_period[t]] = 1;
      }
    }
    segment = changepoint2_segments_rng(
      deaths_ll(cases_true, obs_day, obs_count, obs_deaths, obs_start, obs_n,
                max_dead, cfr, recording_mean, recording_sd, hazard_recording,
                wrong_entry_units(n_days, obs_day, obs_start, obs_n, delay_mean,
                                  delay_sd, weights, is_tail ? tail_weight[1] : 0,
                                  is_tail ? tail_mean[1] : 1, hazard_multiplier, rho),
                eta, wrong_dead[1], rho),
      n_cfr_segments, n_recording_segments,
      n_cfr_segments > 1 ? cfr_p_change[1] : 0.5,
      n_recording_segments > 1 ? recording_p_change[1] : 0.5);
    cfr_n_changes = segment[n_days, 1] - 1;
    recording_n_changes = segment[n_days, 2] - 1;
    for (t in 1:n_days) {
      cfr_day[t] = cfr[segment[t, 1]];
      recording_mean_day[t] = recording_mean[segment[t, 2]];
      recording_sd_day[t] = recording_sd[segment[t, 2]];
      // those shown as dead, those shown as alive who will die, and those
      // not yet reported who will
      deaths_true[t] = binomial_rng(to_int(round(cases_true[t])), cfr_day[t]);
      if (obs_n[t] > 0) {
        int last = obs_start[t] + obs_n[t] - 1;
        real unrecorded = exp(log_undetected_row(
          t, delay_log_survival(n_days, recording_mean_day[t], recording_sd_day[t], 0, 1),
          hazard_recording)[obs_day[last] - t + 1]);
        int cases_all = to_int(round(cases_true[t]));
        int reported = min(obs_count[last], cases_all);
        int shown_dead = min(obs_deaths[last], reported);
        deaths_true[t] = shown_dead
          + binomial_rng(reported - shown_dead,
                         cfr_day[t] * unrecorded
                         / (1 - cfr_day[t] + cfr_day[t] * unrecorded))
          + binomial_rng(cases_all - reported, cfr_day[t]);
      }
    }
  }
}
