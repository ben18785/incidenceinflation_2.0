// A relaxed version of the model in which nothing is discrete, so that it can
// be fitted by Hamiltonian Monte Carlo or by variational methods:
// - true case counts of recent onset days are continuous (those of older onset
//   days are fixed at their latest reported counts);
// - Rt follows a random walk between periods of onset days, in place of
//   changes at unknown times;
// - there is a single reporting delay distribution.
// The model of reporting is otherwise that of the package: a gamma delay, an
// effect of the day of the week of each report, and erroneous records which
// are entered and later removed (summed over by a forward recursion).
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

  // log-probability of the reported counts of one onset day given its true
  // cases. log_undetected[k] is the log-probability that a case remains
  // unreported at the end of day (day_onset + k - 1).
  real observation_logp(real cases_true, array[] int days, array[] int counts,
                      int day_onset, vector log_undetected, real eta, real rho,
                      int max_erroneous) {
    int n = size(days);
    vector[max_erroneous + 1] log_alpha = rep_vector(negative_infinity(),
                                                     max_erroneous + 1);
    if (n <= 1) {
      return 0;
    }
    // erroneous records present at the first report, given its count
    if (days[1] == day_onset || eta == 0) {
      log_alpha[1] = 0;
    } else {
      real mean_erroneous = 0;
      real p_detect = -expm1(log_undetected[days[1] - day_onset + 1]
                             - log_undetected[1]);
      for (k in (day_onset + 1):days[1]) {
        mean_erroneous += (exp(log_undetected[k - day_onset])
                           - exp(log_undetected[k - day_onset + 1]))
                          * pow(1 - rho, days[1] - k);
      }
      mean_erroneous *= eta * cases_true;
      for (e in 0:min(max_erroneous, counts[1])) {
        if (counts[1] - e <= cases_true) {
          log_alpha[e + 1] = poisson_lpmf(e | mean_erroneous)
                             + binomial_real_logp(counts[1] - e, cases_true, p_detect);
        }
      }
      log_alpha -= log_sum_exp(log_alpha);
    }

    for (j in 2:n) {
      real log_before = log_undetected[days[j - 1] - day_onset + 1];
      real log_after = log_undetected[days[j] - day_onset + 1];
      real p_reported = -expm1(log_after - log_before);
      real mean_entered = eta * cases_true * (exp(log_before) - exp(log_after));
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
    vector[n_days] out;
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
    out = log(survival[1:n_days]);
    return out;
  }

  // log-probability that a case of each onset day (rows) remains unreported
  // at the end of each later day (column k is day onset + k - 1), with the
  // daily hazard of reporting multiplied by that day's hazard_multiplier
  matrix log_undetected_matrix(int n_days, real delay_mean, real delay_sd,
                               vector hazard_multiplier) {
    real shape = square(delay_mean / delay_sd);
    real rate = delay_mean / square(delay_sd);
    vector[n_days] log_base;
    matrix[n_days, n_days] out = rep_matrix(0, n_days, n_days);
    log_base = gamma_log_survival(n_days, shape, rate);
    for (t in 1:n_days) {
      out[t, 1] = log_base[1];
      for (k in (t + 1):n_days) {
        out[t, k - t + 1] = out[t, k - t]
                            - hazard_multiplier[k] * (log_base[k - t] - log_base[k - t + 1]);
      }
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

  int<lower=1> n_periods;                   // periods within which Rt is constant
  array[n_days] int<lower=1, upper=n_periods> period;
  array[n_days] int<lower=0, upper=7> weekday_report; // day of week of the first report
                                            // on or after each day (0 if none)
  int<lower=0, upper=1> is_weekday;         // whether to include the effect

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
  real<lower=0> log_Rt_first_sd;            // sd of log Rt in the first period
  real<lower=0> Rt_step_sd;                 // scale of the prior on the sd of steps in log Rt
}
parameters {
  vector<lower=0>[n_sampled] cases_extra;   // true cases above their lower bound
  real<lower=0> delay_mean;
  real<lower=0> delay_sd;
  sum_to_zero_vector[7] log_weekday;
  real<lower=0> eta;
  real<lower=0, upper=1> rho;
  real<lower=0> kappa;
  real log_Rt_first;
  vector[n_periods - 1] Rt_step_raw;
  real<lower=0> Rt_step;
}
transformed parameters {
  vector[n_days] cases_true = cases_fixed;
  vector[n_periods] log_Rt;
  cases_true[sampled_day] = cases_lower[sampled_day] + cases_extra;
  log_Rt[1] = log_Rt_first;
  for (p in 2:n_periods) {
    log_Rt[p] = log_Rt[p - 1] + Rt_step * Rt_step_raw[p - 1];
  }
}
model {
  vector[n_days] hazard_multiplier = rep_vector(1, n_days);
  matrix[n_days, n_days] log_undetected;
  vector[serial_max + n_days] cases_all = append_row(cases_before, cases_true);

  if (is_weekday) {
    for (k in 1:n_days) {
      if (weekday_report[k] > 0) {
        hazard_multiplier[k] = exp(log_weekday[weekday_report[k]]);
      }
    }
  }
  log_undetected = log_undetected_matrix(n_days, delay_mean, delay_sd, hazard_multiplier);

  // priors (gamma distributions given by mean and sd)
  delay_mean ~ gamma(square(delay_mean_mu / delay_mean_sigma),
                     delay_mean_mu / square(delay_mean_sigma));
  delay_sd ~ gamma(square(delay_sd_mu / delay_sd_sigma),
                   delay_sd_mu / square(delay_sd_sigma));
  eta ~ gamma(square(eta_mu / eta_sigma), eta_mu / square(eta_sigma));
  rho ~ beta(rho_alpha, rho_beta);
  kappa ~ gamma(square(kappa_mu / kappa_sigma), kappa_mu / square(kappa_sigma));
  log_weekday ~ normal(0, weekday_sd);
  log_Rt_first ~ normal(0, log_Rt_first_sd);
  Rt_step_raw ~ std_normal();
  Rt_step ~ normal(0, Rt_step_sd);

  // renewal model
  for (t in 1:n_days) {
    real infectiousness = dot_product(
      serial_weights, reverse(cases_all[t:(t + serial_max - 1)]));
    if (infectiousness > 0) {
      target += neg_binomial_real_logp(cases_true[t],
                                     exp(log_Rt[period[t]]) * infectiousness, kappa);
    }
  }

  // reporting
  for (t in 1:n_days) {
    if (obs_n[t] > 1) {
      target += observation_logp(
        cases_true[t],
        obs_day[obs_start[t]:(obs_start[t] + obs_n[t] - 1)],
        obs_count[obs_start[t]:(obs_start[t] + obs_n[t] - 1)],
        t, log_undetected[t]', eta, rho, max_erroneous);
    }
  }
}
generated quantities {
  vector[n_days] Rt = exp(log_Rt[period]);
  vector[7] weekday_effect = exp(log_weekday);
}
