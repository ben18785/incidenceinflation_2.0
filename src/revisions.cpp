#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// binomial probability which is zero outside of support
static inline double dbinom_safe(double x, double size, double prob) {
  if(size < 0 || x < 0 || x > size)
    return 0.0;
  return R::dbinom(x, size, prob, 0);
}

//' Forward recursion for observation process allowing downward revisions
//'
//' Calculates the log-probability of the reported cases for a single onset
//' day, summing over the unobserved count of erroneous records, E, and (if the
//' onset day lies in the footprint of a batch) the unobserved count of records
//' belonging to the batch, X. This is the compiled equivalent of the recursion
//' in \code{observation_process_revisions_logp}, which prepares its inputs.
//'
//' @param cases_true vector of true case counts
//' @param observation_cases vector of reported cases on each reporting day
//' @param p_detect vector giving, for each interval between reporting days,
//' the probability a genuine case not yet recorded is recorded
//' @param p_newly_detected vector giving, for each interval, the probability
//' a genuine case is detected in it
//' @param p_withdrawn vector giving, for each interval, the probability an
//' erroneous record is withdrawn
//' @param eta expected number of spurious records entered per genuine case detected
//' @param mean_entered_extra matrix with a row for each true case count and a
//' column for each interval, giving the expected number of misdated records entered
//' @param log_alpha_initial matrix with a row for each true case count giving
//' the log-probability of each state on the first reporting day. Columns are
//' ordered so that state (e, x) is in column e + (max_erroneous + 1) * x.
//' @param max_erroneous maximum count of erroneous records
//' @param phi intensity of the batch (ignored if max_batch is zero)
//' @param period_detected index of the interval at the end of which the batch
//' is detected; a value exceeding the number of intervals means the batch is
//' not detected within the reporting days
//' @param max_batch maximum count of batch records (zero if no batch)
//'
//' @return a vector of log-probabilities
//' @keywords internal
// [[Rcpp::export]]
NumericVector revisions_logp_cpp(NumericVector cases_true,
                                 NumericVector observation_cases,
                                 NumericVector p_detect,
                                 NumericVector p_newly_detected,
                                 NumericVector p_withdrawn,
                                 double eta,
                                 NumericMatrix mean_entered_extra,
                                 NumericMatrix log_alpha_initial,
                                 int max_erroneous,
                                 double phi,
                                 int period_detected,
                                 int max_batch) {
  const int n_cases = cases_true.size();
  const int n_periods = p_detect.size();
  const int n_e = max_erroneous + 1;
  const int n_x = max_batch + 1;
  const int n_states = n_e * n_x;
  NumericVector logp(n_cases);

  // probabilities of withdrawals, which do not depend on true cases
  std::vector<double> p_m(n_periods * n_e * n_e, 0.0);
  for(int j = 0; j < n_periods; j++)
    for(int e_prev = 0; e_prev < n_e; e_prev++)
      for(int m = 0; m <= e_prev; m++)
        p_m[(j * n_e + e_prev) * n_e + m] = R::dbinom(m, e_prev, p_withdrawn[j], 0);

  std::vector<double> alpha(n_states), alpha_new(n_states);
  std::vector<double> p_entered(n_e), w(n_e * n_e), p_batch(n_x);
  std::vector<double> p_genuine;

  for(int i = 0; i < n_cases; i++) {
    const double cases = cases_true[i];
    double total = 0.0;
    for(int s = 0; s < n_states; s++) {
      alpha[s] = std::exp(log_alpha_initial(i, s));
      total += alpha[s];
    }
    if(!(total > 0.0)) {
      logp[i] = R_NegInf;
      continue;
    }
    double loglik = 0.0;
    bool is_possible = true;

    for(int j = 0; j < n_periods && is_possible; j++) {
      const int period = j + 1;
      const double change = observation_cases[j + 1] - observation_cases[j];
      const double remaining_base = cases - observation_cases[j];
      const double p = p_detect[j];

      // probability erroneous count moves from e_prev to e
      const double mean_entered = eta * cases * p_newly_detected[j] +
        mean_entered_extra(i, j);
      for(int k = 0; k < n_e; k++)
        p_entered[k] = R::dpois(k, mean_entered, 0);
      for(int e_prev = 0; e_prev < n_e; e_prev++) {
        for(int e = 0; e < n_e; e++) {
          double sum = 0.0;
          for(int m = 0; m <= e_prev; m++) {
            const int entered = e - e_prev + m;
            if(entered >= 0 && entered < n_e)
              sum += p_m[(j * n_e + e_prev) * n_e + m] * p_entered[entered];
          }
          w[e_prev * n_e + e] = sum;
        }
      }

      std::fill(alpha_new.begin(), alpha_new.end(), 0.0);
      const bool is_batch = max_batch > 0;

      if(!is_batch || period > period_detected) {
        // no batch records held
        for(int e_prev = 0; e_prev < n_e; e_prev++) {
          if(alpha[e_prev] == 0.0) continue;
          for(int e = 0; e < n_e; e++) {
            const double w_e = w[e_prev * n_e + e];
            if(w_e == 0.0) continue;
            alpha_new[e] += alpha[e_prev] * w_e *
              dbinom_safe(change - (e - e_prev), remaining_base + e_prev, p);
          }
        }
      } else if(period < period_detected) {
        // batch records accumulate
        const double mean_batch = phi * cases * p_newly_detected[j];
        for(int k = 0; k < n_x; k++)
          p_batch[k] = R::dpois(k, mean_batch, 0);
        // genuine cases newly recorded, tabulated by number yet to be recorded
        // (which depends on e_prev + x_prev) and by count
        const int n_r = n_e + n_x - 1;
        const int d_max = (int)change + max_erroneous;
        const int n_d = d_max >= 0 ? d_max + 1 : 0;
        p_genuine.assign(n_r * (n_d > 0 ? n_d : 1), -1.0);
        for(int e_prev = 0; e_prev < n_e; e_prev++) {
          for(int x_prev = 0; x_prev < n_x; x_prev++) {
            const double a = alpha[e_prev + n_e * x_prev];
            if(a == 0.0) continue;
            const int r = e_prev + x_prev;
            for(int e = 0; e < n_e; e++) {
              const double w_e = w[e_prev * n_e + e];
              if(w_e == 0.0) continue;
              const int d_base = (int)change - (e - e_prev);
              for(int x = x_prev; x < n_x; x++) {
                const int d = d_base - (x - x_prev);
                if(d < 0) break;
                double p_d = p_genuine[r * n_d + d];
                if(p_d < 0.0) {
                  p_d = dbinom_safe(d, remaining_base + r, p);
                  p_genuine[r * n_d + d] = p_d;
                }
                alpha_new[e + n_e * x] += a * w_e * p_batch[x - x_prev] * p_d;
              }
            }
          }
        }
      } else {
        // batch detected: all of its records are withdrawn
        for(int e_prev = 0; e_prev < n_e; e_prev++) {
          for(int x_prev = 0; x_prev < n_x; x_prev++) {
            const double a = alpha[e_prev + n_e * x_prev];
            if(a == 0.0) continue;
            for(int e = 0; e < n_e; e++) {
              const double w_e = w[e_prev * n_e + e];
              if(w_e == 0.0) continue;
              alpha_new[e] += a * w_e *
                dbinom_safe(change - (e - e_prev) + x_prev,
                            remaining_base + e_prev + x_prev, p);
            }
          }
        }
      }

      total = 0.0;
      for(int s = 0; s < n_states; s++)
        total += alpha_new[s];
      if(!(total > 0.0)) {
        is_possible = false;
        break;
      }
      loglik += std::log(total);
      for(int s = 0; s < n_states; s++)
        alpha[s] = alpha_new[s] / total;
    }
    logp[i] = is_possible ? loglik : R_NegInf;
  }
  return logp;
}
