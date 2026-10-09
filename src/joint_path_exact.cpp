#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// As JointPath (joint_path.cpp), but the fatal cases dated by onset and those
// dated by death are followed separately: the unobserved state at each report
// is the number of each which has been reported and the numbers of wrong
// records present which are shown as alive and as dead. A case dated by death
// is shown as dead in the first report it appears in.
class JointPathExact {
public:
  int n_reports, fatal, dated, survivors, max_erroneous, n_wrong, r_max, d_max;
  std::vector<int> cases, deaths, wrong_alive, wrong_dead;
  std::vector<double> p, p_fatal, p_dated, g, q, lambda, w, log_factorial;
  double delta;
  std::vector<double> kernel;

  JointPathExact(NumericVector cases_, NumericVector deaths_, int fatal_,
                 int dated_, int survivors_, NumericVector p_,
                 NumericVector p_fatal_, NumericVector p_dated_,
                 NumericVector g_, NumericVector q_, NumericVector lambda_,
                 NumericVector w_, double delta_, int max_erroneous_) :
    n_reports(cases_.size()), fatal(fatal_), dated(dated_), survivors(survivors_),
    max_erroneous(max_erroneous_), delta(delta_) {
    // counts before the first report are zero
    cases.push_back(0);
    deaths.push_back(0);
    int cases_max = 0;
    int deaths_max = 0;
    for(int j = 0; j < n_reports; j++) {
      cases.push_back((int) std::lround(cases_[j]));
      deaths.push_back((int) std::lround(deaths_[j]));
      if(cases[j + 1] > cases_max)
        cases_max = cases[j + 1];
      if(deaths[j + 1] > deaths_max)
        deaths_max = deaths[j + 1];
    }
    p_dated = as<std::vector<double> >(p_dated_);
    p = as<std::vector<double> >(p_);
    p_fatal = as<std::vector<double> >(p_fatal_);
    g = as<std::vector<double> >(g_);
    q = as<std::vector<double> >(q_);
    lambda = as<std::vector<double> >(lambda_);
    w = as<std::vector<double> >(w_);
    for(int total = 0; total <= max_erroneous; total++)
      for(int dead = 0; dead <= total; dead++) {
        wrong_alive.push_back(total - dead);
        wrong_dead.push_back(dead);
      }
    n_wrong = wrong_alive.size();
    r_max = std::min(cases_max, fatal);
    if(r_max < 0)
      r_max = 0;
    // a case dated by death is shown as dead from the first
    d_max = std::min(deaths_max, dated);
    if(d_max < 0)
      d_max = 0;
    int n_max = std::max(std::max(std::max(fatal, dated), survivors), cases_max) +
      max_erroneous + 2;
    log_factorial.resize(n_max + 1);
    log_factorial[0] = 0.0;
    for(int i = 1; i <= n_max; i++)
      log_factorial[i] = log_factorial[i - 1] + std::log((double) i);
    kernel.resize(n_wrong * n_wrong);
  }

  // binomial probability which is zero outside of support
  inline double dbin(int x, int size, double prob) const {
    if(size < 0 || x < 0 || x > size)
      return 0.0;
    if(prob <= 0.0)
      return x == 0 ? 1.0 : 0.0;
    if(prob >= 1.0)
      return x == size ? 1.0 : 0.0;
    return std::exp(log_factorial[size] - log_factorial[x] -
                    log_factorial[size - x] + x * std::log(prob) +
                    (size - x) * std::log1p(-prob));
  }

  inline double dpois(int x, double mean) const {
    if(x < 0)
      return 0.0;
    if(mean <= 0.0)
      return x == 0 ? 1.0 : 0.0;
    return std::exp(x * std::log(mean) - mean - log_factorial[x]);
  }

  // probability of the wrong records removed and entered between two reports
  inline double wrong_term(int j, int from, int removed_alive, int removed_dead,
                           int to) const {
    int entered_alive = wrong_alive[to] - wrong_alive[from] + removed_alive;
    int entered_dead = wrong_dead[to] - wrong_dead[from] + removed_dead;
    if(entered_alive < 0 || entered_dead < 0)
      return 0.0;
    return dbin(removed_alive, wrong_alive[from], w[j - 1]) *
      dbin(removed_dead, wrong_dead[from], w[j - 1]) *
      dpois(entered_alive, lambda[j - 1] * (1.0 - delta)) *
      dpois(entered_dead, lambda[j - 1] * delta);
  }

  // for report j (counted from one), the probability of each change in the
  // wrong records present, summed over those removed
  void set_kernel(int j) {
    for(int from = 0; from < n_wrong; from++)
      for(int to = 0; to < n_wrong; to++) {
        double total = 0.0;
        for(int removed_alive = 0; removed_alive <= wrong_alive[from]; removed_alive++)
          for(int removed_dead = 0; removed_dead <= wrong_dead[from]; removed_dead++)
            total += wrong_term(j, from, removed_alive, removed_dead, to);
        kernel[from * n_wrong + to] = total;
      }
  }

  // probability of the deaths newly shown: those of fatal cases newly
  // reported and of fatal cases reported earlier and shown as alive
  inline double deaths_new(int j, int deaths_new, int fatal_new, int alive) const {
    if(deaths_new < 0)
      return 0.0;
    double total = 0.0;
    int x_min = std::max(0, deaths_new - alive);
    int x_max = std::min(fatal_new, deaths_new);
    for(int x = x_min; x <= x_max; x++)
      total += dbin(x, fatal_new, g[j - 1]) * dbin(deaths_new - x, alive, q[j - 1]);
    return total;
  }

  // probability of the genuine cases and deaths newly reported in report j
  // given the states at the report before and at this one (the wrong records
  // aside). r counts the fatal cases dated by onset which have been reported
  // and d those dated by death. Zero if the states are impossible given the
  // reports.
  inline double genuine(int j, int from, int r_from, int d_from, int to,
                        int r_to, int d_to) const {
    int fatal_new = r_to - r_from;
    int dated_new = d_to - d_from;
    int wrong_from = wrong_alive[from] + wrong_dead[from];
    int wrong_to = wrong_alive[to] + wrong_dead[to];
    int survivors_from = cases[j - 1] - r_from - d_from - wrong_from;
    int survivors_to = cases[j] - r_to - d_to - wrong_to;
    int survivors_new = survivors_to - survivors_from;
    // fatal cases dated by onset which are shown as dead
    int deaths_from = deaths[j - 1] - wrong_dead[from] - d_from;
    int deaths_to = deaths[j] - wrong_dead[to] - d_to;
    if(fatal_new < 0 || dated_new < 0 || survivors_new < 0 || survivors_from < 0 ||
       r_to > fatal || d_to > dated || survivors_to > survivors ||
       deaths_from < 0 || deaths_to < deaths_from || deaths_to > r_to ||
       deaths_from > r_from)
      return 0.0;
    return dbin(fatal_new, fatal - r_from, p_fatal[j - 1]) *
      dbin(dated_new, dated - d_from, p_dated[j - 1]) *
      dbin(survivors_new, survivors - survivors_from, p[j - 1]) *
      deaths_new(j, deaths_to - deaths_from, fatal_new, r_from - deaths_from);
  }

  inline int n_states() const {
    return n_wrong * (d_max + 1) * (r_max + 1);
  }

  inline int index(int j, int state, int d, int r) const {
    return ((j * n_wrong + state) * (d_max + 1) + d) * (r_max + 1) + r;
  }

  // filters forwards; returns the log-probability of the reports
  double forward(std::vector<double>& alpha) {
    const double negligible = 1e-12;
    alpha.assign((n_reports + 1) * n_states(), 0.0);
    alpha[index(0, 0, 0, 0)] = 1.0;
    double logp = 0.0;
    for(int j = 1; j <= n_reports; j++) {
      set_kernel(j);
      int cases_new = cases[j] - cases[j - 1];
      double largest = 0.0;
      for(int i = 0; i < n_states(); i++)
        largest = std::max(largest, alpha[(j - 1) * n_states() + i]);
      for(int from = 0; from < n_wrong; from++) {
        int wrong_from = wrong_alive[from] + wrong_dead[from];
        for(int d = 0; d <= d_max; d++)
          for(int r = 0; r <= r_max; r++) {
            double weight = alpha[index(j - 1, from, d, r)];
            if(weight <= negligible * largest)
              continue;
            for(int to = 0; to < n_wrong; to++) {
              double weight_wrong = kernel[from * n_wrong + to];
              if(weight_wrong <= 0.0)
                continue;
              int wrong_to = wrong_alive[to] + wrong_dead[to];
              int genuine_new = cases_new - (wrong_to - wrong_from);
              for(int a_d = 0; a_d <= std::min(genuine_new, d_max - d); a_d++)
                for(int a = 0; a <= std::min(genuine_new - a_d, r_max - r); a++) {
                  double value = genuine(j, from, r, d, to, r + a, d + a_d);
                  if(value > 0.0)
                    alpha[index(j, to, d + a_d, r + a)] += weight * weight_wrong * value;
                }
            }
          }
      }
      double total = 0.0;
      for(int i = 0; i < n_states(); i++)
        total += alpha[j * n_states() + i];
      if(!(total > 0.0))
        return R_NegInf;
      for(int i = 0; i < n_states(); i++)
        alpha[j * n_states() + i] /= total;
      logp += std::log(total);
    }
    return logp;
  }

  // draws an index with probability proportional to the weights
  static int draw(const std::vector<double>& weights) {
    double total = 0.0;
    for(size_t i = 0; i < weights.size(); i++)
      total += weights[i];
    double u = R::unif_rand() * total;
    double cumulative = 0.0;
    int last = 0;
    for(size_t i = 0; i < weights.size(); i++) {
      if(weights[i] <= 0.0)
        continue;
      cumulative += weights[i];
      last = i;
      if(u < cumulative)
        return i;
    }
    return last;
  }
};

//' Filters forwards over the reports of the cases recorded for one day,
//' following fatal cases dated by onset and by death separately
//'
//' @inheritParams joint_path_forward_cpp
//' @param fatal the number of fatal cases recorded for the day which are
//' dated by onset
//' @param dated the number which are dated by death
//' @param p_fatal for each report, the probability a fatal case dated by
//' onset which was not reported by the report before is reported by it
//' @param p_dated the same for a fatal case dated by death
//' @param g for each report, the probability a fatal case dated by onset
//' which is newly reported in it is shown as dead in it
//'
//' @return as \code{joint_path_forward_cpp}
// [[Rcpp::export]]
List joint_path_exact_forward_cpp(NumericVector cases, NumericVector deaths,
                                  int fatal, int dated, int survivors,
                                  NumericVector p, NumericVector p_fatal,
                                  NumericVector p_dated, NumericVector g,
                                  NumericVector q, NumericVector lambda,
                                  NumericVector w, double delta,
                                  int max_erroneous) {
  JointPathExact path(cases, deaths, fatal, dated, survivors, p, p_fatal,
                      p_dated, g, q, lambda, w, delta, max_erroneous);
  std::vector<double> alpha;
  double logp = path.forward(alpha);
  return List::create(Named("logp") = logp, Named("alpha") = wrap(alpha));
}

//' Samples backwards over the reports of the cases recorded for one day,
//' following fatal cases dated by onset and by death separately
//'
//' @inheritParams joint_path_exact_forward_cpp
//' @param alpha as returned by \code{joint_path_exact_forward_cpp} for the
//' same arguments
//'
//' @return a matrix with a row for each report and columns: the fatal cases
//' dated by onset reported by it, those dated by death reported by it, the
//' wrong records present which are shown as alive and as dead, and the wrong
//' records of each kind removed since the report before
// [[Rcpp::export]]
NumericMatrix joint_path_exact_backward_cpp(NumericVector alpha,
                                            NumericVector cases,
                                            NumericVector deaths, int fatal,
                                            int dated, int survivors,
                                            NumericVector p,
                                            NumericVector p_fatal,
                                            NumericVector p_dated,
                                            NumericVector g, NumericVector q,
                                            NumericVector lambda,
                                            NumericVector w, double delta,
                                            int max_erroneous) {
  JointPathExact path(cases, deaths, fatal, dated, survivors, p, p_fatal,
                      p_dated, g, q, lambda, w, delta, max_erroneous);
  int n_reports = path.n_reports;
  int n_wrong = path.n_wrong;
  int n_r = path.r_max + 1;
  int n_d = path.d_max + 1;
  int n_states = path.n_states();
  NumericMatrix out(n_reports, 6);

  std::vector<double> weights(n_states);
  for(int i = 0; i < n_states; i++)
    weights[i] = alpha[n_reports * n_states + i];
  int chosen = JointPathExact::draw(weights);
  int to = chosen / (n_d * n_r);
  int d_to = (chosen / n_r) % n_d;
  int r_to = chosen % n_r;

  for(int j = n_reports; j >= 1; j--) {
    out(j - 1, 0) = r_to;
    out(j - 1, 1) = d_to;
    out(j - 1, 2) = path.wrong_alive[to];
    out(j - 1, 3) = path.wrong_dead[to];
    path.set_kernel(j);
    for(int from = 0; from < n_wrong; from++)
      for(int d = 0; d < n_d; d++)
        for(int r = 0; r < n_r; r++) {
          double value = 0.0;
          double weight = alpha[path.index(j - 1, from, d, r)];
          if(weight > 0.0 && r <= r_to && d <= d_to)
            value = weight * path.kernel[from * n_wrong + to] *
              path.genuine(j, from, r, d, to, r_to, d_to);
          weights[(from * n_d + d) * n_r + r] = value;
        }
    chosen = JointPathExact::draw(weights);
    int from = chosen / (n_d * n_r);
    int d_from = (chosen / n_r) % n_d;
    int r_from = chosen % n_r;
    std::vector<double> weights_removed;
    std::vector<int> removed_alive, removed_dead;
    for(int alive = 0; alive <= path.wrong_alive[from]; alive++)
      for(int dead = 0; dead <= path.wrong_dead[from]; dead++) {
        weights_removed.push_back(path.wrong_term(j, from, alive, dead, to));
        removed_alive.push_back(alive);
        removed_dead.push_back(dead);
      }
    int removed = JointPathExact::draw(weights_removed);
    out(j - 1, 4) = removed_alive[removed];
    out(j - 1, 5) = removed_dead[removed];
    to = from;
    d_to = d_from;
    r_to = r_from;
  }
  return out;
}
