#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// The reports of the cases recorded for one day in the joint model of cases
// and deaths. The unobserved state at each report is the number of fatal
// cases reported so far and the numbers of wrong records present which are
// shown as alive and as dead.
class JointPath {
public:
  int n_reports, fatal, survivors, max_erroneous, n_wrong, r_max;
  std::vector<int> cases, deaths, wrong_alive, wrong_dead;
  std::vector<double> p, p_fatal, g, q, lambda, w, log_factorial;
  double delta;
  std::vector<double> kernel;

  JointPath(NumericVector cases_, NumericVector deaths_, int fatal_,
            int survivors_, NumericVector p_, NumericVector p_fatal_,
            NumericVector g_, NumericVector q_, NumericVector lambda_,
            NumericVector w_, double delta_, int max_erroneous_) :
    n_reports(cases_.size()), fatal(fatal_), survivors(survivors_),
    max_erroneous(max_erroneous_), delta(delta_) {
    // counts before the first report are zero
    cases.push_back(0);
    deaths.push_back(0);
    int cases_max = 0;
    for(int j = 0; j < n_reports; j++) {
      cases.push_back((int) std::lround(cases_[j]));
      deaths.push_back((int) std::lround(deaths_[j]));
      if(cases[j + 1] > cases_max)
        cases_max = cases[j + 1];
    }
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
    int n_max = std::max(std::max(fatal, survivors), cases_max) + max_erroneous + 2;
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
  // aside). Zero if the states are impossible given the reports.
  inline double genuine(int j, int from, int r_from, int to, int r_to) const {
    int fatal_new = r_to - r_from;
    int wrong_from = wrong_alive[from] + wrong_dead[from];
    int wrong_to = wrong_alive[to] + wrong_dead[to];
    int survivors_from = cases[j - 1] - r_from - wrong_from;
    int survivors_to = cases[j] - r_to - wrong_to;
    int survivors_new = survivors_to - survivors_from;
    int deaths_from = deaths[j - 1] - wrong_dead[from];
    int deaths_to = deaths[j] - wrong_dead[to];
    if(fatal_new < 0 || survivors_new < 0 || survivors_from < 0 ||
       r_to > fatal || survivors_to > survivors ||
       deaths_from < 0 || deaths_to < deaths_from || deaths_to > r_to ||
       deaths_from > r_from)
      return 0.0;
    return dbin(fatal_new, fatal - r_from, p_fatal[j - 1]) *
      dbin(survivors_new, survivors - survivors_from, p[j - 1]) *
      deaths_new(j, deaths_to - deaths_from, fatal_new, r_from - deaths_from);
  }

  inline int index(int j, int state, int r) const {
    return (j * n_wrong + state) * (r_max + 1) + r;
  }

  // filters forwards; returns the log-probability of the reports
  double forward(std::vector<double>& alpha) {
    const double negligible = 1e-12;
    alpha.assign((n_reports + 1) * n_wrong * (r_max + 1), 0.0);
    alpha[index(0, 0, 0)] = 1.0;
    double logp = 0.0;
    for(int j = 1; j <= n_reports; j++) {
      set_kernel(j);
      int cases_new = cases[j] - cases[j - 1];
      double largest = 0.0;
      for(int from = 0; from < n_wrong; from++)
        for(int r = 0; r <= r_max; r++)
          largest = std::max(largest, alpha[index(j - 1, from, r)]);
      for(int from = 0; from < n_wrong; from++) {
        int wrong_from = wrong_alive[from] + wrong_dead[from];
        for(int r = 0; r <= r_max; r++) {
          double weight = alpha[index(j - 1, from, r)];
          if(weight <= negligible * largest)
            continue;
          for(int to = 0; to < n_wrong; to++) {
            double weight_wrong = kernel[from * n_wrong + to];
            if(weight_wrong <= 0.0)
              continue;
            int wrong_to = wrong_alive[to] + wrong_dead[to];
            // genuine cases newly reported
            int genuine_new = cases_new - (wrong_to - wrong_from);
            int a_max = std::min(genuine_new, r_max - r);
            for(int a = 0; a <= a_max; a++) {
              double value = genuine(j, from, r, to, r + a);
              if(value > 0.0)
                alpha[index(j, to, r + a)] += weight * weight_wrong * value;
            }
          }
        }
      }
      double total = 0.0;
      for(int state = 0; state < n_wrong; state++)
        for(int r = 0; r <= r_max; r++)
          total += alpha[index(j, state, r)];
      if(!(total > 0.0))
        return R_NegInf;
      for(int state = 0; state < n_wrong; state++)
        for(int r = 0; r <= r_max; r++)
          alpha[index(j, state, r)] /= total;
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
//' allowing for wrong records
//'
//' The unobserved state at each report is the number of fatal cases reported
//' so far and the numbers of wrong records present which are shown as alive
//' and as dead, of which there are at most max_erroneous in all. States
//' whose probability is negligible given the reports so far are dropped.
//'
//' @param cases the cases in each report
//' @param deaths the deaths in each report
//' @param fatal the number of fatal cases recorded for the day
//' @param survivors the number of cases recorded for the day which are not fatal
//' @param p for each report, the probability a case which is not fatal and
//' was not reported by the report before is reported by it
//' @param p_fatal the same for a fatal case
//' @param g for each report, the probability a fatal case newly reported in
//' it is shown as dead in it
//' @param q for each report, the probability a death not recorded by the
//' report before is recorded by it
//' @param lambda for each report, the expected number of wrong records entered
//' since the report before (for the first report, of those entered which are
//' still present)
//' @param w for each report, the probability a wrong record present at the
//' report before has been removed
//' @param delta the probability a wrong record is shown as dead
//' @param max_erroneous the largest number of wrong records present at once
//'
//' @return a named list of 'logp' (the log-probability of the reports) and
//' 'alpha' (the probabilities of the states at each report given the reports
//' up to it)
// [[Rcpp::export]]
List joint_path_forward_cpp(NumericVector cases, NumericVector deaths,
                            int fatal, int survivors, NumericVector p,
                            NumericVector p_fatal, NumericVector g,
                            NumericVector q, NumericVector lambda,
                            NumericVector w, double delta, int max_erroneous) {
  JointPath path(cases, deaths, fatal, survivors, p, p_fatal, g, q, lambda, w,
                 delta, max_erroneous);
  std::vector<double> alpha;
  double logp = path.forward(alpha);
  return List::create(Named("logp") = logp, Named("alpha") = wrap(alpha));
}

//' Samples backwards over the reports of the cases recorded for one day
//'
//' @inheritParams joint_path_forward_cpp
//' @param alpha as returned by \code{joint_path_forward_cpp} for the same
//' arguments
//'
//' @return a matrix with a row for each report and columns: the fatal cases
//' reported by it, the wrong records present which are shown as alive and as
//' dead, and the wrong records of each kind removed since the report before
// [[Rcpp::export]]
NumericMatrix joint_path_backward_cpp(NumericVector alpha, NumericVector cases,
                                      NumericVector deaths, int fatal,
                                      int survivors, NumericVector p,
                                      NumericVector p_fatal, NumericVector g,
                                      NumericVector q, NumericVector lambda,
                                      NumericVector w, double delta,
                                      int max_erroneous) {
  JointPath path(cases, deaths, fatal, survivors, p, p_fatal, g, q, lambda, w,
                 delta, max_erroneous);
  int n_reports = path.n_reports;
  int n_wrong = path.n_wrong;
  int r_max = path.r_max;
  NumericMatrix out(n_reports, 5);

  std::vector<double> weights(n_wrong * (r_max + 1));
  for(int state = 0; state < n_wrong; state++)
    for(int r = 0; r <= r_max; r++)
      weights[state * (r_max + 1) + r] = alpha[path.index(n_reports, state, r)];
  int chosen = JointPath::draw(weights);
  int to = chosen / (r_max + 1);
  int r_to = chosen % (r_max + 1);

  for(int j = n_reports; j >= 1; j--) {
    out(j - 1, 0) = r_to;
    out(j - 1, 1) = path.wrong_alive[to];
    out(j - 1, 2) = path.wrong_dead[to];
    path.set_kernel(j);
    for(int from = 0; from < n_wrong; from++)
      for(int r = 0; r <= r_max; r++) {
        double value = 0.0;
        double weight = alpha[path.index(j - 1, from, r)];
        if(weight > 0.0 && r <= r_to)
          value = weight * path.kernel[from * n_wrong + to] *
            path.genuine(j, from, r, to, r_to);
        weights[from * (r_max + 1) + r] = value;
      }
    chosen = JointPath::draw(weights);
    int from = chosen / (r_max + 1);
    int r_from = chosen % (r_max + 1);
    // the wrong records removed, given those present before and after
    std::vector<double> weights_removed;
    std::vector<int> removed_alive, removed_dead;
    for(int alive = 0; alive <= path.wrong_alive[from]; alive++)
      for(int dead = 0; dead <= path.wrong_dead[from]; dead++) {
        weights_removed.push_back(path.wrong_term(j, from, alive, dead, to));
        removed_alive.push_back(alive);
        removed_dead.push_back(dead);
      }
    int removed = JointPath::draw(weights_removed);
    out(j - 1, 3) = removed_alive[removed];
    out(j - 1, 4) = removed_dead[removed];
    to = from;
    r_to = r_from;
  }
  return out;
}
