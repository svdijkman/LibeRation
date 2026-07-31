// Native noncompartmental pharmacokinetic calculations.
//
// This implementation is deliberately independent of the optional R
// reference backend.  Keeping the numerical core here avoids copying GPL
// implementation code into LibeRation and gives simulations a small,
// allocation-light path for repeated profiles.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr double kNa = std::numeric_limits<double>::quiet_NaN();

struct IntervalArea {
  double auc;
  double aumc;
};

IntervalArea interval_area(double t0, double c0, double t1, double c1,
                           bool log_down) {
  const double dt = t1 - t0;
  if (!(dt > 0.0)) return {0.0, 0.0};
  if (log_down && c0 > 0.0 && c1 > 0.0 && c1 < c0) {
    const double lambda = std::log(c0 / c1) / dt;
    if (std::abs(lambda) > 1e-12) {
      const double decay = std::exp(-lambda * dt);
      const double auc = c0 * (1.0 - decay) / lambda;
      const double first = c0 * (1.0 - decay * (1.0 + lambda * dt)) /
                           (lambda * lambda);
      return {auc, t0 * auc + first};
    }
  }
  const double slope = (c1 - c0) / dt;
  const double auc = c0 * dt + 0.5 * slope * dt * dt;
  const double aumc = t0 * auc + 0.5 * c0 * dt * dt +
                      slope * dt * dt * dt / 3.0;
  return {auc, aumc};
}

double interpolate(double t0, double c0, double t1, double c1, double at,
                   bool log_down) {
  if (at <= t0) return c0;
  if (at >= t1) return c1;
  const double fraction = (at - t0) / (t1 - t0);
  if (log_down && c0 > 0.0 && c1 > 0.0 && c1 < c0) {
    return c0 * std::exp(std::log(c1 / c0) * fraction);
  }
  return c0 + fraction * (c1 - c0);
}

double partial_auc(const std::vector<double>& time,
                   const std::vector<double>& concentration, double start,
                   double end, bool log_down) {
  if (!(end > start) || time.empty() || start < time.front() ||
      end > time.back()) return kNa;
  std::vector<double> x;
  std::vector<double> y;
  x.push_back(start);
  for (std::size_t i = 0; i + 1 < time.size(); ++i) {
    if (time[i] < start && start < time[i + 1]) {
      y.push_back(interpolate(time[i], concentration[i], time[i + 1],
                              concentration[i + 1], start, log_down));
      break;
    }
    if (time[i] == start) {
      y.push_back(concentration[i]);
      break;
    }
  }
  if (y.empty() && start == time.back()) y.push_back(concentration.back());
  for (std::size_t i = 0; i < time.size(); ++i) {
    if (time[i] > start && time[i] < end) {
      x.push_back(time[i]);
      y.push_back(concentration[i]);
    }
  }
  x.push_back(end);
  bool found_end = false;
  for (std::size_t i = 0; i + 1 < time.size(); ++i) {
    if (time[i] < end && end < time[i + 1]) {
      y.push_back(interpolate(time[i], concentration[i], time[i + 1],
                              concentration[i + 1], end, log_down));
      found_end = true;
      break;
    }
    if (time[i] == end) {
      y.push_back(concentration[i]);
      found_end = true;
      break;
    }
  }
  if (!found_end && end == time.back()) y.push_back(concentration.back());
  if (x.size() != y.size()) return kNa;
  double area = 0.0;
  for (std::size_t i = 1; i < x.size(); ++i) {
    area += interval_area(x[i - 1], y[i - 1], x[i], y[i], log_down).auc;
  }
  return area;
}

struct Regression {
  bool valid = false;
  double slope = kNa;
  double intercept = kNa;
  double r2 = kNa;
  double adjusted_r2 = kNa;
};

Regression log_regression(const std::vector<double>& time,
                          const std::vector<double>& concentration,
                          const std::vector<int>& indices) {
  const int n = static_cast<int>(indices.size());
  if (n < 3) return {};
  double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
  for (int index : indices) {
    if (index < 0 || index >= static_cast<int>(time.size()) ||
        !(concentration[index] > 0.0)) return {};
    const double x = time[index];
    const double y = std::log(concentration[index]);
    sx += x; sy += y; sxx += x * x; sxy += x * y;
  }
  const double denominator = n * sxx - sx * sx;
  if (std::abs(denominator) <= 1e-14) return {};
  const double slope = (n * sxy - sx * sy) / denominator;
  const double intercept = (sy - slope * sx) / n;
  double sse = 0.0, sst = 0.0;
  const double mean = sy / n;
  for (int index : indices) {
    const double y = std::log(concentration[index]);
    const double residual = y - (intercept + slope * time[index]);
    sse += residual * residual;
    const double centered = y - mean;
    sst += centered * centered;
  }
  const double r2 = sst > 0.0 ? 1.0 - sse / sst : 1.0;
  const double adjusted = n > 2 ? 1.0 - (1.0 - r2) * (n - 1.0) / (n - 2.0) : r2;
  return {slope < 0.0, slope, intercept, r2, adjusted};
}

std::vector<int> automatic_terminal(const std::vector<double>& time,
                                    const std::vector<double>& concentration,
                                    int maximum_index) {
  std::vector<int> positive;
  for (int i = std::max(0, maximum_index); i < static_cast<int>(time.size()); ++i) {
    if (concentration[i] > 0.0) positive.push_back(i);
  }
  if (positive.size() < 3) return {};
  std::vector<int> best;
  double best_adjusted = -std::numeric_limits<double>::infinity();
  for (std::size_t start = 0; start + 2 < positive.size(); ++start) {
    std::vector<int> candidate(positive.begin() + start, positive.end());
    Regression fit = log_regression(time, concentration, candidate);
    if (!fit.valid) continue;
    // Prefer additional points when adjusted R2 is practically tied.
    if (fit.adjusted_r2 > best_adjusted + 1e-8 ||
        (std::abs(fit.adjusted_r2 - best_adjusted) <= 1e-8 &&
         candidate.size() > best.size())) {
      best_adjusted = fit.adjusted_r2;
      best = candidate;
    }
  }
  return best;
}

}  // namespace

// [[Rcpp::export(name = ".liberation_nca_profile")]]
Rcpp::List liberation_nca_profile(
    const Rcpp::NumericVector& time_input,
    const Rcpp::NumericVector& concentration_input,
    std::string method, double dose, double tau, std::string route,
    const Rcpp::IntegerVector& terminal_indices,
    const Rcpp::NumericVector& partial_start,
    const Rcpp::NumericVector& partial_end) {
  const int n = time_input.size();
  if (n < 2 || concentration_input.size() != n)
    Rcpp::stop("A profile requires at least two paired observations.");
  std::vector<double> time(n), concentration(n);
  for (int i = 0; i < n; ++i) {
    time[i] = time_input[i]; concentration[i] = concentration_input[i];
    if (!std::isfinite(time[i]) || !std::isfinite(concentration[i]))
      Rcpp::stop("NCA inputs must be finite after preprocessing.");
    if (i > 0 && !(time[i] > time[i - 1]))
      Rcpp::stop("NCA times must be strictly increasing after preprocessing.");
  }
  const bool log_down = method == "lin_up_log_down";
  int maximum_index = static_cast<int>(std::distance(
      concentration.begin(), std::max_element(concentration.begin(), concentration.end())));
  int last_positive = -1;
  for (int i = 0; i < n; ++i) if (concentration[i] > 0.0) last_positive = i;

  double auc_last = 0.0, aumc_last = 0.0;
  if (last_positive > 0) {
    for (int i = 1; i <= last_positive; ++i) {
      IntervalArea area = interval_area(time[i - 1], concentration[i - 1],
                                        time[i], concentration[i], log_down);
      auc_last += area.auc;
      aumc_last += area.aumc;
    }
  } else {
    auc_last = kNa; aumc_last = kNa;
  }
  std::vector<int> terminal;
  if (terminal_indices.size()) {
    terminal.reserve(terminal_indices.size());
    for (int index : terminal_indices) terminal.push_back(index - 1);
  } else {
    terminal = automatic_terminal(time, concentration, maximum_index);
  }
  Regression fit = log_regression(time, concentration, terminal);
  const double lambda_z = fit.valid ? -fit.slope : kNa;
  const double c_last = last_positive >= 0 ? concentration[last_positive] : kNa;
  const double t_last = last_positive >= 0 ? time[last_positive] : kNa;
  const double c_last_pred = fit.valid ? std::exp(fit.intercept + fit.slope * t_last) : kNa;
  const double auc_extra = fit.valid ? c_last / lambda_z : kNa;
  const double auc_inf = fit.valid ? auc_last + auc_extra : kNa;
  const double aumc_extra = fit.valid ? c_last * (t_last / lambda_z + 1.0 / (lambda_z * lambda_z)) : kNa;
  const double aumc_inf = fit.valid ? aumc_last + aumc_extra : kNa;
  const bool has_dose = std::isfinite(dose) && dose > 0.0;
  const bool has_tau = std::isfinite(tau) && tau > 0.0;
  const bool extravascular = route == "extravascular" || route == "oral" || route == "ev";
  const double clearance = has_dose && std::isfinite(auc_inf) ? dose / auc_inf : kNa;
  const double volume_z = std::isfinite(clearance) && std::isfinite(lambda_z) ? clearance / lambda_z : kNa;
  const double mrt = std::isfinite(aumc_inf) && std::isfinite(auc_inf) && auc_inf != 0.0 ? aumc_inf / auc_inf : kNa;
  const double auc_tau = has_tau ? partial_auc(time, concentration, time.front(), time.front() + tau, log_down) : kNa;

  Rcpp::NumericVector partial(partial_start.size(), kNa);
  if (partial_start.size() != partial_end.size())
    Rcpp::stop("Partial-AUC start and end vectors must have equal length.");
  for (R_xlen_t i = 0; i < partial.size(); ++i)
    partial[i] = partial_auc(time, concentration, partial_start[i], partial_end[i], log_down);

  Rcpp::IntegerVector selected(terminal.size());
  for (std::size_t i = 0; i < terminal.size(); ++i) selected[i] = terminal[i] + 1;
  return Rcpp::List::create(
      Rcpp::_["N"] = n,
      Rcpp::_["CMAX"] = concentration[maximum_index],
      Rcpp::_["CMIN"] = *std::min_element(concentration.begin(), concentration.end()),
      Rcpp::_["TMAX"] = time[maximum_index],
      Rcpp::_["CLAST"] = c_last,
      Rcpp::_["CLAST_PRED"] = c_last_pred,
      Rcpp::_["TLAST"] = t_last,
      Rcpp::_["AUCLAST"] = auc_last,
      Rcpp::_["AUMCLAST"] = aumc_last,
      Rcpp::_["LAMBDA_Z"] = lambda_z,
      Rcpp::_["LAMBDA_Z_N"] = static_cast<int>(terminal.size()),
      Rcpp::_["LAMBDA_Z_LOWER"] = terminal.empty() ? kNa : time[terminal.front()],
      Rcpp::_["LAMBDA_Z_UPPER"] = terminal.empty() ? kNa : time[terminal.back()],
      Rcpp::_["R2"] = fit.r2,
      Rcpp::_["R2_ADJ"] = fit.adjusted_r2,
      Rcpp::_["HALF_LIFE"] = fit.valid ? std::log(2.0) / lambda_z : kNa,
      Rcpp::_["AUCINF_OBS"] = auc_inf,
      Rcpp::_["AUC_EXTRAP_PERCENT"] = fit.valid && auc_inf > 0.0 ? 100.0 * auc_extra / auc_inf : kNa,
      Rcpp::_["AUMCINF_OBS"] = aumc_inf,
      Rcpp::_["MRT"] = mrt,
      Rcpp::_[extravascular ? "CL_F" : "CL"] = clearance,
      Rcpp::_[extravascular ? "VZ_F" : "VZ"] = volume_z,
      Rcpp::_["CMAX_DOSE_NORM"] = has_dose ? concentration[maximum_index] / dose : kNa,
      Rcpp::_["AUCINF_DOSE_NORM"] = has_dose && std::isfinite(auc_inf) ? auc_inf / dose : kNa,
      Rcpp::_["AUC_TAU"] = auc_tau,
      Rcpp::_["CAVG"] = has_tau && std::isfinite(auc_tau) ? auc_tau / tau : kNa,
      Rcpp::_["FLUCTUATION_PERCENT"] = has_tau && std::isfinite(auc_tau) && auc_tau != 0.0 ?
        100.0 * (concentration[maximum_index] - c_last) / (auc_tau / tau) : kNa,
      Rcpp::_["partial_auc"] = partial,
      Rcpp::_["terminal_indices"] = selected);
}
