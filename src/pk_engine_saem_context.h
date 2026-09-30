// Internal implementation fragment included inside namespace liberation by
// pk_engine.cpp. It intentionally has no include guard or namespace declaration.
// Contains the persistent fixed-ETA collection used by the compatibility SAEM M-step.

class SaemFixedEtaCollection {
 public:
  SaemFixedEtaCollection(
      const Rcpp::List& tape_pointers, const Rcpp::NumericMatrix& eta,
      int n_theta, int n_sigma, int n_omega)
      : n_theta_(n_theta), n_eta_(eta.ncol()), n_sigma_(n_sigma),
        n_omega_(n_omega) {
    if (tape_pointers.size() != eta.nrow() || tape_pointers.size() < 1 ||
        n_theta < 0 || n_sigma < 0 || n_omega < 0 || eta.ncol() < 0) {
      throw std::invalid_argument(
        "Persistent SAEM fixed-ETA inputs are inconsistent.");
    }
    domain_ = static_cast<std::size_t>(
      n_theta_ + n_eta_ + n_sigma_ + n_omega_);
    tapes_.reserve(static_cast<std::size_t>(tape_pointers.size()));
    points_.resize(static_cast<std::size_t>(tape_pointers.size()));
    for (int subject = 0; subject < tape_pointers.size(); ++subject) {
      Rcpp::XPtr<ObjectiveTape> tape(tape_pointers[subject]);
      if (tape->domain_names.size() != domain_) {
        throw std::invalid_argument(
          "A persistent SAEM objective tape has an inconsistent domain.");
      }
      tapes_.push_back(tape.get());
      std::vector<double>& point = points_[static_cast<std::size_t>(subject)];
      point.assign(domain_, 0.0);
      for (int effect = 0; effect < n_eta_; ++effect) {
        const double value = eta(subject, effect);
        if (!std::isfinite(value)) {
          throw std::invalid_argument("Persistent SAEM ETAs must be finite.");
        }
        point[static_cast<std::size_t>(n_theta_ + effect)] = value;
      }
    }
  }

  Rcpp::List evaluate(
      const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega) {
    validate_parameters(theta, sigma, omega);
    Rcpp::NumericVector values(static_cast<R_xlen_t>(tapes_.size()));
    Rcpp::NumericMatrix gradients(
      static_cast<int>(tapes_.size()), static_cast<int>(domain_));
    const std::vector<double> weight(1U, 1.0);
    std::ostringstream messages;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      std::vector<double>& point = points_[subject];
      std::copy(theta.begin(), theta.end(), point.begin());
      std::copy(
        sigma.begin(), sigma.end(), point.begin() + n_theta_ + n_eta_);
      std::copy(
        omega.begin(), omega.end(),
        point.begin() + n_theta_ + n_eta_ + n_sigma_);
      ObjectiveTape& tape = *tapes_[subject];
      const std::vector<double> value = tape.fun.Forward(0, point, messages);
      require_unchanged_path(tape.fun, "persistent SAEM fixed-ETA objective");
      if (value.size() != 1U || !std::isfinite(value[0])) {
        throw std::domain_error(
          "A persistent SAEM fixed-ETA objective was non-finite.");
      }
      values[static_cast<R_xlen_t>(subject)] = value[0];
      const std::vector<double> derivative = tape.fun.Reverse(1, weight);
      require_unchanged_path(tape.fun, "persistent SAEM fixed-ETA gradient");
      if (derivative.size() != domain_) {
        throw std::runtime_error(
          "A persistent SAEM tape returned an invalid gradient length.");
      }
      for (std::size_t column = 0; column < domain_; ++column) {
        gradients(static_cast<int>(subject), static_cast<int>(column)) =
          derivative[column];
      }
      if ((subject + 1U) % 256U == 0U) Rcpp::checkUserInterrupt();
    }
    ++evaluations_;
    return Rcpp::List::create(
      Rcpp::Named("value") = values,
      Rcpp::Named("gradient") = gradients,
      Rcpp::Named("evaluations") = evaluations_);
  }

  Rcpp::List evaluate_aggregate(
      const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega) {
    validate_parameters(theta, sigma, omega);
    double value_total = 0.0;
    Rcpp::NumericVector gradient_total(static_cast<R_xlen_t>(domain_));
    const std::vector<double> weight(1U, 1.0);
    std::ostringstream messages;
    for (std::size_t subject = 0; subject < tapes_.size(); ++subject) {
      std::vector<double>& point = points_[subject];
      std::copy(theta.begin(), theta.end(), point.begin());
      std::copy(
        sigma.begin(), sigma.end(), point.begin() + n_theta_ + n_eta_);
      std::copy(
        omega.begin(), omega.end(),
        point.begin() + n_theta_ + n_eta_ + n_sigma_);
      ObjectiveTape& tape = *tapes_[subject];
      const std::vector<double> value = tape.fun.Forward(0, point, messages);
      require_unchanged_path(
        tape.fun, "persistent aggregate SAEM fixed-ETA objective");
      if (value.size() != 1U || !std::isfinite(value[0])) {
        throw std::domain_error(
          "A persistent aggregate SAEM fixed-ETA objective was non-finite.");
      }
      value_total += value[0];
      const std::vector<double> derivative = tape.fun.Reverse(1, weight);
      require_unchanged_path(
        tape.fun, "persistent aggregate SAEM fixed-ETA gradient");
      if (derivative.size() != domain_) {
        throw std::runtime_error(
          "A persistent aggregate SAEM tape returned an invalid gradient length.");
      }
      for (std::size_t column = 0; column < domain_; ++column) {
        gradient_total[static_cast<R_xlen_t>(column)] += derivative[column];
      }
      if ((subject + 1U) % 256U == 0U) Rcpp::checkUserInterrupt();
    }
    ++evaluations_;
    if (!tapes_.empty()) {
      gradient_total.attr("names") = Rcpp::wrap(tapes_.front()->domain_names);
    }
    return Rcpp::List::create(
      Rcpp::Named("value") = value_total,
      Rcpp::Named("gradient") = gradient_total,
      Rcpp::Named("evaluations") = evaluations_);
  }

 private:
  int n_theta_ = 0;
  int n_eta_ = 0;
  int n_sigma_ = 0;
  int n_omega_ = 0;
  std::size_t domain_ = 0U;
  int evaluations_ = 0;
  std::vector<ObjectiveTape*> tapes_;
  std::vector<std::vector<double>> points_;

  void validate_parameters(
      const Rcpp::NumericVector& theta, const Rcpp::NumericVector& sigma,
      const Rcpp::NumericVector& omega) const {
    if (theta.size() != n_theta_ || sigma.size() != n_sigma_ ||
        omega.size() != n_omega_) {
      throw std::invalid_argument(
        "Persistent SAEM population parameter dimensions changed.");
    }
    for (double value : theta) if (!std::isfinite(value)) {
      throw std::invalid_argument("Persistent SAEM THETAs must be finite.");
    }
    for (double value : sigma) if (!std::isfinite(value)) {
      throw std::invalid_argument("Persistent SAEM SIGMAs must be finite.");
    }
    for (double value : omega) if (!std::isfinite(value)) {
      throw std::invalid_argument("Persistent SAEM OMEGAs must be finite.");
    }
  }
};

// Persistent complete-data expectation shared by ITS, IMP and SAEM.  ETA
// support points and their normalized weights may be replaced (ITS/IMP) or
// advanced through the exact Robbins--Monro recurrence (SAEM), while the
// subject tapes, dynamic inputs and point buffers remain resident.  All
// reductions deliberately retain subject-then-support order so the
// compatibility path does not acquire scheduler-dependent floating-point
// results.
