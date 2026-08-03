struct HmmDecodeObservation {
  int row = -1;
  std::vector<double> initial;
  std::vector<std::vector<double>> transition;
  std::vector<double> log_emission;
  std::vector<double> filtered;
  double log_scale = 0.0;
};

struct HmmDecodeSequence {
  int subject = 0;
  int sequence = 1;
  std::vector<HmmDecodeObservation> observations;
};

double hmm_log_probability(double value) {
  return value > 0.0 ? std::log(value) :
    -std::numeric_limits<double>::infinity();
}

double hmm_log_sum_exp(const std::vector<double>& values) {
  if (values.empty()) {
    throw std::invalid_argument("HMM log-sum-exp requires values.");
  }
  const double maximum = *std::max_element(values.begin(), values.end());
  if (!std::isfinite(maximum)) return maximum;
  double total = 0.0;
  for (double value : values) total += std::exp(value - maximum);
  return maximum + std::log(total);
}

Rcpp::List hmm_filter(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  if (!engine.hmm_enabled) {
    throw std::invalid_argument("The model does not define a hidden Markov likelihood.");
  }
  if (!engine.mixture_probabilities.empty()) {
    throw std::invalid_argument(
      "Hidden-state filtering for an additional finite-mixture layer is not yet available.");
  }
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector mdv = data["MDV"];
  Rcpp::NumericVector dv = data["DV"];
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  if (eta.nrow() != n_subjects) {
    throw std::invalid_argument("The ETA matrix does not match the HMM dataset subjects.");
  }
  const std::vector<double> theta_values = Rcpp::as<std::vector<double>>(theta);
  std::vector<double> eta_values;
  eta_values.reserve(static_cast<std::size_t>(eta.size()));
  for (int row = 0; row < eta.nrow(); ++row) {
    for (int column = 0; column < eta.ncol(); ++column) {
      eta_values.push_back(eta(row, column));
    }
  }
  const std::vector<double> sigma_values = Rcpp::as<std::vector<double>>(sigma);
  const std::vector<double> prediction = simulate_analytical_t(
    engine, data, theta_values, eta_values, sigma_values);
  const int rows = data.nrows();
  const int states = engine.hmm_states;
  Rcpp::NumericMatrix probability(rows, states);
  Rcpp::NumericMatrix smoothed(rows, states);
  std::fill(probability.begin(), probability.end(), NA_REAL);
  std::fill(smoothed.begin(), smoothed.end(), NA_REAL);
  Rcpp::IntegerVector filtered_state(rows, NA_INTEGER);
  Rcpp::IntegerVector smoothed_state(rows, NA_INTEGER);
  Rcpp::IntegerVector viterbi_state(rows, NA_INTEGER);
  Rcpp::NumericVector row_nll(rows, NA_REAL);
  const bool has_dvid = data.containsElementNamed("DVID");
  int previous_subject = -1;
  std::unordered_map<int, double> previous_outcome;
  std::unordered_map<int, double> previous_outcome_time;
  std::unordered_map<int, std::vector<double>> filtered;
  std::unordered_map<int, std::size_t> sequence_lookup;
  std::vector<HmmDecodeSequence> sequences;
  double total_nll = 0.0;
  for (int row = 0; row < rows; ++row) {
    const int subject = subjects[row] - 1;
    if (subject != previous_subject) {
      previous_outcome.clear();
      previous_outcome_time.clear();
      filtered.clear();
      sequence_lookup.clear();
    }
    const int dvid = has_dvid ?
      static_cast<int>(row_optional(data, "DVID", row, 1.0)) : 1;
    const int sequence = engine.hmm_by_dvid ? dvid : 1;
    if (evid[row] == 0 && mdv[row] == 0 && std::isfinite(dv[row])) {
      const auto previous = previous_outcome.find(sequence);
      const bool first = previous == previous_outcome.end();
      const double previous_value = first ? dv[row] : previous->second;
      const double previous_time = first ? row_optional(data, "TIME", row, 0.0) :
        previous_outcome_time.at(sequence);
      const int mixture_number = static_cast<int>(
        row_optional(data, "MIXNUM", row, 1.0));
      HmmRowComponents<double> components;
      const double contribution = hmm_row_nll_t(
        engine, data, row, subject, theta_values, eta_values, eta.ncol(),
        sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
        first, previous_value, previous_time, filtered[sequence], &components);
      row_nll[row] = contribution;
      total_nll += contribution;
      int most_likely = 0;
      for (int hidden = 0; hidden < states; ++hidden) {
        const double value = filtered[sequence][static_cast<std::size_t>(hidden)];
        probability(row, hidden) = value;
        if (value > filtered[sequence][static_cast<std::size_t>(most_likely)]) {
          most_likely = hidden;
        }
      }
      filtered_state[row] = most_likely + 1;
      if (first) {
        sequence_lookup[sequence] = sequences.size();
        HmmDecodeSequence track;
        track.subject = subject + 1;
        track.sequence = sequence;
        sequences.push_back(std::move(track));
      }
      const auto track = sequence_lookup.find(sequence);
      if (track == sequence_lookup.end()) {
        throw std::logic_error("Hidden Markov decoding sequence was not initialized.");
      }
      HmmDecodeObservation observation;
      observation.row = row;
      observation.initial = std::move(components.initial);
      observation.transition = std::move(components.transition);
      observation.log_emission = std::move(components.log_emission);
      if (engine.hmm_emission_scale != "log") {
        for (int hidden = 0; hidden < states; ++hidden) {
          if (!(components.emission[static_cast<std::size_t>(hidden)] > 0.0)) {
            observation.log_emission[static_cast<std::size_t>(hidden)] =
              -std::numeric_limits<double>::infinity();
          }
        }
      }
      observation.filtered = filtered[sequence];
      observation.log_scale = -0.5 * contribution;
      sequences[track->second].observations.push_back(std::move(observation));
      previous_outcome[sequence] = dv[row];
      previous_outcome_time[sequence] = row_optional(data, "TIME", row, 0.0);
    }
    previous_subject = subject;
  }

  Rcpp::IntegerVector summary_subject(static_cast<R_xlen_t>(sequences.size()));
  Rcpp::IntegerVector summary_sequence(static_cast<R_xlen_t>(sequences.size()));
  Rcpp::IntegerVector summary_observations(static_cast<R_xlen_t>(sequences.size()));
  Rcpp::NumericVector summary_log_likelihood(static_cast<R_xlen_t>(sequences.size()));
  Rcpp::NumericVector summary_viterbi_joint(static_cast<R_xlen_t>(sequences.size()));
  Rcpp::NumericVector summary_viterbi_posterior(static_cast<R_xlen_t>(sequences.size()));

  for (std::size_t sequence_index = 0; sequence_index < sequences.size();
       ++sequence_index) {
    const HmmDecodeSequence& sequence = sequences[sequence_index];
    const std::size_t length = sequence.observations.size();
    if (length == 0U) continue;
    const std::size_t state_count = static_cast<std::size_t>(states);

    std::vector<std::vector<double>> log_beta(
      length, std::vector<double>(state_count, 0.0));
    for (std::size_t current = length - 1U; current > 0U; --current) {
      const HmmDecodeObservation& next = sequence.observations[current];
      for (std::size_t from = 0; from < state_count; ++from) {
        std::vector<double> terms(state_count);
        for (std::size_t to = 0; to < state_count; ++to) {
          terms[to] = hmm_log_probability(next.transition[from][to]) +
            next.log_emission[to] + log_beta[current][to];
        }
        log_beta[current - 1U][from] =
          hmm_log_sum_exp(terms) - next.log_scale;
      }
    }
    for (std::size_t current = 0; current < length; ++current) {
      const HmmDecodeObservation& observation = sequence.observations[current];
      std::vector<double> log_probability(state_count);
      for (std::size_t hidden = 0; hidden < state_count; ++hidden) {
        log_probability[hidden] = hmm_log_probability(observation.filtered[hidden]) +
          log_beta[current][hidden];
      }
      const double normalizer = hmm_log_sum_exp(log_probability);
      int most_likely = 0;
      for (std::size_t hidden = 0; hidden < state_count; ++hidden) {
        const double value = std::exp(log_probability[hidden] - normalizer);
        smoothed(observation.row, static_cast<int>(hidden)) = value;
        if (value > smoothed(observation.row, most_likely)) {
          most_likely = static_cast<int>(hidden);
        }
      }
      smoothed_state[observation.row] = most_likely + 1;
    }

    std::vector<std::vector<double>> delta(
      length, std::vector<double>(state_count,
        -std::numeric_limits<double>::infinity()));
    std::vector<std::vector<int>> back_pointer(
      length, std::vector<int>(state_count, -1));
    const HmmDecodeObservation& first = sequence.observations.front();
    for (std::size_t hidden = 0; hidden < state_count; ++hidden) {
      delta[0][hidden] = hmm_log_probability(first.initial[hidden]) +
        first.log_emission[hidden];
    }
    for (std::size_t current = 1; current < length; ++current) {
      const HmmDecodeObservation& observation = sequence.observations[current];
      for (std::size_t to = 0; to < state_count; ++to) {
        int best_from = 0;
        double best = delta[current - 1U][0] +
          hmm_log_probability(observation.transition[0][to]);
        for (std::size_t from = 1; from < state_count; ++from) {
          const double candidate = delta[current - 1U][from] +
            hmm_log_probability(observation.transition[from][to]);
          if (candidate > best) {
            best = candidate;
            best_from = static_cast<int>(from);
          }
        }
        delta[current][to] = best + observation.log_emission[to];
        back_pointer[current][to] = best_from;
      }
    }
    int final_state = 0;
    for (int hidden = 1; hidden < states; ++hidden) {
      if (delta[length - 1U][static_cast<std::size_t>(hidden)] >
          delta[length - 1U][static_cast<std::size_t>(final_state)]) {
        final_state = hidden;
      }
    }
    std::vector<int> path(length, final_state);
    for (std::size_t current = length - 1U; current > 0U; --current) {
      path[current - 1U] = back_pointer[current]
        [static_cast<std::size_t>(path[current])];
    }
    for (std::size_t current = 0; current < length; ++current) {
      viterbi_state[sequence.observations[current].row] = path[current] + 1;
    }

    double sequence_log_likelihood = 0.0;
    for (const HmmDecodeObservation& observation : sequence.observations) {
      sequence_log_likelihood += observation.log_scale;
    }
    const double viterbi_log_joint =
      delta[length - 1U][static_cast<std::size_t>(final_state)];
    summary_subject[static_cast<R_xlen_t>(sequence_index)] = sequence.subject;
    summary_sequence[static_cast<R_xlen_t>(sequence_index)] = sequence.sequence;
    summary_observations[static_cast<R_xlen_t>(sequence_index)] =
      static_cast<int>(length);
    summary_log_likelihood[static_cast<R_xlen_t>(sequence_index)] =
      sequence_log_likelihood;
    summary_viterbi_joint[static_cast<R_xlen_t>(sequence_index)] =
      viterbi_log_joint;
    summary_viterbi_posterior[static_cast<R_xlen_t>(sequence_index)] =
      viterbi_log_joint - sequence_log_likelihood;
  }

  probability.attr("dimnames") = Rcpp::List::create(
    R_NilValue, Rcpp::wrap(engine.hmm_state_names));
  smoothed.attr("dimnames") = Rcpp::List::create(
    R_NilValue, Rcpp::wrap(engine.hmm_state_names));
  auto state_labels = [&](const Rcpp::IntegerVector& state) {
    Rcpp::CharacterVector labels(rows, NA_STRING);
    for (int row = 0; row < rows; ++row) {
      if (state[row] != NA_INTEGER) {
        labels[row] = engine.hmm_state_names[
          static_cast<std::size_t>(state[row] - 1)];
      }
    }
    return labels;
  };
  const Rcpp::CharacterVector filtered_label = state_labels(filtered_state);
  const Rcpp::CharacterVector smoothed_label = state_labels(smoothed_state);
  const Rcpp::CharacterVector viterbi_label = state_labels(viterbi_state);
  const Rcpp::DataFrame sequence_summary = Rcpp::DataFrame::create(
    Rcpp::Named(".ID_INDEX") = summary_subject,
    Rcpp::Named("HMM_SEQUENCE") = summary_sequence,
    Rcpp::Named("OBSERVATIONS") = summary_observations,
    Rcpp::Named("LOG_LIKELIHOOD") = summary_log_likelihood,
    Rcpp::Named("VITERBI_LOG_JOINT") = summary_viterbi_joint,
    Rcpp::Named("VITERBI_LOG_POSTERIOR") = summary_viterbi_posterior);
  return Rcpp::List::create(
    Rcpp::Named("filtered") = probability,
    Rcpp::Named("smoothed") = smoothed,
    Rcpp::Named("state") = filtered_state,
    Rcpp::Named("state_label") = filtered_label,
    Rcpp::Named("filtered_state") = filtered_state,
    Rcpp::Named("filtered_state_label") = filtered_label,
    Rcpp::Named("smoothed_state") = smoothed_state,
    Rcpp::Named("smoothed_state_label") = smoothed_label,
    Rcpp::Named("viterbi_state") = viterbi_state,
    Rcpp::Named("viterbi_state_label") = viterbi_label,
    Rcpp::Named("row_nll") = row_nll,
    Rcpp::Named("nll") = total_nll,
    Rcpp::Named("log_likelihood") = -0.5 * total_nll,
    Rcpp::Named("sequence_summary") = sequence_summary,
    Rcpp::Named("states") = engine.hmm_state_names);
}

struct KalmanDecodeObservation {
  int row = -1;
  Vector predicted_mean;
  Matrix predicted_covariance;
  Vector filtered_mean;
  Matrix filtered_covariance;
  Matrix transition;
  Matrix smoother_cross_covariance;
  std::vector<Vector> particle_states;
  std::vector<double> particle_weights;
  std::vector<int> particle_ancestors;
  std::vector<int> particle_regimes;
  double innovation = NA_REAL;
  double innovation_variance = NA_REAL;
  double row_nll = NA_REAL;
};

struct KalmanDecodeSequence {
  std::vector<KalmanDecodeObservation> observations;
};

Rcpp::List kalman_filter(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma) {
  if (!engine.kalman_enabled) {
    throw std::invalid_argument("The model does not define a Kalman state-space likelihood.");
  }
  if (!engine.mixture_probabilities.empty()) {
    throw std::invalid_argument(
      "Kalman decoding for an additional finite-mixture layer is not yet available.");
  }
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector mdv = data["MDV"];
  Rcpp::NumericVector dv = data["DV"];
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  if (eta.nrow() != n_subjects) {
    throw std::invalid_argument("The ETA matrix does not match the Kalman dataset subjects.");
  }
  const std::vector<double> theta_values = Rcpp::as<std::vector<double>>(theta);
  std::vector<double> eta_values;
  eta_values.reserve(static_cast<std::size_t>(eta.size()));
  for (int row = 0; row < eta.nrow(); ++row) {
    for (int column = 0; column < eta.ncol(); ++column) {
      eta_values.push_back(eta(row, column));
    }
  }
  const std::vector<double> sigma_values = Rcpp::as<std::vector<double>>(sigma);
  const std::vector<double> prediction = simulate_analytical_t(
    engine, data, theta_values, eta_values, sigma_values);
  const int rows = data.nrows();
  const int states = engine.kalman_states;
  Rcpp::NumericMatrix predicted_mean(rows, states);
  Rcpp::NumericMatrix filtered_mean(rows, states);
  Rcpp::NumericMatrix smoothed_mean(rows, states);
  Rcpp::NumericMatrix filtered_variance(rows, states);
  Rcpp::NumericMatrix smoothed_variance(rows, states);
  std::fill(predicted_mean.begin(), predicted_mean.end(), NA_REAL);
  std::fill(filtered_mean.begin(), filtered_mean.end(), NA_REAL);
  std::fill(smoothed_mean.begin(), smoothed_mean.end(), NA_REAL);
  std::fill(filtered_variance.begin(), filtered_variance.end(), NA_REAL);
  std::fill(smoothed_variance.begin(), smoothed_variance.end(), NA_REAL);
  Rcpp::NumericVector innovation(rows, NA_REAL);
  Rcpp::NumericVector innovation_variance(rows, NA_REAL);
  Rcpp::NumericVector row_nll(rows, NA_REAL);
  Rcpp::NumericMatrix filtered_regime(
    rows, engine.switching_enabled ? engine.switching_regimes : 0);
  Rcpp::NumericMatrix smoothed_regime(
    rows, engine.switching_enabled ? engine.switching_regimes : 0);
  std::fill(filtered_regime.begin(), filtered_regime.end(), NA_REAL);
  std::fill(smoothed_regime.begin(), smoothed_regime.end(), NA_REAL);
  const bool has_dvid = data.containsElementNamed("DVID");
  int previous_subject = -1;
  std::unordered_map<int, double> previous_outcome;
  std::unordered_map<int, double> previous_outcome_time;
  std::unordered_map<int, KalmanFilterState<double>> filtered;
  std::unordered_map<int, ParticleFilterState<double>> particle_filter;
  std::unordered_map<int, std::size_t> sequence_lookup;
  std::vector<KalmanDecodeSequence> sequences;
  double total_nll = 0.0;
  for (int row = 0; row < rows; ++row) {
    const int subject = subjects[row] - 1;
    if (subject != previous_subject) {
      previous_outcome.clear();
      previous_outcome_time.clear();
      filtered.clear();
      particle_filter.clear();
      sequence_lookup.clear();
    }
    const int dvid = has_dvid ?
      static_cast<int>(row_optional(data, "DVID", row, 1.0)) : 1;
    const int sequence = engine.kalman_by_dvid ? dvid : 1;
    if (evid[row] == 0 && mdv[row] == 0 && std::isfinite(dv[row])) {
      const auto previous = previous_outcome.find(sequence);
      const bool first = previous == previous_outcome.end();
      const double previous_value = first ? dv[row] : previous->second;
      const double previous_time = first ? row_optional(data, "TIME", row, 0.0) :
        previous_outcome_time.at(sequence);
      const int mixture_number = static_cast<int>(
        row_optional(data, "MIXNUM", row, 1.0));
      KalmanRowComponents<double> components;
      const double contribution = engine.kalman_filter_type == "particle" ?
        particle_row_nll_t(
          engine, data, row, subject, theta_values, eta_values, eta.ncol(),
          sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
          first, previous_value, previous_time, particle_filter[sequence], &components) :
        kalman_row_nll_t(
          engine, data, row, subject, theta_values, eta_values, eta.ncol(),
          sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
          first, previous_value, previous_time, filtered[sequence], &components);
      if (first) {
        sequence_lookup[sequence] = sequences.size();
        sequences.emplace_back();
      }
      const auto track = sequence_lookup.find(sequence);
      if (track == sequence_lookup.end()) {
        throw std::logic_error("Kalman decoding sequence was not initialized.");
      }
      KalmanDecodeObservation observation;
      observation.row = row;
      observation.predicted_mean = components.predicted_mean;
      observation.predicted_covariance = components.predicted_covariance;
      observation.filtered_mean = components.filtered_mean;
      observation.filtered_covariance = components.filtered_covariance;
      observation.transition = components.transition;
      observation.smoother_cross_covariance = components.smoother_cross_covariance;
      observation.particle_states = components.particle_states;
      observation.particle_weights = components.particle_weights;
      observation.particle_ancestors = components.particle_ancestors;
      observation.particle_regimes = components.particle_regimes;
      observation.innovation = components.innovation;
      observation.innovation_variance = components.innovation_variance;
      observation.row_nll = contribution;
      sequences[track->second].observations.push_back(std::move(observation));
      for (int state = 0; state < states; ++state) {
        predicted_mean(row, state) = components.predicted_mean[state];
        filtered_mean(row, state) = components.filtered_mean[state];
        filtered_variance(row, state) = components.filtered_covariance(state, state);
      }
      if (engine.switching_enabled) {
        for (int regime = 0; regime < engine.switching_regimes; ++regime) {
          double probability = 0.0;
          for (std::size_t particle = 0;
               particle < components.particle_regimes.size(); ++particle) {
            if (components.particle_regimes[particle] == regime &&
                particle < components.particle_weights.size()) {
              probability += components.particle_weights[particle];
            }
          }
          filtered_regime(row, regime) = probability;
        }
      }
      innovation[row] = components.innovation;
      innovation_variance[row] = components.innovation_variance;
      row_nll[row] = contribution;
      total_nll += contribution;
      previous_outcome[sequence] = dv[row];
      previous_outcome_time[sequence] = row_optional(data, "TIME", row, 0.0);
    }
    previous_subject = subject;
  }
  for (const KalmanDecodeSequence& sequence : sequences) {
    const std::size_t length = sequence.observations.size();
    if (!length) continue;
    std::vector<Vector> smooth_mean(length);
    std::vector<Matrix> smooth_covariance(length);
    if (engine.kalman_filter_type == "particle") {
      std::vector<std::vector<double>> smooth_weight(length);
      smooth_weight[length - 1U] = sequence.observations[length - 1U].particle_weights;
      for (std::size_t offset = 1U; offset < length; ++offset) {
        const std::size_t current = length - offset - 1U;
        const KalmanDecodeObservation& next = sequence.observations[current + 1U];
        const std::size_t particles = sequence.observations[current].particle_states.size();
        smooth_weight[current].assign(particles, 0.0);
        for (std::size_t particle = 0; particle < next.particle_ancestors.size(); ++particle) {
          const int ancestor = next.particle_ancestors[particle];
          if (ancestor >= 0 && static_cast<std::size_t>(ancestor) < particles &&
              particle < smooth_weight[current + 1U].size()) {
            smooth_weight[current][static_cast<std::size_t>(ancestor)] +=
              smooth_weight[current + 1U][particle];
          }
        }
        const double total = std::accumulate(
          smooth_weight[current].begin(), smooth_weight[current].end(), 0.0);
        if (total > 0.0) for (double& weight : smooth_weight[current]) weight /= total;
      }
      for (std::size_t current = 0; current < length; ++current) {
        const auto& cloud = sequence.observations[current].particle_states;
        const auto& weights = smooth_weight[current];
        smooth_mean[current] = Vector::Zero(states);
        for (std::size_t particle = 0; particle < cloud.size(); ++particle) {
          smooth_mean[current] += weights[particle] * cloud[particle];
        }
        smooth_covariance[current] = Matrix::Zero(states, states);
        for (std::size_t particle = 0; particle < cloud.size(); ++particle) {
          const Vector centered = cloud[particle] - smooth_mean[current];
          smooth_covariance[current] += weights[particle] * centered * centered.transpose();
        }
        if (engine.switching_enabled) {
          const int row = sequence.observations[current].row;
          const auto& regimes = sequence.observations[current].particle_regimes;
          for (int regime = 0; regime < engine.switching_regimes; ++regime) {
            double probability = 0.0;
            for (std::size_t particle = 0; particle < regimes.size() &&
                 particle < weights.size(); ++particle) {
              if (regimes[particle] == regime) probability += weights[particle];
            }
            smoothed_regime(row, regime) = probability;
          }
        }
      }
    } else {
      smooth_mean[length - 1U] = sequence.observations[length - 1U].filtered_mean;
      smooth_covariance[length - 1U] =
        sequence.observations[length - 1U].filtered_covariance;
      for (std::size_t offset = 1U; offset < length; ++offset) {
        const std::size_t current = length - offset - 1U;
        const KalmanDecodeObservation& observation = sequence.observations[current];
        const KalmanDecodeObservation& next = sequence.observations[current + 1U];
        Eigen::LDLT<Matrix> decomposition(next.predicted_covariance);
        if (decomposition.info() != Eigen::Success) {
          throw std::runtime_error("RTS smoother predicted covariance factorization failed.");
        }
        const Matrix smoother_gain = next.smoother_cross_covariance *
          decomposition.solve(Matrix::Identity(states, states));
        smooth_mean[current] = observation.filtered_mean + smoother_gain *
          (smooth_mean[current + 1U] - next.predicted_mean);
        smooth_covariance[current] = observation.filtered_covariance + smoother_gain *
          (smooth_covariance[current + 1U] - next.predicted_covariance) *
            smoother_gain.transpose();
        smooth_covariance[current] = 0.5 *
          (smooth_covariance[current] + smooth_covariance[current].transpose()).eval();
      }
    }
    for (std::size_t current = 0; current < length; ++current) {
      const int row = sequence.observations[current].row;
      for (int state = 0; state < states; ++state) {
        smoothed_mean(row, state) = smooth_mean[current][state];
        smoothed_variance(row, state) = smooth_covariance[current](state, state);
      }
    }
  }
  const Rcpp::List dimnames = Rcpp::List::create(
    R_NilValue, Rcpp::wrap(engine.kalman_state_names));
  predicted_mean.attr("dimnames") = dimnames;
  filtered_mean.attr("dimnames") = dimnames;
  smoothed_mean.attr("dimnames") = dimnames;
  filtered_variance.attr("dimnames") = dimnames;
  smoothed_variance.attr("dimnames") = dimnames;
  if (engine.switching_enabled) {
    const Rcpp::List regime_dimnames = Rcpp::List::create(
      R_NilValue, Rcpp::wrap(engine.switching_regime_names));
    filtered_regime.attr("dimnames") = regime_dimnames;
    smoothed_regime.attr("dimnames") = regime_dimnames;
  }
  return Rcpp::List::create(
    Rcpp::Named("predicted_mean") = predicted_mean,
    Rcpp::Named("filtered_mean") = filtered_mean,
    Rcpp::Named("smoothed_mean") = smoothed_mean,
    Rcpp::Named("filtered_variance") = filtered_variance,
    Rcpp::Named("smoothed_variance") = smoothed_variance,
    Rcpp::Named("innovation") = innovation,
    Rcpp::Named("innovation_variance") = innovation_variance,
    Rcpp::Named("row_nll") = row_nll,
    Rcpp::Named("filtered_regime") = filtered_regime,
    Rcpp::Named("smoothed_regime") = smoothed_regime,
    Rcpp::Named("regimes") = engine.switching_regime_names,
    Rcpp::Named("nll") = total_nll,
    Rcpp::Named("log_likelihood") = -0.5 * total_nll,
    Rcpp::Named("states") = engine.kalman_state_names,
    Rcpp::Named("filter") = engine.kalman_filter_type,
    Rcpp::Named("smoother") = engine.kalman_filter_type == "particle" ?
      "genealogical" : "RTS");
}

Matrix covariance_sampling_root(const Matrix& covariance,
                                const std::string& context) {
  const Matrix symmetric = 0.5 * (covariance + covariance.transpose()).eval();
  auto eigen = libertad::detail::self_adjoint_eigen(symmetric, true);
  if (eigen.info != Eigen::Success) {
    throw std::runtime_error(context + " eigen decomposition failed.");
  }
  const double scale = std::max(eigen.values.cwiseAbs().maxCoeff(), 1.0);
  if (eigen.values.minCoeff() < -1e-10 * scale) {
    throw std::domain_error(context + " is not positive semidefinite.");
  }
  return eigen.vectors * eigen.values.cwiseMax(0.0).cwiseSqrt().asDiagonal();
}

Rcpp::NumericVector kalman_simulate(
    const ModelEngine& engine, const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& theta, const Rcpp::NumericMatrix& eta,
    const Rcpp::NumericVector& sigma,
    const Rcpp::NumericMatrix& process_normals,
    const Rcpp::NumericVector& observation_normals) {
  if (!engine.kalman_enabled) {
    throw std::invalid_argument("The model does not define a Kalman state-space likelihood.");
  }
  const int rows = data.nrows();
  const int states = engine.kalman_states;
  const int process_columns = states *
    (engine.kalman_dynamics == "sde" ? engine.kalman_sde_substeps : 1);
  if (process_normals.nrow() != rows || process_normals.ncol() != process_columns ||
      observation_normals.size() != rows) {
    throw std::invalid_argument("Kalman simulation normal draws have inconsistent dimensions.");
  }
  Rcpp::IntegerVector subjects = data[".ID_INDEX"];
  Rcpp::IntegerVector evid = data["EVID"];
  Rcpp::IntegerVector mdv = data["MDV"];
  Rcpp::NumericVector input_dv = data["DV"];
  Rcpp::NumericVector simulated = Rcpp::clone(input_dv);
  int n_subjects = 0;
  for (int value : subjects) n_subjects = std::max(n_subjects, value);
  if (eta.nrow() != n_subjects) {
    throw std::invalid_argument("The ETA matrix does not match the Kalman simulation subjects.");
  }
  const std::vector<double> theta_values = Rcpp::as<std::vector<double>>(theta);
  std::vector<double> eta_values;
  eta_values.reserve(static_cast<std::size_t>(eta.size()));
  for (int row = 0; row < eta.nrow(); ++row) {
    for (int column = 0; column < eta.ncol(); ++column) eta_values.push_back(eta(row, column));
  }
  const std::vector<double> sigma_values = Rcpp::as<std::vector<double>>(sigma);
  const std::vector<double> prediction = simulate_analytical_t(
    engine, data, theta_values, eta_values, sigma_values);
  const bool has_dvid = data.containsElementNamed("DVID");
  int previous_subject = -1;
  std::unordered_map<int, Vector> latent;
  std::unordered_map<int, int> latent_regime;
  std::unordered_map<int, double> previous_time;
  std::unordered_map<int, double> previous_dv;
  for (int row = 0; row < rows; ++row) {
    const int subject = subjects[row] - 1;
    if (subject != previous_subject) {
      latent.clear();
      latent_regime.clear();
      previous_time.clear();
      previous_dv.clear();
    }
    if (evid[row] != 0 || mdv[row] != 0) {
      previous_subject = subject;
      continue;
    }
    const int dvid = has_dvid ?
      static_cast<int>(row_optional(data, "DVID", row, 1.0)) : 1;
    const int sequence = engine.kalman_by_dvid ? dvid : 1;
    const bool first = latent.find(sequence) == latent.end();
    const double previous_value = first ? 0.0 : previous_dv.at(sequence);
    const double previous_observation_time = first ?
      row_optional(data, "TIME", row, 0.0) : previous_time.at(sequence);
    const int mixture_number = static_cast<int>(row_optional(data, "MIXNUM", row, 1.0));
    Vector normal(states);
    for (int state = 0; state < states; ++state) normal[state] = process_normals(row, state);
    if (engine.kalman_filter_type != "linear") {
      const Vector evaluation_state = first ? Vector::Zero(states) : latent[sequence];
      NonlinearStateOutputs<double> output = nonlinear_state_outputs_t(
        engine, data, row, subject, theta_values, eta_values, eta.ncol(),
        sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
        first, previous_value, previous_observation_time, evaluation_state);
      if (engine.switching_enabled) {
        const SwitchingStateOutputs<double> switching = switching_state_raw_outputs_t(
          engine, data, row, subject, theta_values, eta_values, eta.ncol(),
          sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
          first, previous_value, previous_observation_time, evaluation_state);
        std::vector<double> probability;
        if (first) {
          probability = switching.initial;
        } else {
          probability.resize(static_cast<std::size_t>(engine.switching_regimes));
          const int from = latent_regime.at(sequence);
          for (int to = 0; to < engine.switching_regimes; ++to) {
            probability[static_cast<std::size_t>(to)] = switching.regime_transition(from, to);
          }
        }
        latent_regime[sequence] = sample_switching_regime(
          probability, state_space_uniform(
            engine.kalman_seed + 65537 * subject, row, sequence, 0, first ? 31 : 32));
        const NonlinearStateOutputs<double>& local =
          switching.regime[static_cast<std::size_t>(latent_regime[sequence])];
        output.transition = local.transition;
        output.process_covariance = local.process_covariance;
        output.observation = local.observation;
        output.observation_variance = local.observation_variance;
      }
      if (!(output.observation_variance >= 0.0)) {
        throw std::domain_error("Nonlinear state-space simulation observation variance must be non-negative.");
      }
      if (first) {
        latent[sequence] = output.initial_mean + covariance_sampling_root(
          output.initial_covariance, "State-space initial covariance") * normal;
      } else if (engine.kalman_dynamics == "sde") {
        Vector current = latent[sequence];
        const double current_time = row_optional(data, "TIME", row, previous_observation_time);
        const double interval = std::max(0.0, current_time - previous_observation_time);
        const double step = interval / static_cast<double>(engine.kalman_sde_substeps);
        for (int substep = 0; substep < engine.kalman_sde_substeps; ++substep) {
          const NonlinearStateOutputs<double> local = engine.switching_enabled ?
            switching_state_raw_outputs_t(
              engine, data, row, subject, theta_values, eta_values, eta.ncol(),
              sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
              false, previous_value, previous_observation_time + substep * step,
              current).regime[static_cast<std::size_t>(latent_regime[sequence])] :
            nonlinear_state_raw_outputs_t(
              engine, data, row, subject, theta_values, eta_values, eta.ncol(),
              sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
              false, previous_value, previous_observation_time + substep * step, current);
          Vector increment(states);
          for (int state = 0; state < states; ++state) {
            increment[state] = std::sqrt(step) *
              process_normals(row, substep * states + state);
          }
          Vector next = current + step * local.transition +
            local.process_covariance * increment;
          if (engine.kalman_sde_method == "milstein") {
            for (int first_state = 0; first_state < states; ++first_state) {
              for (int second_state = 0; second_state < states; ++second_state) {
                if (first_state != second_state &&
                    std::abs(local.process_covariance(first_state, second_state)) > 1e-12) {
                  throw std::domain_error(
                    "Milstein simulation currently requires diagonal diffusion.");
                }
              }
              const double derivative_step = engine.kalman_jacobian_step *
                std::max(1.0, std::abs(current[first_state]));
              Vector plus = current;
              Vector minus = current;
              plus[first_state] += derivative_step;
              minus[first_state] -= derivative_step;
              const double upper = engine.switching_enabled ?
                switching_state_raw_outputs_t(
                  engine, data, row, subject, theta_values, eta_values, eta.ncol(),
                  sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
                  false, previous_value, previous_observation_time + substep * step,
                  plus).regime[static_cast<std::size_t>(latent_regime[sequence])].
                    process_covariance(first_state, first_state) :
                nonlinear_state_raw_outputs_t(
                  engine, data, row, subject, theta_values, eta_values, eta.ncol(),
                  sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
                  false, previous_value, previous_observation_time + substep * step,
                  plus).process_covariance(first_state, first_state);
              const double lower = engine.switching_enabled ?
                switching_state_raw_outputs_t(
                  engine, data, row, subject, theta_values, eta_values, eta.ncol(),
                  sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
                  false, previous_value, previous_observation_time + substep * step,
                  minus).regime[static_cast<std::size_t>(latent_regime[sequence])].
                    process_covariance(first_state, first_state) :
                nonlinear_state_raw_outputs_t(
                  engine, data, row, subject, theta_values, eta_values, eta.ncol(),
                  sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
                  false, previous_value, previous_observation_time + substep * step,
                  minus).process_covariance(first_state, first_state);
              const double derivative = (upper - lower) / (2.0 * derivative_step);
              const double diffusion = local.process_covariance(first_state, first_state);
              next[first_state] += 0.5 * diffusion * derivative *
                (increment[first_state] * increment[first_state] - step);
            }
          }
          current = next;
        }
        latent[sequence] = current;
      } else {
        latent[sequence] = output.transition + covariance_sampling_root(
          output.process_covariance, "State-space process covariance") * normal;
      }
      const NonlinearStateOutputs<double> observed = engine.switching_enabled ?
        switching_state_raw_outputs_t(
          engine, data, row, subject, theta_values, eta_values, eta.ncol(),
          sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
          first, previous_value, previous_observation_time,
          latent[sequence]).regime[static_cast<std::size_t>(latent_regime[sequence])] :
        nonlinear_state_outputs_t(
          engine, data, row, subject, theta_values, eta_values, eta.ncol(),
          sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
          first, previous_value, previous_observation_time, latent[sequence]);
      const double baseline = engine.kalman_prediction_baseline ?
        prediction[static_cast<std::size_t>(row)] : 0.0;
      simulated[row] = baseline + observed.observation +
        std::sqrt(observed.observation_variance) * observation_normals[row];
      previous_time[sequence] = row_optional(data, "TIME", row, 0.0);
      previous_dv[sequence] = simulated[row];
      previous_subject = subject;
      continue;
    }
    const std::vector<double> output = evaluate_error_outputs_t(
      engine, data, row, subject, theta_values, eta_values, eta.ncol(),
      sigma_values, mixture_number, prediction[static_cast<std::size_t>(row)],
      first, previous_value, previous_observation_time, engine.kalman_outputs,
      "Kalman state-space simulation");
    std::size_t cursor = 0U;
    Vector initial_mean(states);
    for (int state = 0; state < states; ++state) initial_mean[state] = output[cursor++];
    auto read_matrix = [&]() {
      Matrix matrix(states, states);
      for (int matrix_row = 0; matrix_row < states; ++matrix_row) {
        for (int matrix_column = 0; matrix_column < states; ++matrix_column) {
          matrix(matrix_row, matrix_column) = output[cursor++];
        }
      }
      return matrix;
    };
    const Matrix initial_covariance = read_matrix();
    const Matrix transition = read_matrix();
    const Matrix process_covariance = read_matrix();
    Vector observation(states);
    for (int state = 0; state < states; ++state) observation[state] = output[cursor++];
    const double observation_variance = output[cursor++];
    if (!(observation_variance >= 0.0)) {
      throw std::domain_error("Kalman simulation observation variance must be non-negative.");
    }
    if (first) {
      latent[sequence] = initial_mean + covariance_sampling_root(
        initial_covariance, "Kalman initial covariance") * normal;
    } else {
      latent[sequence] = transition * latent[sequence] + covariance_sampling_root(
        process_covariance, "Kalman process covariance") * normal;
    }
    const double baseline = engine.kalman_prediction_baseline ?
      prediction[static_cast<std::size_t>(row)] : 0.0;
    simulated[row] = baseline + observation.dot(latent[sequence]) +
      std::sqrt(observation_variance) * observation_normals[row];
    previous_time[sequence] = row_optional(data, "TIME", row, 0.0);
    previous_dv[sequence] = simulated[row];
    previous_subject = subject;
  }
  return simulated;
}

#include "hmc_sampler.h"

SEXP population_objective_create_api(
    SEXP engine_pointer, const Rcpp::List& subject_data,
    const Rcpp::List& primary_tape_pointers,
    const Rcpp::List& curvature_tape_pointers,
    const Rcpp::List& config) {
  Rcpp::XPtr<PopulationObjective> pointer(
    new PopulationObjective(
      engine_pointer, subject_data, primary_tape_pointers,
      curvature_tape_pointers, config),
    true
  );
  pointer.attr("class") = Rcpp::CharacterVector::create(
    "liberation_population_objective_ptr", "externalptr");
  pointer.attr("keepers") = Rcpp::List::create(
    engine_pointer, subject_data, primary_tape_pointers,
    curvature_tape_pointers);
  return pointer;
}

double population_objective_value_api(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  Rcpp::XPtr<PopulationObjective> objective(pointer);
  return objective->value(encoded);
}

Rcpp::NumericVector population_objective_gradient_api(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  Rcpp::XPtr<PopulationObjective> objective(pointer);
  return objective->gradient(encoded);
}

Rcpp::NumericMatrix population_objective_hessian_api(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  Rcpp::XPtr<PopulationObjective> objective(pointer);
  return objective->hessian(encoded);
}

Rcpp::List population_objective_state_api(
    SEXP pointer, const Rcpp::NumericVector& encoded) {
  Rcpp::XPtr<PopulationObjective> objective(pointer);
  return objective->state(encoded);
}

Rcpp::List population_objective_telemetry_api(SEXP pointer) {
  Rcpp::XPtr<PopulationObjective> objective(pointer);
  return objective->telemetry();
}
