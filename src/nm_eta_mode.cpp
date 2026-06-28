#include "nm_eta_mode.h"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using nm_subject::subject_nll;

namespace nm_eta {

SubjectData list_to_subject(const List& subj) {
  SubjectData s;
  s.time = subj["time"];
  s.amt = subj["amt"];
  s.rate = subj["rate"];
  s.f1 = subj["f1"];
  s.cmt = subj["cmt"];
  s.evid = subj["evid"];
  s.ss = subj["ss"];
  s.ii = subj["ii"];
  s.dv = subj["dv"];
  s.obs_idx = subj["obs_idx"];
  s.extras = nm_subject::extras_from_list(subj);
  return s;
}

namespace {

double eval_subject_nll(
    const SubjectData& subj,
    const NumericVector& eta,
    const NumericVector& theta,
    const NumericVector& omega,
    const NumericVector& sigma,
    const CharacterVector& pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss) {
  return subject_nll(
      subj.time, subj.amt, subj.rate, subj.f1, subj.cmt, subj.evid, subj.ss, subj.ii,
      subj.dv, subj.obs_idx, eta, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, true,
      List(), IntegerVector(), CharacterVector(), subj.extras);
}

void sort_simplex(
    std::vector<NumericVector>& verts,
    std::vector<double>& fvals) {
  const int m = static_cast<int>(verts.size());
  std::vector<int> ord(m);
  for (int i = 0; i < m; ++i) {
    ord[i] = i;
  }
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
    return fvals[a] < fvals[b];
  });
  std::vector<NumericVector> v2(m);
  std::vector<double> f2(m);
  for (int i = 0; i < m; ++i) {
    v2[i] = verts[ord[i]];
    f2[i] = fvals[ord[i]];
  }
  verts.swap(v2);
  fvals.swap(f2);
}

NumericVector centroid(const std::vector<NumericVector>& verts, int n) {
  const int m = static_cast<int>(verts.size());
  NumericVector c(n);
  for (int i = 0; i < m - 1; ++i) {
    for (int j = 0; j < n; ++j) {
      c[j] += verts[i][j];
    }
  }
  const double inv = 1.0 / static_cast<double>(m - 1);
  for (int j = 0; j < n; ++j) {
    c[j] *= inv;
  }
  return c;
}

}  // namespace

NumericVector find_eta_mode(
    const SubjectData& subj,
    const NumericVector& theta,
    const NumericVector& omega,
    const NumericVector& sigma,
    const CharacterVector& pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    SEXP eta_init,
    int max_iter) {
  const int n_eta = omega.size();
  NumericVector eta(n_eta);
  if (eta_init != R_NilValue) {
    eta = as<NumericVector>(eta_init);
    if (eta.size() != n_eta) {
      stop("eta_init length mismatch");
    }
  } else {
    std::fill(eta.begin(), eta.end(), 0.0);
  }
  if (n_eta == 0) {
    return eta;
  }

  auto nll_fn = [&](const NumericVector& e) {
    return eval_subject_nll(
        subj, e, theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss);
  };

  // Nelder-Mead (derivative-free; robust when finite-difference gradients mislead).
  const int m = n_eta + 1;
  std::vector<NumericVector> verts(m);
  std::vector<double> fvals(m);
  verts[0] = clone(eta);
  for (int i = 0; i < n_eta; ++i) {
    verts[i + 1] = clone(eta);
    const double step = (std::abs(eta[i]) > 1e-8) ? 0.05 * std::abs(eta[i]) : 0.05;
    verts[i + 1][i] += step;
  }
  for (int i = 0; i < m; ++i) {
    fvals[i] = nll_fn(verts[i]);
  }

  const double alpha = 1.0;
  const double beta = 2.0;
  const double gamma = 0.5;
  const double shrink = 0.5;
  const double ftol = 1e-8;
  const double xtol = 1e-6;

  for (int iter = 0; iter < max_iter; ++iter) {
    sort_simplex(verts, fvals);
    double fspread = fvals[m - 1] - fvals[0];
    double xspread = 0.0;
    for (int j = 0; j < n_eta; ++j) {
      double d = std::abs(verts[m - 1][j] - verts[0][j]);
      if (d > xspread) {
        xspread = d;
      }
    }
    if (fspread < ftol && xspread < xtol) {
      break;
    }

    NumericVector xo = centroid(verts, n_eta);
    NumericVector xr(n_eta);
    for (int j = 0; j < n_eta; ++j) {
      xr[j] = xo[j] + alpha * (xo[j] - verts[m - 1][j]);
    }
    const double fr = nll_fn(xr);

    if (fr >= fvals[0] && fr < fvals[m - 2]) {
      verts[m - 1] = xr;
      fvals[m - 1] = fr;
      continue;
    }

    if (fr < fvals[0]) {
      NumericVector xe(n_eta);
      for (int j = 0; j < n_eta; ++j) {
        xe[j] = xo[j] + beta * (xr[j] - xo[j]);
      }
      const double fe = nll_fn(xe);
      if (fe < fr) {
        verts[m - 1] = xe;
        fvals[m - 1] = fe;
      } else {
        verts[m - 1] = xr;
        fvals[m - 1] = fr;
      }
      continue;
    }

    NumericVector xc(n_eta);
    if (fr < fvals[m - 1]) {
      for (int j = 0; j < n_eta; ++j) {
        xc[j] = xo[j] + gamma * (xr[j] - xo[j]);
      }
    } else {
      for (int j = 0; j < n_eta; ++j) {
        xc[j] = xo[j] - gamma * (xo[j] - verts[m - 1][j]);
      }
    }
    const double fc = nll_fn(xc);
    if (fc < fvals[m - 1]) {
      verts[m - 1] = xc;
      fvals[m - 1] = fc;
      continue;
    }

    for (int i = 1; i < m; ++i) {
      for (int j = 0; j < n_eta; ++j) {
        verts[i][j] = verts[0][j] + shrink * (verts[i][j] - verts[0][j]);
      }
      fvals[i] = nll_fn(verts[i]);
    }
  }

  sort_simplex(verts, fvals);
  return verts[0];
}

}  // namespace nm_eta

// [[Rcpp::export]]
NumericVector nm_fit_eta_mode_cpp(
    List subject,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    SEXP eta_init = R_NilValue,
    int max_iter = 200L) {
  nm_eta::SubjectData subj = nm_eta::list_to_subject(subject);
  return nm_eta::find_eta_mode(
      subj, theta, omega, sigma, pred_lines,
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
      eta_init, max_iter);
}

// [[Rcpp::export]]
NumericMatrix nm_fit_all_eta_cpp(
    List subjects,
    NumericVector theta,
    NumericVector omega,
    NumericVector sigma,
    CharacterVector pred_lines,
    int advan,
    int trans,
    int obs_cmp,
    int dose_cmp,
    int n_transit,
    bool use_ode,
    int model_ss,
    SEXP eta_init = R_NilValue,
    int max_iter = 200L) {
  const int n_sub = subjects.size();
  const int n_eta = omega.size();
  NumericMatrix eta_out(n_sub, n_eta);
  NumericMatrix eta_warm;
  bool has_warm = false;
  if (eta_init != R_NilValue) {
    eta_warm = as<NumericMatrix>(eta_init);
    if (eta_warm.nrow() != n_sub || eta_warm.ncol() != n_eta) {
      stop("eta_init dimension mismatch");
    }
    has_warm = true;
  }
  for (int s = 0; s < n_sub; ++s) {
    NumericVector init_eta(n_eta);
    if (has_warm) {
      for (int j = 0; j < n_eta; ++j) {
        init_eta[j] = eta_warm(s, j);
      }
    }
    NumericVector row = nm_eta::find_eta_mode(
        nm_eta::list_to_subject(subjects[s]),
        theta, omega, sigma, pred_lines,
        advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss,
        has_warm ? wrap(init_eta) : R_NilValue,
        max_iter);
    for (int j = 0; j < n_eta; ++j) {
      eta_out(s, j) = row[j];
    }
  }
  return eta_out;
}
