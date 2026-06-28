#include <Rcpp.h>
#include <cmath>
#include <cstdlib>
#include <string>
#include <cctype>
#include "nm_pk_core.h"
#include "nm_pred_expr.h"
using namespace Rcpp;

namespace {

bool is_ode_des_line(const std::string& lhs, const std::string& rhs) {
  std::string ul = nm_pred::upper_copy(lhs);
  if (ul.rfind("DADT", 0) == 0) return true;
  if (ul == "F") return true;
  if (rhs.find("A(") != std::string::npos) return true;
  return false;
}

nm_pk::PkParams list_to_pkparams(List pred, nm_pk::PkParams base) {
  nm_pk::PkParams p = base;
  CharacterVector nms = pred.names();
  for (int i = 0; i < pred.size(); ++i) {
    std::string nm = as<std::string>(nms[i]);
    double v = as<double>(pred[i]);
    for (char& c : nm) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (nm == "KA") p.ka = v;
    else if (nm == "CL") p.cl = v;
    else if (nm == "VC") p.vc = v;
    else if (nm == "VP") p.vp = v;
    else if (nm == "VP2") p.vp2 = v;
    else if (nm == "Q2" || nm == "Q") p.q2 = v;
    else if (nm == "Q3") p.q3 = v;
    else if (nm == "Q4") p.q4 = v;
    else if (nm == "K10" || nm == "K") p.k10 = v;
    else if (nm == "K12") p.k12 = v;
    else if (nm == "K21") p.k21 = v;
    else if (nm == "K13") p.k13 = v;
    else if (nm == "K31") p.k31 = v;
    else if (nm == "K14") p.k24 = v;
    else if (nm == "K41") p.k42 = v;
    else if (nm == "K20") p.k20 = v;
    else if (nm == "K23") p.k23 = v;
    else if (nm == "K32") p.k32 = v;
    else if (nm == "K24") p.k24 = v;
    else if (nm == "K42") p.k42 = v;
    else if (nm == "KTR") p.ktr = v;
    else if (nm == "KE") p.k10 = v;
    else if (nm == "V") {
      if (p.v1 <= 0.0) p.v1 = v;
      if (p.vc <= 0.0) p.vc = v;
    }
    else if (nm == "KMF") { p.kmf = v; p.has_metabolite = true; }
    else if (nm == "KME") { p.kme = v; p.has_metabolite = true; }
    else if (nm == "V1") p.v1 = v;
    else if (nm == "V2") p.v2 = v;
    else if (nm == "V3") p.v3 = v;
    else if (nm == "V4") p.v4 = v;
    else if (nm == "VSS") p.vss = v;
    else if (nm == "AOB") p.aob = v;
    else if (nm == "ALPHA") p.pk_alpha = v;
    else if (nm == "BETA") p.pk_beta = v;
    else if (nm == "GAMMA") p.pk_gamma = v;
    else if (nm.size() > 1 && nm[0] == 'S') {
      const int idx = std::atoi(nm.c_str() + 1);
      if (idx >= 1 && idx <= nm_pk::kMaxScale) {
        p.scale[static_cast<size_t>(idx - 1)] = v;
      }
    }
    else if (nm.size() > 1 && nm[0] == 'F') {
      const int idx = std::atoi(nm.c_str() + 1);
      if (idx >= 1 && idx <= nm_pk::kMaxScale) {
        p.f[static_cast<size_t>(idx - 1)] = v;
      }
    }
    else if (nm == "ALAG1" || nm == "ALAG") p.alag1 = v;
    else if (nm == "VMAX") p.vmax = v;
    else if (nm == "KM") p.km = v;
    else if (nm == "VM") {
      p.vm = v;
      if (p.v1 <= 0.0) p.v1 = v;
    }
    else {
      // allow sequential references to earlier PRED assignments
    }
  }
  return p;
}

}  // namespace

// [[Rcpp::export]]
bool nm_pred_expr_check_cpp(CharacterVector pred_lines) {
  Rcpp::NumericVector th(10), et(10);
  for (int i = 0; i < 10; ++i) { th[i] = 1.0; et[i] = 0.0; }
  nm_pred::ExprEnv env;
  env.theta = th;
  env.eta = et;
  for (int i = 0; i < pred_lines.size(); ++i) {
    std::string line = as<std::string>(pred_lines[i]);
    std::string rhs = nm_pred::pred_rhs(line);
    if (rhs.empty()) return false;
    try {
      double v = nm_pred::ExprParser::eval_with_env(rhs, env);
      env.vars[nm_pred::upper_copy(nm_pred::pred_lhs(line))] = v;
    } catch (...) {
      return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
List nm_eval_pred_cpp(CharacterVector pred_lines,
                      NumericVector theta,
                      NumericVector eta,
                      List covariates = List(),
                      CharacterVector des_lines = CharacterVector()) {
  List out(pred_lines.size());
  CharacterVector names(pred_lines.size());
  nm_pred::ExprEnv env;
  env.theta = theta;
  env.eta = eta;
  if (covariates.size() > 0) {
    CharacterVector cnames = covariates.names();
    for (int i = 0; i < covariates.size(); ++i) {
      std::string nm = as<std::string>(cnames[i]);
      env.vars[nm_pred::upper_copy(nm)] = as<double>(covariates[i]);
    }
  }
  for (int i = 0; i < pred_lines.size(); ++i) {
    std::string line = as<std::string>(pred_lines[i]);
    std::string lhs = nm_pred::pred_lhs(line);
    std::string rhs = nm_pred::pred_rhs(line);
    names[i] = lhs;
    double v = nm_pred::ExprParser::eval_with_env(rhs, env);
    env.vars[nm_pred::upper_copy(lhs)] = v;
    out[i] = v;
  }
  out.names() = names;
  if (des_lines.size() > 0) {
    for (int i = 0; i < des_lines.size(); ++i) {
      std::string line = as<std::string>(des_lines[i]);
      std::string lhs = nm_pred::pred_lhs(line);
      std::string rhs = nm_pred::pred_rhs(line);
      if (lhs.empty() || rhs.empty()) continue;
      if (is_ode_des_line(lhs, rhs)) continue;
      double v = nm_pred::ExprParser::eval_with_env(rhs, env);
      env.vars[nm_pred::upper_copy(lhs)] = v;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector nm_pk_route_r(
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss,
    NumericVector time, NumericVector amt, NumericVector rate,
    NumericVector f1, IntegerVector cmt, IntegerVector evid,
    IntegerVector ss, NumericVector ii,
    List pk_params,
    NumericVector s1 = NumericVector(),
    NumericVector s2 = NumericVector(),
    NumericVector s3 = NumericVector(),
    NumericVector s4 = NumericVector(),
    NumericMatrix scale_mat = NumericMatrix(),
    bool use_data_scale = false,
    NumericMatrix f_mat = NumericMatrix(),
    bool use_data_f = false) {
  nm_pk::SubjectEvents ev;
  ev.time = time;
  ev.amt = amt;
  ev.rate = rate;
  ev.f1 = f1;
  ev.cmt = cmt;
  ev.evid = evid;
  ev.ss = ss;
  ev.ii = ii;
  ev.s1 = s1;
  ev.s2 = s2;
  ev.s3 = s3;
  ev.s4 = s4;
  ev.scale_mat = scale_mat;
  ev.use_data_scale = use_data_scale;
  ev.f_mat = f_mat;
  ev.use_data_f = use_data_f;
  nm_pk::PkParams p = list_to_pkparams(pk_params, nm_pk::PkParams());
  return nm_pk_route_cpp(advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, ev, p);
}

// [[Rcpp::export]]
List nm_pk_route_detail_r(
    int advan, int trans, int obs_cmp, int dose_cmp,
    int n_transit, bool use_ode, int model_ss, int n_state,
    NumericVector time, NumericVector amt, NumericVector rate,
    NumericVector f1, IntegerVector cmt, IntegerVector evid,
    IntegerVector ss, NumericVector ii,
    List pk_params,
    NumericVector s1 = NumericVector(),
    NumericVector s2 = NumericVector(),
    NumericVector s3 = NumericVector(),
    NumericVector s4 = NumericVector(),
    NumericMatrix scale_mat = NumericMatrix(),
    bool use_data_scale = false,
    NumericMatrix f_mat = NumericMatrix(),
    bool use_data_f = false) {
  nm_pk::SubjectEvents ev;
  ev.time = time;
  ev.amt = amt;
  ev.rate = rate;
  ev.f1 = f1;
  ev.cmt = cmt;
  ev.evid = evid;
  ev.ss = ss;
  ev.ii = ii;
  ev.s1 = s1;
  ev.s2 = s2;
  ev.s3 = s3;
  ev.s4 = s4;
  ev.scale_mat = scale_mat;
  ev.use_data_scale = use_data_scale;
  ev.f_mat = f_mat;
  ev.use_data_f = use_data_f;
  nm_pk::PkParams p = list_to_pkparams(pk_params, nm_pk::PkParams());
  return nm_pk_route_detail_cpp(
      advan, trans, obs_cmp, dose_cmp, n_transit, use_ode, model_ss, ev, p, n_state);
}

// [[Rcpp::export]]
bool nm_cpp_advan_supported(int advan, int trans) {
  if (advan >= 1 && advan <= 4) return true;
  if (advan == 11 || advan == 12) return true;
  if (advan >= 6 && advan <= 13) return true;
  if (advan == 9 || advan == 10) return true;
  return false;
}
