#include "nm_pk_ad.h"
#include "nm_pk_pkadvan.h"
#include <R.h>
#include <cmath>
#include <vector>
#include <limits>

using namespace Rcpp;
using nm_pk::pkadvan::step_1_oral;
using nm_pk::pkadvan::step_1_iv_bolus;
using nm_pk::pkadvan::step_2_oral;

namespace {

typedef void (*add_grad_sidecar_fn)(SEXP, SEXP);

add_grad_sidecar_fn add_grad_sidecar_ptr() {
  static add_grad_sidecar_fn fn = NULL;
  if (fn == NULL) {
    fn = (add_grad_sidecar_fn) R_GetCCallable("LibeRtAD", "add_grad_sidecar");
    if (fn == NULL) {
      Rf_error("LibeRtAD add_grad_sidecar callable not found (is LibeRtAD loaded?)");
    }
  }
  return fn;
}

void ad_add_grad_sidecar(Environment node, SEXP increment) {
  add_grad_sidecar_ptr()(node, increment);
}

Environment pkg_env() {
  return Environment::namespace_env("LibeRtAD");
}

Environment newvar_sexp_pk(std::string myname, SEXP myval, std::string myop,
                           List myparents, List meta) {
  Environment Var = pkg_env()["Variable"];
  Function new_var = Var["new"];
  NumericVector z = as<NumericVector>(myval);
  NumericVector grad(z.size());
  std::fill(grad.begin(), grad.end(), 0.0);
  Environment res = new_var(
    Named("name") = myname,
    Named("value") = myval,
    Named("grad") = grad,
    Named("op") = myop,
    Named("parents") = myparents
  );
  res["meta"] = meta;
  return res;
}

double node_scalar(SEXP x) {
  if (Rf_inherits(x, "Variable") || Rf_inherits(x, "Constant")) {
    return as<NumericVector>(Environment(x)["value"])[0];
  }
  return as<NumericVector>(x)[0];
}

NumericVector align_grad_pk(NumericVector gv, R_xlen_t n) {
  if (gv.size() == n) {
    return gv;
  }
  if (gv.size() == 1) {
    NumericVector out(n);
    std::fill(out.begin(), out.end(), gv[0]);
    return out;
  }
  stop("Gradient size mismatch");
}

void oral2_step_vjp(
    double dt, double ka, double k20, double k23, double k32,
    double pg, double pa2, double pa3,
    double lg, double la2, double la3,
    double& ng, double& na2, double& na3,
    double& gka, double& gk20, double& gk23, double& gk32) {
  const double eps = 1e-7;
  double g0, a20, a30;
  step_2_oral(dt, ka, k20, k23, k32, pg, pa2, pa3, g0, a20, a30);

  auto contrib = [&](double g, double a2, double a3) {
    ng += lg * (g - g0) / eps;
    na2 += la2 * (a2 - a20) / eps;
    na3 += la3 * (a3 - a30) / eps;
  };

  {
    double g, a2, a3;
    step_2_oral(dt, ka, k20, k23, k32, pg + eps, pa2, pa3, g, a2, a3);
    contrib(g, a2, a3);
  }
  {
    double g, a2, a3;
    step_2_oral(dt, ka, k20, k23, k32, pg, pa2 + eps, pa3, g, a2, a3);
    contrib(g, a2, a3);
  }
  {
    double g, a2, a3;
    step_2_oral(dt, ka, k20, k23, k32, pg, pa2, pa3 + eps, g, a2, a3);
    contrib(g, a2, a3);
  }
  {
    double g, a2, a3;
    step_2_oral(dt, ka + eps, k20, k23, k32, pg, pa2, pa3, g, a2, a3);
    gka += lg * (g - g0) / eps + la2 * (a2 - a20) / eps + la3 * (a3 - a30) / eps;
  }
  {
    double g, a2, a3;
    step_2_oral(dt, ka, k20 + eps, k23, k32, pg, pa2, pa3, g, a2, a3);
    gk20 += lg * (g - g0) / eps + la2 * (a2 - a20) / eps + la3 * (a3 - a30) / eps;
  }
  {
    double g, a2, a3;
    step_2_oral(dt, ka, k20, k23 + eps, k32, pg, pa2, pa3, g, a2, a3);
    gk23 += lg * (g - g0) / eps + la2 * (a2 - a20) / eps + la3 * (a3 - a30) / eps;
  }
  {
    double g, a2, a3;
    step_2_oral(dt, ka, k20, k23, k32 + eps, pg, pa2, pa3, g, a2, a3);
    gk32 += lg * (g - g0) / eps + la2 * (a2 - a20) / eps + la3 * (a3 - a30) / eps;
  }
}

void oral1_step_vjp(
    double dt, double ka, double k10,
    double pg, double pc,
    double lg, double lc,
    double& ng, double& nc,
    double& gka, double& gk10) {
  const double eps = 1e-7;
  double g0, c0;
  step_1_oral(dt, ka, k10, pg, pc, g0, c0);

  auto contrib = [&](double g, double c) {
    ng += lg * (g - g0) / eps;
    nc += lc * (c - c0) / eps;
  };

  {
    double g, c;
    step_1_oral(dt, ka, k10, pg + eps, pc, g, c);
    contrib(g, c);
  }
  {
    double g, c;
    step_1_oral(dt, ka, k10, pg, pc + eps, g, c);
    contrib(g, c);
  }
  {
    double g, c;
    step_1_oral(dt, ka + eps, k10, pg, pc, g, c);
    gka += lg * (g - g0) / eps + lc * (c - c0) / eps;
  }
  {
    double g, c;
    step_1_oral(dt, ka, k10 + eps, pg, pc, g, c);
    gk10 += lg * (g - g0) / eps + lc * (c - c0) / eps;
  }
}

void bolus1_step_vjp(
    double dt, double k10, double pa, double la,
    double& na, double& gk10) {
  const double eps = 1e-7;
  double a0;
  step_1_iv_bolus(dt, k10, pa, a0);
  {
    double a;
    step_1_iv_bolus(dt, k10, pa + eps, a);
    na += la * (a - a0) / eps;
  }
  {
    double a;
    step_1_iv_bolus(dt, k10 + eps, pa, a);
    gk10 += la * (a - a0) / eps;
  }
}

}  // namespace

// [[Rcpp::export]]
Environment pk_oral2_trans4_block_var(
    NumericVector time,
    NumericVector amt,
    NumericVector f1,
    IntegerVector evid,
    SEXP ka_node,
    SEXP k20_node,
    SEXP k23_node,
    SEXP k32_node,
    SEXP vc_node,
    SEXP vp_node,
    int obs_cmp) {
  const int n = time.size();
  const double ka = node_scalar(ka_node);
  const double k20 = node_scalar(k20_node);
  const double k23 = node_scalar(k23_node);
  const double k32 = node_scalar(k32_node);
  const double vc = node_scalar(vc_node);
  const double vp = node_scalar(vp_node);

  std::vector<double> gut(n, 0.0), a2(n, 0.0), a3(n, 0.0);
  NumericVector conc(n);

  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = time[i] - time[i - 1];
      step_2_oral(dt, ka, k20, k23, k32,
                  gut[i - 1], a2[i - 1], a3[i - 1],
                  gut[i], a2[i], a3[i]);
    }
    if (evid[i] == 1) {
      gut[i] += amt[i] * f1[i];
    }
    if (obs_cmp == 1) {
      conc[i] = gut[i] / vc;
    } else if (obs_cmp == 3) {
      conc[i] = a3[i] / vp;
    } else {
      conc[i] = a2[i] / vc;
    }
  }

  List meta = List::create(
    Named("time") = time,
    Named("amt") = amt,
    Named("f1") = f1,
    Named("evid") = evid,
    Named("obs_cmp") = obs_cmp,
    Named("gut") = NumericVector(gut.begin(), gut.end()),
    Named("a2") = NumericVector(a2.begin(), a2.end()),
    Named("a3") = NumericVector(a3.begin(), a3.end()),
    Named("ka") = ka,
    Named("k20") = k20,
    Named("k23") = k23,
    Named("k32") = k32,
    Named("vc") = vc,
    Named("vp") = vp
  );

  return newvar_sexp_pk(
    "pk_oral2_block",
    conc,
    "pk_oral2_trans4",
    List::create(ka_node, k20_node, k23_node, k32_node, vc_node, vp_node),
    meta
  );
}

// [[Rcpp::export]]
Environment pk_oral1_block_var(
    NumericVector time,
    NumericVector amt,
    NumericVector f1,
    IntegerVector evid,
    SEXP ka_node,
    SEXP k10_node,
    SEXP v_node,
    int obs_cmp) {
  const int n = time.size();
  const double ka = node_scalar(ka_node);
  const double k10 = node_scalar(k10_node);
  const double v = node_scalar(v_node);

  std::vector<double> gut(n, 0.0), cen(n, 0.0);
  NumericVector conc(n);

  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = time[i] - time[i - 1];
      step_1_oral(dt, ka, k10, gut[i - 1], cen[i - 1], gut[i], cen[i]);
    }
    if (evid[i] == 1) {
      gut[i] += amt[i] * f1[i];
    }
    conc[i] = (obs_cmp == 1L) ? gut[i] / v : cen[i] / v;
  }

  List meta = List::create(
    Named("time") = time,
    Named("amt") = amt,
    Named("f1") = f1,
    Named("evid") = evid,
    Named("obs_cmp") = obs_cmp,
    Named("gut") = NumericVector(gut.begin(), gut.end()),
    Named("cen") = NumericVector(cen.begin(), cen.end()),
    Named("ka") = ka,
    Named("k10") = k10,
    Named("v") = v
  );

  return newvar_sexp_pk(
    "pk_oral1_block",
    conc,
    "pk_oral1",
    List::create(ka_node, k10_node, v_node),
    meta
  );
}

// [[Rcpp::export]]
Environment pk_bolus1_block_var(
    NumericVector time,
    NumericVector amt,
    IntegerVector evid,
    SEXP k10_node,
    SEXP v_node,
    int obs_cmp) {
  const int n = time.size();
  const double k10 = node_scalar(k10_node);
  const double v = node_scalar(v_node);

  std::vector<double> a1(n, 0.0);
  NumericVector conc(n);

  for (int i = 0; i < n; ++i) {
    if (evid[i] == 1) {
      a1[i] = amt[i];
    } else if (i > 0) {
      const double dt = time[i] - time[i - 1];
      step_1_iv_bolus(dt, k10, a1[i - 1], a1[i]);
    }
    conc[i] = a1[i] / v;
  }

  List meta = List::create(
    Named("time") = time,
    Named("amt") = amt,
    Named("evid") = evid,
    Named("obs_cmp") = obs_cmp,
    Named("a1") = NumericVector(a1.begin(), a1.end()),
    Named("k10") = k10,
    Named("v") = v
  );

  return newvar_sexp_pk(
    "pk_bolus1_block",
    conc,
    "pk_bolus1",
    List::create(k10_node, v_node),
    meta
  );
}

namespace nm_pk_ad {

void reverse_pk_oral2_trans4(Environment myvar, SEXP mygrad) {
  List meta = myvar["meta"];
  NumericVector time = meta["time"];
  IntegerVector evid = meta["evid"];
  int obs_cmp = meta["obs_cmp"];
  NumericVector gut = meta["gut"];
  NumericVector a2 = meta["a2"];
  NumericVector a3 = meta["a3"];
  const double ka = meta["ka"];
  const double k20 = meta["k20"];
  const double k23 = meta["k23"];
  const double k32 = meta["k32"];
  const double vc = meta["vc"];
  const double vp = meta["vp"];

  NumericVector upstream = align_grad_pk(as<NumericVector>(mygrad), gut.size());
  const int n = upstream.size();

  std::vector<double> lg(n, 0.0), la2(n, 0.0), la3(n, 0.0);
  for (int i = 0; i < n; ++i) {
    if (obs_cmp == 1) {
      lg[i] += upstream[i] / vc;
    } else if (obs_cmp == 3) {
      la3[i] += upstream[i] / vp;
    } else {
      la2[i] += upstream[i] / vc;
    }
  }

  double gka = 0.0, gk20 = 0.0, gk23 = 0.0, gk32 = 0.0, gvc = 0.0, gvp = 0.0;

  for (int i = 0; i < n; ++i) {
    if (obs_cmp == 1) {
      gvc += -upstream[i] * gut[i] / (vc * vc);
    } else if (obs_cmp == 3) {
      gvp += -upstream[i] * a3[i] / (vp * vp);
    } else {
      gvc += -upstream[i] * a2[i] / (vc * vc);
    }
  }

  for (int i = n - 1; i >= 0; --i) {
    if (i > 0) {
      const double dt = time[i] - time[i - 1];
      double ng = 0.0, na2 = 0.0, na3 = 0.0;
      oral2_step_vjp(
        dt, ka, k20, k23, k32,
        gut[i - 1], a2[i - 1], a3[i - 1],
        lg[i], la2[i], la3[i],
        ng, na2, na3,
        gka, gk20, gk23, gk32
      );
      lg[i - 1] += ng;
      la2[i - 1] += na2;
      la3[i - 1] += na3;
    }
  }

  List parents = myvar["parents"];
  ad_add_grad_sidecar(parents[0], NumericVector::create(gka));
  ad_add_grad_sidecar(parents[1], NumericVector::create(gk20));
  ad_add_grad_sidecar(parents[2], NumericVector::create(gk23));
  ad_add_grad_sidecar(parents[3], NumericVector::create(gk32));
  ad_add_grad_sidecar(parents[4], NumericVector::create(gvc));
  ad_add_grad_sidecar(parents[5], NumericVector::create(gvp));
}

void reverse_pk_oral1(Environment myvar, SEXP mygrad) {
  List meta = myvar["meta"];
  NumericVector time = meta["time"];
  int obs_cmp = meta["obs_cmp"];
  NumericVector gut = meta["gut"];
  NumericVector cen = meta["cen"];
  const double ka = meta["ka"];
  const double k10 = meta["k10"];
  const double v = meta["v"];

  NumericVector upstream = align_grad_pk(as<NumericVector>(mygrad), gut.size());
  const int n = upstream.size();

  std::vector<double> lg(n, 0.0), lc(n, 0.0);
  for (int i = 0; i < n; ++i) {
    if (obs_cmp == 1L) {
      lg[i] += upstream[i] / v;
    } else {
      lc[i] += upstream[i] / v;
    }
  }

  double gka = 0.0, gk10 = 0.0, gv = 0.0;
  for (int i = 0; i < n; ++i) {
    if (obs_cmp == 1L) {
      gv += -upstream[i] * gut[i] / (v * v);
    } else {
      gv += -upstream[i] * cen[i] / (v * v);
    }
  }

  for (int i = n - 1; i >= 0; --i) {
    if (i > 0) {
      const double dt = time[i] - time[i - 1];
      double ng = 0.0, nc = 0.0;
      oral1_step_vjp(
        dt, ka, k10,
        gut[i - 1], cen[i - 1],
        lg[i], lc[i],
        ng, nc,
        gka, gk10
      );
      lg[i - 1] += ng;
      lc[i - 1] += nc;
    }
  }

  List parents = myvar["parents"];
  ad_add_grad_sidecar(parents[0], NumericVector::create(gka));
  ad_add_grad_sidecar(parents[1], NumericVector::create(gk10));
  ad_add_grad_sidecar(parents[2], NumericVector::create(gv));
}

void reverse_pk_bolus1(Environment myvar, SEXP mygrad) {
  List meta = myvar["meta"];
  NumericVector time = meta["time"];
  NumericVector a1 = meta["a1"];
  const double k10 = meta["k10"];
  const double v = meta["v"];

  NumericVector upstream = align_grad_pk(as<NumericVector>(mygrad), a1.size());
  const int n = upstream.size();

  std::vector<double> la(n, 0.0);
  for (int i = 0; i < n; ++i) {
    la[i] += upstream[i] / v;
  }

  double gk10 = 0.0, gv = 0.0;
  for (int i = 0; i < n; ++i) {
    gv += -upstream[i] * a1[i] / (v * v);
  }

  for (int i = n - 1; i >= 0; --i) {
    if (i > 0) {
      const double dt = time[i] - time[i - 1];
      double na = 0.0;
      bolus1_step_vjp(dt, k10, a1[i - 1], la[i], na, gk10);
      la[i - 1] += na;
    }
  }

  List parents = myvar["parents"];
  ad_add_grad_sidecar(parents[0], NumericVector::create(gk10));
  ad_add_grad_sidecar(parents[1], NumericVector::create(gv));
}

void replay_pk_oral2_trans4(Environment node) {
  List meta = node["meta"];
  List parents = node["parents"];
  Environment tmp = pk_oral2_trans4_block_var(
    meta["time"], meta["amt"], meta["f1"], meta["evid"],
    parents[0], parents[1], parents[2], parents[3], parents[4], parents[5],
    as<int>(meta["obs_cmp"])
  );
  node.assign("value", tmp["value"]);
  node.assign("meta", tmp["meta"]);
}

void replay_pk_oral1(Environment node) {
  List meta = node["meta"];
  List parents = node["parents"];
  Environment tmp = pk_oral1_block_var(
    meta["time"], meta["amt"], meta["f1"], meta["evid"],
    parents[0], parents[1], parents[2],
    as<int>(meta["obs_cmp"])
  );
  node.assign("value", tmp["value"]);
  node.assign("meta", tmp["meta"]);
}

void replay_pk_bolus1(Environment node) {
  List meta = node["meta"];
  List parents = node["parents"];
  Environment tmp = pk_bolus1_block_var(
    meta["time"], meta["amt"], meta["evid"],
    parents[0], parents[1],
    as<int>(meta["obs_cmp"])
  );
  node.assign("value", tmp["value"]);
  node.assign("meta", tmp["meta"]);
}

}  // namespace nm_pk_ad
