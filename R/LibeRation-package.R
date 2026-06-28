#' LibeRation: NONMEM-Style Nonlinear Mixed-Effects Modelling
#'
#' LibeRation fits population PK/PD models with FO, FOCE, FOCEI, SAEM,
#' Laplace, IMP, and full Bayesian MCMC. Models follow NONMEM conventions
#' (ADVAN/TRANS, THETA/OMEGA/SIGMA, $PK/$ERROR). C++ PK solvers and
#' LibeRtAD-powered gradients are used when available.
#'
#' @section Core workflow:
#' * Build or import a model: [nm_model()], [nm_read_nonmem()], [nm_ctl_template()]
#' * Load data: [nm_dataset()]
#' * Estimate: [nm_est()]
#' * Goodness-of-fit: [predict()], [residuals()], [nm_etab()], [nm_add_cwres()]
#' * Inference: [nm_cov_step()], [nm_fit_standard_errors()], [nm_bootstrap_se()]
#'
#' @section Project workspace (Pirana-style):
#' * [nm_workspace_init()], [nm_workspace_create_project()],
#'   [nm_workspace_new_version()], [nm_workspace_save_run()]
#' * GUI: [liberation_shiny()]
#'
#' @section Synthetic examples:
#' * [nm_synthetic_catalog()], [nm_synthetic_theo()], [nm_synthetic_dataset()]
#'
#' @section Background jobs:
#' * [nm_job_submit()], [nm_job_status()], [nm_job_result()] (requires **callr**)
#'
#' @section Package options:
#' * `LibeRation.workspace` — workspace root directory
#' * `LibeRation.job_dir` — async job storage directory
#' * `LibeRation.profile` — enable estimation profiling
#' * `LibeRtAD.*` — passed through to the AD backend (see **LibeRtAD**)
#'
#' @seealso Package **LibeRtAD** for the automatic differentiation engine.
"_PACKAGE"
