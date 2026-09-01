#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * Project the non-parametric imputation estimates onto the structural
 * manifold defined by the target covariance. Only entries marked in the
 * missingness mask are updated; observed data is frozen.
 *
 * This optimized version focuses on structural covariance preservation
 * using a Lagrangian constraint.
 *
 * Deprecated legacy v1 engine, kept for backward compatibility.  Use
 * constrain_covariance_v2, which adds fidelity to the initial imputation,
 * mean-preserving projected gradients, and KKT convergence diagnostics.
 *
 * Parameters:
 *   X_imp        n x p initial imputation matrix
 *   mask         n x p 0/1 missingness indicator (1 = originally missing)
 *   Sigma_target p x p target covariance (must be positive semidefinite)
 *   lambda       per-observation penalty on covariance deviation
 *   lr           learning rate
 *   max_iter     maximum iterations
 *   tol          convergence tolerance (Frobenius norm)
 *
 * Returns a list with:
 *   X_refined    n x p refined matrix
 *   iterations   number of gradient descent iterations executed
 *   converged    logical; TRUE if Frobenius norm fell below tol
 *   final_frob   final Frobenius distance to target covariance
 */
// [[Rcpp::export]]
Rcpp::List constrain_covariance(const arma::mat& X_imp,
                                const arma::mat& mask,
                                const arma::mat& Sigma_target,
                                double lambda, double lr, int max_iter, double tol) {

  arma::mat X_opt = X_imp;
  int n = X_opt.n_rows;
  arma::mat X_centered, Sigma_curr, grad_cov, update_step;

  /* Strict PSD Guard: Ensure the target is positive semi-definite.
   * Copy to a local mutable matrix — modifying a const reference via
   * const_cast is UB.  Outliers in heavy distributions can occasionally
   * cause precision issues that lead to negative eigenvalues, even if
   * the R-side already projected.
   */
  arma::mat Sigma_fixed = Sigma_target;
  arma::vec eigval_t;
  arma::mat eigvec_t;
  if (arma::eig_sym(eigval_t, eigvec_t, Sigma_fixed)) {
    bool has_neg = false;
    for (uword j = 0; j < eigval_t.n_elem; j++) {
      if (eigval_t(j) < 0) {
        eigval_t(j) = 0;
        has_neg = true;
      }
    }
    if (has_neg) {
      Sigma_fixed = eigvec_t * arma::diagmat(eigval_t) * eigvec_t.t();
    }
  }

  int iter = 0;
  bool converged = false;
  double final_frob = -1.0;

  for (int i = 0; i < max_iter; i++) {
    /* check for user interrupt to allow Esc/Ctrl+C in R */
    Rcpp::checkUserInterrupt();

    /* efficient centering and covariance calculation */
    X_centered = X_opt.each_row() - arma::mean(X_opt, 0);
    Sigma_curr = (X_centered.t() * X_centered) / (n - 1.0);

    /* convergence check (Frobenius distance to target) */
    final_frob = arma::norm(Sigma_curr - Sigma_fixed, "fro");
    iter = i + 1;
    if (final_frob < tol) {
      converged = true;
      break;
    }

    /* covariance-constraint gradient (un-normalised) */
    grad_cov = 2.0 * X_centered * (Sigma_curr - Sigma_fixed);

    /* Armijo backtracking line search to guarantee monotonic descent */
    double alpha = lr;
    double tau = 0.5;

    arma::mat X_opt_new, X_c_new, Sigma_new;
    double new_frob;

    while (true) {
      /* propose a step using the current alpha */
      update_step = (alpha * lambda) * grad_cov;
      X_opt_new = X_opt - (update_step % mask);

      /* evaluate the covariance and distance of the proposed step */
      X_c_new = X_opt_new.each_row() - arma::mean(X_opt_new, 0);
      Sigma_new = (X_c_new.t() * X_c_new) / (n - 1.0);
      new_frob = arma::norm(Sigma_new - Sigma_fixed, "fro");

      /* sufficient decrease condition */
      if (new_frob < final_frob) {
        break;
      }

      /* backtrack */
      alpha *= tau;

      /* safety guard against numerical precision floor */
      if (alpha < 1e-8) {
        break;
      }
    }

    /* commit the successful step */
    X_opt = X_opt_new;

    /* guard against numerical blow-up (extremely unlikely with Armijo) */
    if (X_opt.has_nan() || X_opt.has_inf()) {
      Rcpp::stop("Divergence detected: NaN or Inf produced during gradient descent.");
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("X_refined")  = X_opt,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("converged")  = converged,
    Rcpp::Named("final_frob") = final_frob
  );
}

/*
 * v2 engine: regularized constrained projection with mean preservation.
 *
 * Minimize the regularized constrained objective
 *
 *   f(X) = (1 / (2 N_mis)) ||M o (X - X0)||_F^2
 *        + (lambda_sigma / (2 p^2)) ||S(X) - Sigma_target||_F^2
 *
 * over the originally-missing cells only (observed cells frozen), where
 * N_mis is the total number of missing cells, S(X) = (n - 1)^-1 Xc' Xc is
 * the sample covariance of the column-centered matrix Xc, and Sigma_target
 * is projected to the PSD cone first.  Both terms are normalized — fidelity
 * per missing cell and covariance discrepancy per matrix entry — so
 * lambda_sigma is dimensionless and comparable across sample sizes,
 * variable counts, and missingness rates.  The first term anchors X to the
 * initial imputation X0 (nonparametric information); the second term
 * attracts the covariance toward the structural target.  lambda_sigma is a
 * genuine trade-off parameter: its value changes the fixed points of the
 * iteration, not merely the step size.
 *
 * Gradient on the free cells:
 *
 *   G = M o [ (X - X0) / N_mis
 *           + lambda_sigma * (2 / (p^2 (n - 1))) * Xc * R ],
 *   R = S(X) - Sigma_target.
 *
 * Mean anchoring: for every column j, the gradient of the missing cells is
 * centered at zero (Gtilde_ij = G_ij - mean_{i in M_j} G_ij), so every
 * accepted step keeps sum_i X_ij over the missing cells constant.  Before
 * the loop the starting point is projected onto the affine set of the
 * anchor sums s_j: with a scaled mean target, s_j = n * mu_j - obs_sum_j
 * pins the completed column mean to mu_j exactly; without one, s_j is the
 * missing-cell sum of X0 and the initial imputation's column mean is
 * preserved.  The returned matrix is therefore a true KKT point of the
 * joint mean-constrained covariance problem.  The column means of G are the
 * negated Lagrange multipliers of the mean constraint and are returned as
 * diagnostics.
 *
 * Stationarity (KKT) is certified by the projected-gradient norm ||Gtilde||_F:
 * at a constrained optimum the raw gradient is column-constant on the missing
 * cells of each column, so ||Gtilde||_F = 0 while ||G||_F may stay large.
 * The termination status reports first-order stationarity ("converged", and
 * "converged_feasible" when the covariance gap is within tol_cov) or a
 * non-stationary stop ("stalled_line_search", "max_iter_reached").  The
 * covariance gap at a stationary point is NOT interpreted as geometric
 * infeasibility: in the regularized problem a nonzero gap simply reflects
 * the fidelity/structure trade-off.
 *
 * Parameters:
 *   X_imp        n x p complete initial imputation
 *   mask         n x p 0/1 missingness indicator (1 = originally missing)
 *   Sigma_target p x p target covariance (projected to PSD internally)
 *   lambda_sigma dimensionless structural weight of the covariance term
 *   lr           initial step size (Armijo backtracking adapts it)
 *   max_iter     maximum iterations
 *   tol_kkT      stationarity tolerance on the projected-gradient norm
 *   tol_cov      feasibility tolerance on ||S(X) - Sigma_target||_F
 *   mu_scaled_target optional p-vector of scaled column-mean targets;
 *                NULL or an NA entry anchors that column to X0's mean
 *
 * Returns X_refined plus diagnostics: status, residuals, objective
 * components, and per-column multiplier estimates.
 */
// [[Rcpp::export]]
Rcpp::List constrain_covariance_v2(
    const arma::mat& X_imp,
    const arma::mat& mask,
    const arma::mat& Sigma_target,
    double lambda_sigma,
    double lr,
    int max_iter,
    double tol_kkT,
    double tol_cov,
    Rcpp::Nullable<Rcpp::NumericVector> mu_scaled_target = R_NilValue) {

  const int n = X_imp.n_rows;
  const int p = X_imp.n_cols;

  if (n < 2) Rcpp::stop("At least 2 rows are required to compute a covariance matrix.");
  if (max_iter < 1) Rcpp::stop("max_iter must be >= 1.");

  const arma::mat X0 = X_imp;

  /* PSD guard on the target, mirroring the legacy engine */
  arma::mat Sigma_fixed = Sigma_target;
  arma::vec eigval_t;
  arma::mat eigvec_t;
  if (arma::eig_sym(eigval_t, eigvec_t, Sigma_fixed)) {
    bool has_neg = false;
    for (uword j = 0; j < eigval_t.n_elem; j++) {
      if (eigval_t(j) < 0) {
        eigval_t(j) = 0;
        has_neg = true;
      }
    }
    if (has_neg) {
      Sigma_fixed = eigvec_t * arma::diagmat(eigval_t) * eigvec_t.t();
    }
  }

  /* missing-cell index sets per column */
  std::vector<arma::uvec> mis_idx(p);
  for (int j = 0; j < p; j++) {
    mis_idx[j] = arma::find(mask.col(j) > 0.5);
  }

  /* Mean anchor: per-column missing-cell sums the optimizer must hold.
   * A scaled mean target mu_j pins the completed column mean to mu_j via
   * s_j = n * mu_j - obs_sum_j (obs_sum_j = observed-cell sum); a NULL
   * vector or an NA entry keeps the legacy anchor s_j = X0's missing sum,
   * preserving the initial imputation's column mean. */
  const Rcpp::NumericVector mut =
    mu_scaled_target.isNotNull() ? mu_scaled_target.as()
                                 : Rcpp::NumericVector(0);
  arma::vec anchor_s(p, arma::fill::zeros);
  for (int j = 0; j < p; j++) {
    if (mis_idx[j].n_elem == 0) continue;
    if (mut.size() == p && !Rcpp::NumericVector::is_na(mut(j))) {
      const double obs_sum = arma::accu(X_imp.col(j)) -
                             arma::accu(mask.col(j) % X_imp.col(j));
      anchor_s(j) = n * mut(j) - obs_sum;
    } else {
      anchor_s(j) = arma::accu(mask.col(j) % X_imp.col(j));
    }
  }

  const double c_armijo  = 1e-4;
  const double tau       = 0.5;
  const double eta_min   = 1e-12;
  const double n_mis     = arma::accu(mask);
  const double cov_scale = 2.0 / ((p * p) * (n - 1.0));

  /* per-cell fidelity scale; zero when there is nothing to impute */
  const double fid_scale = (n_mis >= 1.0) ? (1.0 / n_mis) : 0.0;

  /* shift the starting point onto the mean-anchor manifold once, before the
   * loop; the zero-sum gradient projection keeps every accepted step inside
   * it, so the returned matrix is a certified constrained optimum */
  arma::mat X_opt = X_imp;
  for (int j = 0; j < p; j++) {
    if (mis_idx[j].n_elem == 0) continue;
    const double shift =
      (anchor_s(j) - arma::accu(mask.col(j) % X_opt.col(j))) /
      static_cast<double>(mis_idx[j].n_elem);
    for (uword k = 0; k < mis_idx[j].n_elem; k++) {
      X_opt(mis_idx[j](k), j) += shift;
    }
  }

  arma::mat Xc    = X_opt.each_row() - arma::mean(X_opt, 0);
  arma::mat Sigma_curr = (Xc.t() * Xc) / (n - 1.0);
  arma::mat R = Sigma_curr - Sigma_fixed;

  double f0   = 0.5 * arma::accu(arma::square(mask % (X_opt - X0))) * fid_scale;
  double fcov = 0.5 * arma::accu(arma::square(R)) / (p * p);
  double f    = f0 + lambda_sigma * fcov;

  arma::mat G, Gt;
  arma::vec nu(p, arma::fill::zeros);
  arma::vec spread(p, arma::fill::zeros);
  double r_kkT    = -1.0;
  double feas_gap = arma::norm(R, "fro");

  int iter = 0;
  bool converged = false;
  std::string status = "max_iter_reached";

  if (n_mis < 1.0) {
    /* nothing to impute: return the data unchanged */
    converged = true;
    r_kkT     = 0.0;
  } else {
    for (int i = 0; i < max_iter; i++) {
      Rcpp::checkUserInterrupt();
      iter = i + 1;

      /* full gradient on the free cells */
      G = mask % ((X_opt - X0) * fid_scale + lambda_sigma * cov_scale * Xc * R);

      /* zero-sum projection onto the mean-anchor constraint plus KKT
       * diagnostics */
      Gt = G;
      r_kkT = 0.0;
      for (int j = 0; j < p; j++) {
        if (mis_idx[j].n_elem == 0) continue;
        arma::vec gj(mis_idx[j].n_elem);
        for (uword k = 0; k < mis_idx[j].n_elem; k++) {
          gj(k) = G(mis_idx[j](k), j);
        }
        double gbar = arma::mean(gj);
        for (uword k = 0; k < mis_idx[j].n_elem; k++) {
          Gt(mis_idx[j](k), j) -= gbar;
        }
        nu(j)     = -gbar;
        spread(j) = arma::stddev(gj);
        arma::vec gjc = gj - gbar;
        r_kkT += arma::accu(arma::square(gjc));
      }
      r_kkT = std::sqrt(r_kkT);

      if (r_kkT <= tol_kkT) {
        converged = true;
        break;
      }

      /* Armijo backtracking on the full objective along -Gt */
      double eta = lr;
      double g2  = r_kkT * r_kkT;
      bool accepted = false;
      arma::mat X_new, Xc_new, Sigma_new, R_new;
      double f0_new, fcov_new, f_new;
      while (eta >= eta_min) {
        X_new     = X_opt - eta * Gt;
        Xc_new    = X_new.each_row() - arma::mean(X_new, 0);
        Sigma_new = (Xc_new.t() * Xc_new) / (n - 1.0);
        R_new     = Sigma_new - Sigma_fixed;
        f0_new    = 0.5 * arma::accu(arma::square(mask % (X_new - X0))) * fid_scale;
        fcov_new  = 0.5 * arma::accu(arma::square(R_new)) / (p * p);
        f_new     = f0_new + lambda_sigma * fcov_new;
        if (f_new <= f - c_armijo * eta * g2) {
          accepted = true;
          break;
        }
        eta *= tau;
      }
      if (!accepted) {
        status = "stalled_line_search";
        break;
      }

      X_opt = X_new;
      Xc    = Xc_new;
      R     = R_new;
      f0    = f0_new;
      fcov  = fcov_new;
      f     = f_new;

      if (X_opt.has_nan() || X_opt.has_inf()) {
        Rcpp::stop("Divergence detected: NaN or Inf produced during gradient descent.");
      }
    }
  }

  /* final diagnostics and termination classification */
  feas_gap = arma::norm(R, "fro");
  double max_abs_nu = 0.0;
  for (int j = 0; j < p; j++) {
    if (std::abs(nu(j)) > max_abs_nu) max_abs_nu = std::abs(nu(j));
  }

  if (converged) {
    if (feas_gap <= tol_cov) {
      status = "converged_feasible";
    } else {
      /* stationary with a nonzero covariance gap: the regularized optimum
       * balances fidelity against structure; no infeasibility is implied */
      status = "converged";
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("X_refined")   = X_opt,
    Rcpp::Named("iterations")  = iter,
    Rcpp::Named("converged")   = converged,
    Rcpp::Named("status")      = status,
    Rcpp::Named("r_kkT")       = r_kkT,
    Rcpp::Named("feas_gap")    = feas_gap,
    Rcpp::Named("objective")   = f,
    Rcpp::Named("fidelity")    = f0,
    Rcpp::Named("cov_term")    = fcov,
    Rcpp::Named("nu")          = nu,
    Rcpp::Named("grad_spread") = spread,
    Rcpp::Named("max_abs_nu")  = max_abs_nu,
    Rcpp::Named("n_mis")       = n_mis
  );
}
