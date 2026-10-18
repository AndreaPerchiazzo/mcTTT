.post_tau_numeric <- function(a, b1, b2, TimeF, tau, s0, t_fin) {
  L <- t_fin - s0
  if (L <= 0) return(0)
  U <- min(L, 60 / a)                      
  f <- function(r) exp(-a * r) * ((TimeF - t_fin + r) / (TimeF - tau))^b2
  b1 * stats::integrate(f, 0, U, rel.tol = 1e-10, subdivisions = 500L)$value
}

integral1 <- function(a, b1, b2, TimeF, t_in, t_fin, tau) {
  if (tau >= TimeF) stop("integral1: needed tau < TimeF")
  nn    <- max(length(t_in), length(t_fin))
  t_in  <- rep_len(t_in, nn)
  t_fin <- rep_len(t_fin, nn)
  out   <- numeric(nn)
  ok    <- t_fin > t_in
  if (!any(ok)) return(out)
  
  
  s1  <- pmin(t_fin, tau)
  pre <- ok & (t_in < s1)
  out[pre] <- b1 / a * exp(-a * (t_fin[pre] - s1[pre])) *
    (-expm1(-a * (s1[pre] - t_in[pre])))
  
  post <- ok & (t_fin > tau)
  if (any(post)) {
    s0  <- pmax(t_in[post], tau)
    tf  <- t_fin[post]
    lam <- a * (TimeF - tau)
    sh  <- b2 + 1
    lo  <- (TimeF - tf) / (TimeF - tau)     
    hi  <- (TimeF - s0) / (TimeF - tau)     
    ld  <- numeric(length(tf))
    up  <- lam * lo > sh                     
    if (any(up)) {
      l1 <- pgamma(lo[up], shape = sh, rate = lam, lower.tail = FALSE, log.p = TRUE)
      l2 <- pgamma(hi[up], shape = sh, rate = lam, lower.tail = FALSE, log.p = TRUE)
      ld[up] <- l1 + log1p(-exp(l2 - l1))
    }
    if (any(!up)) {
      l2 <- pgamma(hi[!up], shape = sh, rate = lam, log.p = TRUE)
      l1 <- pgamma(lo[!up], shape = sh, rate = lam, log.p = TRUE)
      ld[!up] <- l2 + log1p(-exp(l1 - l2))
    }
    lv  <- log(b1) + log(TimeF - tau) + a * (TimeF - tf) +
      lgamma(sh) - sh * log(lam) + ld
    val <- exp(lv)
    val[hi <= lo] <- 0
    bad <- !is.finite(val) | val < 0
    if (any(bad)) {
      idx <- which(bad)
      val[idx] <- vapply(idx, function(k)
        .post_tau_numeric(a, b1, b2, TimeF, tau, s0[k], tf[k]), numeric(1))
    }
    out[post] <- out[post] + val
  }
  out
}

SRCDM_m_t_s <- function(a, theta1, TimeF, t_in, t_fin, X_in, tau) {
  exp(-a * (t_fin - t_in)) * X_in +
    a * integral1(a, theta1[1], theta1[2], TimeF, t_in, t_fin, tau)
}
SRCDM_v_t_s <- function(a, theta2, TimeF, t_in, t_fin, tau) {
  integral1(2 * a, theta2[1]^2, 2 * theta2[2], TimeF, t_in, t_fin, tau)
}


par_names_full <- function(K) {
  c(paste0("a", 1:K), paste0("b_minus", 1:K), paste0("b_plus", 1:K),
    paste0("theta_minus", 1:K), paste0("theta_plus", 1:K))
}

default_bounds <- function(K) {
  lo <- c(a = 1e-7, b_minus = 1e-4, b_plus = 0.1, theta_minus = 1e-4, theta_plus = 0.1)
  up <- c(a = 50,   b_minus = 1,    b_plus = 2,   theta_minus = 1,    theta_plus = 2.5)
  list(lower = rep(lo, each = K), upper = rep(up, each = K))
}

make_spec <- function(K, start, fix_b_plus = NA, fix_theta_plus = NA, fix = NULL,
                      lower = NULL, upper = NULL) {
  nm <- par_names_full(K)
  if (length(start) != 5 * K) stop("start needs length 5K = ", 5 * K)
  start <- setNames(as.numeric(start), nm)
  fixed <- setNames(rep(NA_real_, 5 * K), nm)
  if (is.null(fix_b_plus))     fix_b_plus     <- NA
  if (is.null(fix_theta_plus)) fix_theta_plus <- NA
  fixed[paste0("b_plus", 1:K)]     <- rep_len(fix_b_plus, K)
  fixed[paste0("theta_plus", 1:K)] <- rep_len(fix_theta_plus, K)
  if (!is.null(fix)) {
    if (is.null(names(fix)) || any(!names(fix) %in% nm))
      stop("fix: admissible names = ", paste(nm, collapse = ", "))
    fixed[names(fix)] <- fix
  }
  b <- default_bounds(K)
  lower <- setNames(if (is.null(lower)) b$lower else rep_len(lower, 5 * K), nm)
  upper <- setNames(if (is.null(upper)) b$upper else rep_len(upper, 5 * K), nm)
  free  <- is.na(fixed)
  st    <- pmin(pmax(start, lower), upper)
  list(K = K, names = nm, free = free, fixed = fixed,
       start_free = st[free], lower = lower[free], upper = upper[free],
       lower_full = lower, upper_full = upper)
}

full_from_free <- function(par_free, spec) {
  full <- spec$fixed
  full[spec$free] <- par_free
  full
}

unpack <- function(full, K) {
  list(a      = full[1:K],
       theta1 = matrix(full[K + 1:(2 * K)], K, 2),        
       theta2 = matrix(full[3 * K + 1:(2 * K)], K, 2))    
}


bridge_moments <- function(full, X, Ti, TimeF, tau) {
  K  <- length(TimeF)
  p  <- unpack(full, K)
  n  <- length(X)
  if (length(Ti) != n) stop("X and Ti must have the same length")
  if (any(diff(Ti) <= 0)) stop("Ti must be strictly increasing.")
  if (max(Ti) >= min(TimeF)) stop("All observations must precede min(TimeF)")
  tin <- Ti[-n]; tfin <- Ti[-1]; xin <- X[-n]; xfin <- X[-1]
  N  <- n - 1
  k  <- s2 <- logf <- matrix(NA_real_, N, K)
  for (j in 1:K) {
    a <- p$a[j]; TT <- TimeF[j]; tj <- tau[j]
    th1 <- p$theta1[j, ]; th2 <- p$theta2[j, ]
    m_t  <- SRCDM_m_t_s(a, th1, TT, tin, tfin, xin, tj)
    m_T  <- SRCDM_m_t_s(a, th1, TT, tin, TT,   xin, tj)
    v_t  <- SRCDM_v_t_s(a, th2, TT, tin, tfin, tj)
    v_Tt <- SRCDM_v_t_s(a, th2, TT, tfin, TT,  tj)       
    e2   <- exp(-2 * a * (TT - tfin))
    v_T  <- e2 * v_t + v_Tt                               
    k[, j]  <- m_t - exp(-a * (TT - tfin)) * v_t / v_T * m_T
    s2[, j] <- v_t * v_Tt / v_T                           
    logf[, j] <- dnorm(xfin, mean = k[, j], sd = sqrt(s2[, j]), log = TRUE)
  }
  list(k = k, s2 = s2, logf = logf)
}

logdens_matrix <- function(full, X, Ti, TimeF, tau) {
  bridge_moments(full, X, Ti, TimeF, tau)$logf
}


rdcm_loglik_full <- function(full, X, Ti, TimeF, tau) {
  sum(logdens_matrix(full, X, Ti, TimeF, tau))
}

rdcm_negloglik <- function(par_free, spec, X, Ti, TimeF, tau, verbose = FALSE) {
  full <- full_from_free(par_free, spec)
  ll <- tryCatch(rdcm_loglik_full(full, X, Ti, TimeF, tau), error = function(e) -Inf)
  if (verbose) cat("\n", signif(c(ll, par_free), 6))
  if (!is.finite(ll)) return(1e10)
  -ll
}


fit_rdcm <- function(X, Ti, TimeF, tau, start,
                     fix_b_plus = NA, fix_theta_plus = NA, fix = NULL,
                     lower = NULL, upper = NULL, control = list(), verbose = FALSE) {
  if (length(TimeF) != 1 || length(tau) != 1) stop("fit_rdcm: TimeF e tau devono essere scalari")
  spec <- make_spec(1, start, fix_b_plus, fix_theta_plus, fix, lower, upper)
  if (is.null(control$parscale)) control$parscale <- pmax(abs(spec$start_free), 1e-3)
  if (is.null(control$maxit))    control$maxit    <- 500
  opt <- optim(spec$start_free, rdcm_negloglik, method = "L-BFGS-B",
               lower = spec$lower, upper = spec$upper, control = control,
               spec = spec, X = X, Ti = Ti, TimeF = TimeF, tau = tau, verbose = verbose)
  full <- full_from_free(opt$par, spec)
  ll   <- -opt$value
  kpar <- sum(spec$free)
  n_obs <- length(X) - 1
  structure(list(model = "RDCM", K = 1, par = full, spec = spec, loglik = ll,
                 n_par = kpar, AIC = -2 * ll + 2 * kpar, BIC = -2 * ll + log(n_obs) * kpar,
                 convergence = opt$convergence, message = opt$message, opt = opt,
                 X = X, Ti = Ti, TimeF = TimeF, tau = tau, control = control),
            class = "rdcm_fit")
}

rdcm_residuals <- function(fit) {
  bm <- bridge_moments(fit$par, fit$X, fit$Ti, fit$TimeF, fit$tau)
  as.numeric((fit$X[-1] - bm$k[, 1]) / sqrt(bm$s2[, 1]))
}

simulate_rdcm <- function(full, Ti, x0, TimeF, tau) {
  n <- length(Ti); X <- numeric(n); X[1] <- x0
  for (i in 2:n) {
    bm <- bridge_moments(full, X[(i - 1):i], Ti[(i - 1):i], TimeF, tau)
    X[i] <- rnorm(1, bm$k[1, 1], sqrt(bm$s2[1, 1]))
  }
  X
}


logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

forward_backward_log <- function(logE, pi0, Q) {
  N <- nrow(logE); K <- ncol(logE)
  lQ <- log(Q)
  la <- lb <- matrix(-Inf, N, K)
  la[1, ] <- log(pi0) + logE[1, ]                           
  if (N > 1) for (i in 2:N) for (j in 1:K)
    la[i, j] <- logE[i, j] + logsumexp(la[i - 1, ] + lQ[, j])    
  ll <- logsumexp(la[N, ])
  lb[N, ] <- 0
  if (N > 1) for (i in (N - 1):1) for (h in 1:K)
    lb[i, h] <- logsumexp(lQ[h, ] + logE[i + 1, ] + lb[i + 1, ])  
  gamma <- exp(la + lb - ll)
  gamma <- gamma / rowSums(gamma)
  xi_sum <- matrix(0, K, K)                                  
  if (N > 1) for (i in 1:(N - 1)) {
    lx <- outer(la[i, ], logE[i + 1, ] + lb[i + 1, ], "+") + lQ - ll
    xi_sum <- xi_sum + exp(lx)
  }
  P_hat <- xi_sum / rowSums(xi_sum)
  list(loglik = ll, gamma = gamma, P_hat = P_hat, log_alpha = la, log_beta = lb)
}

srdcm_loglik_full <- function(full, pi0, Q, X, Ti, TimeF, tau) {
  forward_backward_log(logdens_matrix(full, X, Ti, TimeF, tau), pi0, Q)$loglik
}

.mstep_obj <- function(par_free, spec, gamma, X, Ti, TimeF, tau) {
  full <- full_from_free(par_free, spec)
  E <- tryCatch(logdens_matrix(full, X, Ti, TimeF, tau), error = function(e) NULL)
  if (is.null(E)) return(1e10)
  v <- -sum(gamma * E)
  if (!is.finite(v)) return(1e10)
  v
}


fit_srdcm <- function(X, Ti, TimeF, tau, start, pi0, Q,
                      fix_b_plus = NA, fix_theta_plus = NA, fix = NULL,
                      lower = NULL, upper = NULL,
                      max_iter = 200, tol = 1e-7, mstep_maxit = 50,
                      control = list(), verbose = TRUE) {
  K <- length(TimeF)
  if (length(tau) != K) stop("tau e TimeF must have length K")
  spec <- make_spec(K, start, fix_b_plus, fix_theta_plus, fix, lower, upper)
  if (is.null(control$parscale)) control$parscale <- pmax(abs(spec$start_free), 1e-3)
  control$maxit <- mstep_maxit
  par_free <- spec$start_free
  ll_trace <- numeric(0)
  converged <- FALSE
  for (it in 1:max_iter) {
    full <- full_from_free(par_free, spec)
    fb <- forward_backward_log(logdens_matrix(full, X, Ti, TimeF, tau), pi0, Q)
    ll_trace <- c(ll_trace, fb$loglik)
    if (verbose) cat(sprintf("EM iter %3d   loglik = %.6f\n", it, fb$loglik))
    if (it > 1 && abs(ll_trace[it] - ll_trace[it - 1]) < tol * (1 + abs(ll_trace[it]))) {
      converged <- TRUE; break
    }
    opt <- optim(par_free, .mstep_obj, method = "L-BFGS-B",
                 lower = spec$lower, upper = spec$upper, control = control,
                 spec = spec, gamma = fb$gamma, X = X, Ti = Ti, TimeF = TimeF, tau = tau)
    par_free <- opt$par
    pi0 <- fb$gamma[1, ]
    Q   <- fb$P_hat
  }
  full <- full_from_free(par_free, spec)
  fb <- forward_backward_log(logdens_matrix(full, X, Ti, TimeF, tau), pi0, Q)
  kpar <- sum(spec$free) + K * (K - 1) + (K - 1)
  n_obs <- length(X) - 1
  structure(list(model = "SRDCM", K = K, par = full, spec = spec, pi0 = pi0, Q = Q,
                 loglik = fb$loglik, n_par = kpar,
                 AIC = -2 * fb$loglik + 2 * kpar, BIC = -2 * fb$loglik + log(n_obs) * kpar,
                 gamma = fb$gamma, decoded = max.col(fb$gamma, ties.method = "first"),
                 ll_trace = ll_trace, converged = converged, iterations = length(ll_trace),
                 X = X, Ti = Ti, TimeF = TimeF, tau = tau, control = control,
                 mstep_maxit = mstep_maxit, tol = tol),
            class = "srdcm_fit")
}

srdcm_residuals <- function(fit) {
  bm <- bridge_moments(fit$par, fit$X, fit$Ti, fit$TimeF, fit$tau)
  j  <- fit$decoded
  idx <- cbind(seq_along(j), j)
  as.numeric((fit$X[-1] - bm$k[idx]) / sqrt(bm$s2[idx]))
}

simulate_srdcm <- function(full, pi0, Q, Ti, x0, TimeF, tau) {
  K <- length(TimeF); n <- length(Ti)
  S <- integer(n - 1)
  S[1] <- sample.int(K, 1, prob = pi0)
  if (n > 2) for (i in 2:(n - 1)) S[i] <- sample.int(K, 1, prob = Q[S[i - 1], ])
  X <- numeric(n); X[1] <- x0
  for (i in 2:n) {
    bm <- bridge_moments(full, X[(i - 1):i], Ti[(i - 1):i], TimeF, tau)
    j <- S[i - 1]
    X[i] <- rnorm(1, bm$k[1, j], sqrt(bm$s2[1, j]))
  }
  list(X = X, S = S)
}

bootstrap_rdcm <- function(fit, B = 1000, seed = 1, jitter = 0,
                           save_file = NULL, save_every = 10, verbose = TRUE) {
  set.seed(seed)
  sp <- fit$spec
  est <- matrix(NA_real_, B, sum(sp$free), dimnames = list(NULL, names(sp$start_free)))
  conv <- rep(NA_integer_, B)
  for (b in 1:B) {
    Xb <- simulate_rdcm(fit$par, fit$Ti, fit$X[1], fit$TimeF, fit$tau)
    st <- fit$par
    st[sp$free] <- st[sp$free] + rnorm(sum(sp$free), 0, jitter)
    fb <- fit_rdcm(Xb, fit$Ti, fit$TimeF, fit$tau, start = st,
                   fix = sp$fixed[!sp$free], lower = sp$lower_full, upper = sp$upper_full,
                   control = fit$control)
    est[b, ] <- fb$par[sp$free]; conv[b] <- fb$convergence
    if (verbose && b %% save_every == 0) cat("bootstrap", b, "/", B, "\n")
    if (!is.null(save_file) && b %% save_every == 0)
      saveRDS(list(est = est[1:b, , drop = FALSE], conv = conv[1:b]), save_file)
  }
  ok <- conv == 0
  list(estimates = est, convergence = conv, n_converged = sum(ok),
       sd = apply(est[ok, , drop = FALSE], 2, sd))
}

bootstrap_srdcm <- function(fit, B = 1000, seed = 1, jitter = 0,
                            save_file = NULL, save_every = 10, verbose = TRUE) {
  set.seed(seed)
  sp <- fit$spec; K <- fit$K
  qn <- as.vector(outer(1:K, 1:K, function(h, k) paste0("p", h, k)))
  est <- matrix(NA_real_, B, sum(sp$free) + K * K,
                dimnames = list(NULL, c(names(sp$start_free), qn)))
  conv <- rep(NA, B)
  for (b in 1:B) {
    sim <- simulate_srdcm(fit$par, fit$pi0, fit$Q, fit$Ti, fit$X[1], fit$TimeF, fit$tau)
    st <- fit$par
    st[sp$free] <- st[sp$free] + rnorm(sum(sp$free), 0, jitter)
    fb <- fit_srdcm(sim$X, fit$Ti, fit$TimeF, fit$tau, start = st, pi0 = fit$pi0, Q = fit$Q,
                    fix = sp$fixed[!sp$free], lower = sp$lower_full, upper = sp$upper_full,
                    mstep_maxit = fit$mstep_maxit, tol = fit$tol, verbose = FALSE)
    est[b, ] <- c(fb$par[sp$free], as.vector(fb$Q)); conv[b] <- fb$converged
    if (verbose && b %% save_every == 0) cat("bootstrap", b, "/", B, "\n")
    if (!is.null(save_file) && b %% save_every == 0)
      saveRDS(list(est = est[1:b, , drop = FALSE], conv = conv[1:b]), save_file)
  }
  ok <- conv %in% TRUE
  list(estimates = est, converged = conv, n_converged = sum(ok),
       sd = apply(est[ok, , drop = FALSE], 2, sd))
}

param_table <- function(fit, boot = NULL) {
  out <- data.frame(param = names(fit$par), estimate = unname(fit$par),
                    status = ifelse(fit$spec$free, "free", "fixed"), row.names = NULL)
  if (!is.null(boot)) {
    out$sd <- NA_real_
    m <- match(names(boot$sd), out$param)
    out$sd[m[!is.na(m)]] <- boot$sd[!is.na(m)]
  }
  out
}

print.rdcm_fit <- function(x, ...) {
  cat(x$model, " | loglik =", format(x$loglik, digits = 10),
      "| free params =", x$n_par, "| AIC =", format(x$AIC, digits = 8),
      "| convergence =", x$convergence, "\n")
  print(param_table(x), digits = 6)
  invisible(x)
}

print.srdcm_fit <- function(x, ...) {
  cat(x$model, "K =", x$K, "| loglik =", format(x$loglik, digits = 10),
      "| params =", x$n_par, "| AIC =", format(x$AIC, digits = 8),
      "| EM converged =", x$converged, "in", x$iterations, "iter\n")
  print(param_table(x), digits = 6)
  cat("Q =\n"); print(round(x$Q, 4))
  invisible(x)
}
