# =============================================================================
# core_oneshot.R
# Shared core for the shelf-life paper (P96_FOOD).
# Implements the non-parametric and semi-parametric methods of
#   Zhu, X. and Balakrishnan, N. (2022). One-shot device test data analysis
#   using non-parametric and semi-parametric inferential methods and
#   applications. Reliability Engineering and System Safety 221, 108319.
# Base R only. Every other script sources this file.
#
# Notation (as in the paper)
#   K        number of test environments (storage temperatures), k = 1..K
#   IT_j     distinct inspection times, j = 1..J (union over environments)
#   N_kj     number of units tested in environment k at IT_j
#   r_kj     number of failures (unacceptable units) among them
#   S(IT_j)  survival probability; S_1 = baseline (reference environment k = 1)
#   a_j      S(IT_j) / S(IT_{j-1}),  S(IT_0) = 1
#   m_kj     latent number of failures in (IT_{j-1}, IT_j]
#   lambda_k proportional hazard rate, S_k(t) = S_1(t)^lambda_k, lambda_1 = 1
#   beta     log-link coefficient, lambda_k = exp{(x_k - x_1) beta}
# =============================================================================

# ---------------------------------------------------------------- settings ---
eps_loglik <- 1e-12   # absolute change in observed log-likelihood
eps_par    <- 1e-9    # largest relative change in the parameters
max_iter   <- 2e6     # EM iteration cap (EM is sublinear at active constraints)
z975       <- qnorm(0.975)
var_scale  <- 1e3
var_scale_tex <- "10^{3}"

# ------------------------------------------------------------ data helpers ---
# A one-shot data set is a data.frame with columns
#   env (1..K), x (stress covariate of env), IT, N, r
os_data <- function(env, x, IT, N, r) {
  d <- data.frame(env = as.integer(env), x = as.numeric(x),
                  IT = as.numeric(IT), N = as.numeric(N), r = as.numeric(r))
  d <- d[order(d$env, d$IT), ]
  rownames(d) <- NULL
  d
}

xlogy <- function(x, y) {
  y <- rep_len(y, length(x)); out <- numeric(length(x)); i <- x != 0
  out[i] <- x[i] * log(y[i]); out
}

# =============================================================================
# Section 2: one test environment (K = 1)
# =============================================================================

# Incomplete log-likelihood (1), without the binomial constant
loglik_np <- function(S, N, r) sum(xlogy(r, 1 - S) + xlogy(N - r, S))

# Lemma 1 / independent route: weighted isotonic regression (PAVA) of the
# failure proportions r_j / N_j with weights N_j gives the NPMLE of
# F(IT_j) = 1 - S(IT_j) for current-status data (Ayer et al., 1955).
pava_np <- function(N, r) {
  p <- r / N; w <- N
  blk_p <- p; blk_w <- w; blk_n <- rep(1, length(p))
  i <- 1
  while (i < length(blk_p)) {
    if (blk_p[i] > blk_p[i + 1] + 1e-15) {
      pw <- (blk_p[i] * blk_w[i] + blk_p[i + 1] * blk_w[i + 1]) /
            (blk_w[i] + blk_w[i + 1])
      blk_p[i] <- pw; blk_w[i] <- blk_w[i] + blk_w[i + 1]
      blk_n[i] <- blk_n[i] + blk_n[i + 1]
      blk_p <- blk_p[-(i + 1)]; blk_w <- blk_w[-(i + 1)]; blk_n <- blk_n[-(i + 1)]
      if (i > 1) i <- i - 1
    } else i <- i + 1
  }
  1 - rep(blk_p, blk_n)
}

# EM algorithm of Section 2.1, equations (3) and (4)
em_np <- function(IT, N, r, S0 = NULL) {
  J <- length(IT); Ntot <- sum(N)
  S <- if (is.null(S0)) 1 - (1:J) / (J + 1) else S0
  ll_old <- loglik_np(S, N, r)
  for (h in 1:max_iter) {
    Sext <- c(1, S, 0)                       # S(IT_0) ... S(IT_{J+1})
    dS   <- Sext[1:(J + 1)] - Sext[2:(J + 2)]  # Delta S(IT_i), i = 1..J+1
    m <- numeric(J + 1)
    for (j in 1:J) {
      if (S[j] > 0)  m[(j + 1):(J + 1)] <- m[(j + 1):(J + 1)] +
                       (N[j] - r[j]) * dS[(j + 1):(J + 1)] / S[j]
      if (S[j] < 1)  m[1:j] <- m[1:j] + r[j] * dS[1:j] / (1 - S[j])
    }
    S_new <- rev(cumsum(rev(m)))[2:(J + 1)] / Ntot     # equation (4)
    ll_new <- loglik_np(S_new, N, r)
    rel <- max(abs(S_new - S) / pmax(abs(S), 1e-300))
    S <- S_new
    if (abs(ll_new - ll_old) < eps_loglik && rel < eps_par) break
    ll_old <- ll_new
  }
  list(S = S, loglik = ll_new, iter = h, m = m)
}

# Section 2.2-2.3: variance (6) and the three intervals (7), (8), (9)
ci_S <- function(S, se, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  plain  <- cbind(S - z * se, S + z * se)
  U      <- z * se / (S * log(S))              # log-log, equation (8)
  a <- S^exp(-U); b <- S^exp(U)
  loglog <- cbind(pmin(a, b), pmax(a, b))
  V      <- exp(z * se / (S * (1 - S)))        # logit, equation (9)
  logit  <- cbind(S / (S + (1 - S) * V), S / (S + (1 - S) / V))
  list(plain = plain, loglog = loglog, logit = logit)
}
se_np <- function(S, N) sqrt(S * (1 - S) / N)

# =============================================================================
# Section 3: Nelson-Aalen based estimation
# =============================================================================
na_np <- function(IT, N, r, S0 = NULL) {
  J <- length(IT); Ntot <- sum(N)
  S <- if (is.null(S0)) 1 - (1:J) / (J + 1) else S0
  for (h in 1:max_iter) {
    Sext <- c(1, S, 0); dS <- Sext[1:(J + 1)] - Sext[2:(J + 2)]
    m <- numeric(J + 1)
    for (j in 1:J) {
      if (S[j] > 0) m[(j + 1):(J + 1)] <- m[(j + 1):(J + 1)] +
                      (N[j] - r[j]) * dS[(j + 1):(J + 1)] / S[j]
      if (S[j] < 1) m[1:j] <- m[1:j] + r[j] * dS[1:j] / (1 - S[j])
    }
    risk  <- Ntot - c(0, cumsum(m[1:(J - 1)]))
    S_new <- exp(-cumsum(m[1:J] / risk))
    rel <- max(abs(S_new - S) / pmax(abs(S), 1e-300))
    S <- S_new
    if (rel < eps_par) break
  }
  list(S = S, iter = h)
}

# Non-parametric bootstrap: units are resampled within each cell
boot_cells <- function(N, r) rbinom(length(N), N, r / N)

# =============================================================================
# Section 4: proportional hazards model, K environments
# =============================================================================
# Build the internal structure: union grid of inspection times, the grid
# index of every cell, and the covariate contrast delta_k = x_k - x_1.
ph_setup <- function(d) {
  grid <- sort(unique(d$IT))
  K <- max(d$env)
  xk <- tapply(d$x, d$env, function(v) v[1])
  list(d = d, grid = grid, J = length(grid), K = K,
       g = match(d$IT, grid), xk = as.numeric(xk),
       delta = as.numeric(xk - xk[1]))
}

lambda_of <- function(st, lam_free = NULL, beta = NULL) {
  if (!is.null(beta)) exp(st$delta * beta) else c(1, lam_free)
}

# Incomplete log-likelihood (10)-(11)
loglik_ph <- function(st, S1, lam) {
  d <- st$d
  Sk <- S1[st$g]^lam[d$env]
  sum(xlogy(d$r, 1 - Sk) + xlogy(d$N - d$r, Sk))
}

# E-step: expected failures m_kj in each interval of the union grid,
# equation (3) applied environment by environment
estep_ph <- function(st, S1, lam) {
  J <- st$J; d <- st$d
  m <- matrix(0, st$K, J + 1)
  for (k in 1:st$K) {
    Sk <- S1^lam[k]; Sext <- c(1, Sk, 0)
    dS <- Sext[1:(J + 1)] - Sext[2:(J + 2)]
    cells <- which(d$env == k)
    for (c in cells) {
      j <- st$g[c]
      if (Sk[j] > 0) m[k, (j + 1):(J + 1)] <- m[k, (j + 1):(J + 1)] +
                       (d$N[c] - d$r[c]) * dS[(j + 1):(J + 1)] / Sk[j]
      if (Sk[j] < 1) m[k, 1:j] <- m[k, 1:j] + d$r[c] * dS[1:j] / (1 - Sk[j])
    }
  }
  m
}

# M-step pieces. With u_j = -log a_j >= 0, the complete log-likelihood (12) is
#   sum_k sum_j [ m_kj log(1 - exp(-lambda_k u_j)) - lambda_k u_j R_kj ],
#   R_kj = sum_{i > j} m_ki.
# Score in u_j : sum_k m_kj lambda_k / (exp(lambda_k u_j) - 1) - sum_k lambda_k R_kj
solve_u <- function(mj, Rj, lam) {
  A <- sum(lam * Rj)
  if (sum(mj) <= 0) return(0)
  if (A <= 0) return(Inf)
  f <- function(lu) { u <- exp(lu); sum(mj * lam / expm1(lam * u)) - A }
  lo <- -30; hi <- 10
  while (f(hi) > 0) hi <- hi + 10
  while (f(lo) < 0 && lo > -700) lo <- lo - 30
  if (f(lo) < 0) return(0)
  exp(uniroot(f, c(lo, hi), tol = 1e-15)$root)
}
# Score in lambda_k : sum_j m_kj u_j / (exp(lambda_k u_j) - 1) - sum_j u_j R_kj
score_lam_k <- function(lam, mk, Rk, u) {
  ok <- u > 0 & is.finite(u)
  sum(mk[ok] * u[ok] / expm1(lam * u[ok])) - sum(u[ok] * Rk[ok])
}
solve_lam <- function(mk, Rk, u) {
  f <- function(ll) score_lam_k(exp(ll), mk, Rk, u)
  lo <- -20; hi <- 20
  while (f(lo) < 0 && lo > -200) lo <- lo - 20
  while (f(hi) > 0 && hi < 200) hi <- hi + 20
  exp(uniroot(f, c(lo, hi), tol = 1e-15)$root)
}
solve_beta <- function(m, R, u, delta, beta0) {
  f <- function(b) {
    lam <- exp(delta * b)
    sum(sapply(seq_along(delta), function(k)
      delta[k] * lam[k] * score_lam_k(lam[k], m[k, ], R[k, ], u)))
  }
  w <- 1 / max(abs(delta)); lo <- beta0 - w; hi <- beta0 + w
  while (f(lo) < 0) lo <- lo - w
  while (f(hi) > 0) hi <- hi + w
  uniroot(f, c(lo, hi), tol = 1e-15)$root
}

# EM (ECM: one conditional maximisation cycle per iteration) for the PH model.
# link = "none": lambda_2..lambda_K free;  link = "log": lambda_k = exp(delta_k beta)
em_ph <- function(d, link = c("none", "log"), S1_0 = NULL, lam_0 = NULL,
                  beta_0 = 0, it_max = max_iter) {
  link <- match.arg(link)
  st <- ph_setup(d); J <- st$J; K <- st$K
  S1 <- if (is.null(S1_0)) 1 - (1:J) / (J + 1) else S1_0
  lam <- if (link == "log") exp(st$delta * beta_0) else
         if (is.null(lam_0)) rep(1, K) else c(1, lam_0)
  beta <- beta_0
  ll_old <- loglik_ph(st, S1, lam)
  for (h in 1:it_max) {
    m <- estep_ph(st, S1, lam)
    R <- t(apply(m, 1, function(v) rev(cumsum(rev(v)))[-1]))   # K x J
    if (K == 1) R <- matrix(R, nrow = 1)
    # (i) a_j given lambda, equation (16)
    u <- sapply(1:J, function(j) solve_u(m[, j], R[, j], lam))
    # (ii) lambda given a, equation (13), or beta, equation (15)
    if (link == "none") {
      lam_new <- if (K == 1) 1 else
        c(1, sapply(2:K, function(k) solve_lam(m[k, 1:J], R[k, ], u)))
    } else {
      beta <- solve_beta(m[, 1:J, drop = FALSE], R, u, st$delta, beta)
      lam_new <- exp(st$delta * beta)
    }
    S1_new <- exp(-cumsum(u))
    ll_new <- loglik_ph(st, S1_new, lam_new)
    par_old <- c(S1, lam); par_new <- c(S1_new, lam_new)
    rel <- max(abs(par_new - par_old) / pmax(abs(par_old), 1e-300))
    S1 <- S1_new; lam <- lam_new
    if (abs(ll_new - ll_old) < eps_loglik && rel < eps_par) break
    ll_old <- ll_new
  }
  out <- list(st = st, link = link, S1 = S1, lam = lam, loglik = ll_new,
              iter = h, converged = h < it_max)
  if (link == "log") out$beta <- beta
  out
}

# EM followed, when EM has not met the stopping rule within it_em iterations
# (sublinear convergence at active monotonicity constraints), by direct
# maximisation started at the EM point. The two routes are compared and the
# larger log-likelihood is kept; the relative difference is stored.
fit_ph <- function(d, link = c("none", "log"), it_em = 20000, ...) {
  link <- match.arg(link)
  f <- em_ph(d, link, it_max = it_em, ...)
  if (f$converged) { f$route <- "EM"; return(f) }
  u <- pmax(-diff(c(0, log(f$S1))), 1e-12)
  start <- c(log(u), if (link == "log") f$beta else log(f$lam[-1]))
  g <- direct_ph(d, link, start = start)
  f$route <- "EM + direct"
  f$rel_diff_em_direct <- max(abs(c(g$S1, g$lam) - c(f$S1, f$lam)) / abs(c(g$S1, g$lam)))
  if (g$loglik >= f$loglik) {
    f$S1 <- g$S1; f$lam <- g$lam; f$loglik <- g$loglik
    if (link == "log") f$beta <- g$beta
  }
  f
}

# Direct maximisation of (11), an independent route to the same MLE.
# S1_j = exp(-sum_{i<=j} exp(phi_i)); log lambda_k or beta unconstrained.
direct_ph <- function(d, link = c("none", "log"), start = NULL) {
  link <- match.arg(link)
  st <- ph_setup(d); J <- st$J; K <- st$K
  unpack <- function(th) {
    S1 <- exp(-cumsum(exp(th[1:J])))
    lam <- if (link == "log") exp(st$delta * th[J + 1]) else c(1, exp(th[J + seq_len(K - 1)]))
    list(S1 = S1, lam = lam)
  }
  nll <- function(th) { p <- unpack(th); -loglik_ph(st, p$S1, p$lam) }
  if (is.null(start)) start <- c(log(rep(1 / J, J)),
                                 if (link == "log") 0 else rep(0, K - 1))
  fit <- optim(start, nll, method = "BFGS",
               control = list(reltol = 1e-16, maxit = 1e5))
  fit <- optim(fit$par, nll, method = "Nelder-Mead",
               control = list(reltol = 1e-16, maxit = 2e5))
  fit <- optim(fit$par, nll, method = "BFGS",
               control = list(reltol = 1e-16, maxit = 1e5))
  p <- unpack(fit$par)
  list(S1 = p$S1, lam = p$lam, loglik = -fit$value,
       beta = if (link == "log") fit$par[J + 1] else NULL)
}

# Expected (Fisher) information of the paper's Appendix, written cell by cell:
# for a binomial cell with failure probability p = 1 - S_1(IT_j)^lambda_k,
#   I = N (dp/dtheta)(dp/dtheta)^T / {p (1 - p)}.
# Parameters: theta = (S_1(IT_1..IT_J), lambda_2..lambda_K)  or  (S_1, beta).
dp_dtheta <- function(fit, c) {
  st <- fit$st; d <- st$d; J <- st$J; K <- st$K
  k <- d$env[c]; j <- st$g[c]; s <- fit$S1[j]; lk <- fit$lam[k]
  np <- if (fit$link == "log") J + 1 else J + K - 1
  gr <- numeric(np)
  gr[j] <- -lk * s^(lk - 1)
  dl <- -s^lk * log(s)                     # dp / dlambda_k
  if (fit$link == "none") { if (k > 1) gr[J + k - 1] <- dl }
  else gr[J + 1] <- dl * st$delta[k] * lk
  gr
}
info_ph <- function(fit) {
  st <- fit$st; d <- st$d
  np <- if (fit$link == "log") st$J + 1 else st$J + st$K - 1
  I <- matrix(0, np, np)
  for (c in seq_len(nrow(d))) {
    k <- d$env[c]; s <- fit$S1[st$g[c]]^fit$lam[k]; p <- 1 - s
    if (p < 1e-12 || p > 1 - 1e-12) next     # cell on the boundary: no information
    gr <- dp_dtheta(fit, c)
    I <- I + d$N[c] * outer(gr, gr) / (p * (1 - p))
  }
  I
}
par_names <- function(fit) {
  st <- fit$st
  c(sprintf("S1(%g)", st$grid),
    if (fit$link == "log") "beta" else sprintf("lambda%d", 2:st$K))
}
par_vec <- function(fit) c(fit$S1, if (fit$link == "log") fit$beta else fit$lam[-1])

# Observed information by finite differences of the incomplete log-likelihood
obs_info_ph <- function(fit, h = 1e-5) {
  st <- fit$st; th <- par_vec(fit); J <- st$J
  ll <- function(t) {
    S1 <- t[1:J]
    lam <- if (fit$link == "log") exp(st$delta * t[J + 1]) else c(1, t[-(1:J)])
    loglik_ph(st, S1, lam)
  }
  np <- length(th); H <- matrix(0, np, np)
  for (a in 1:np) for (b in a:np) {
    ea <- eb <- numeric(np); ea[a] <- h * max(1, abs(th[a])); eb[b] <- h * max(1, abs(th[b]))
    H[a, b] <- H[b, a] <- (ll(th + ea + eb) - ll(th + ea - eb) -
                           ll(th - ea + eb) + ll(th - ea - eb)) / (4 * ea[a] * eb[b])
  }
  -H
}
score_ph <- function(fit, h = 1e-6) {
  st <- fit$st; th <- par_vec(fit); J <- st$J
  ll <- function(t) {
    lam <- if (fit$link == "log") exp(st$delta * t[J + 1]) else c(1, t[-(1:J)])
    loglik_ph(st, t[1:J], lam)
  }
  sapply(seq_along(th), function(a) {
    e <- numeric(length(th)); e[a] <- h * max(1, abs(th[a]))
    (ll(th + e) - ll(th - e)) / (2 * e[a])
  })
}

# Summary table of a PH fit: MLE, se, three intervals for S_1, plain for others
summary_ph <- function(fit, level = 0.95) {
  I <- info_ph(fit)
  # parameters on the boundary of Omega (S_1 = 1 or 0) are not estimated with
  # an asymptotic variance; they are reported with se = NA
  J <- fit$st$J
  free <- c(fit$S1 < 1 - 1e-10 & fit$S1 > 1e-10, rep(TRUE, nrow(I) - J))
  V <- matrix(NA, nrow(I), ncol(I)); V[free, free] <- solve(I[free, free])
  se <- sqrt(diag(V))
  th <- par_vec(fit); J <- fit$st$J; z <- qnorm(1 - (1 - level) / 2)
  ciS <- ci_S(th[1:J], se[1:J], level)
  tab <- data.frame(par = par_names(fit), mle = th, se = se,
                    plain_lo = th - z * se, plain_hi = th + z * se,
                    loglog_lo = NA, loglog_hi = NA, logit_lo = NA, logit_hi = NA)
  tab[1:J, c("loglog_lo", "loglog_hi")] <- ciS$loglog
  tab[1:J, c("logit_lo", "logit_hi")]   <- ciS$logit
  list(tab = tab, V = V, I = I)
}

# Section 4.3: goodness of fit. Pearson statistic with the binomial variance
# N S (1 - S); df = number of cells - number of free parameters (a nominal
# reference: boundary estimates and active order constraints are not deducted).
gof_ph <- function(fit) {
  st <- fit$st; d <- st$d
  S <- fit$S1[st$g]^fit$lam[d$env]
  e <- d$N * (1 - S)
  # a cell with a fitted probability of 0 or 1 (boundary estimate) reproduces
  # its count exactly at the MLE and contributes zero; this avoids 0/0 when
  # the EM is stopped early (as in the Shiny application)
  v <- d$N * S * (1 - S)
  stat <- sum(ifelse(v > 1e-12, (e - d$r)^2 / pmax(v, 1e-300), 0))
  # the version printed as (17) in the paper, without N in the denominator,
  # applied to proportions: sum (e/N - r/N)^2 / {S (1 - S)}
  w <- S * (1 - S)
  stat_prop <- sum(ifelse(w > 1e-12, (e / d$N - d$r / d$N)^2 / pmax(w, 1e-300), 0))
  np <- if (fit$link == "log") st$J + 1 else st$J + st$K - 1
  df <- nrow(d) - np
  list(stat = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE),
       stat_prop = stat_prop, p_prop = pchisq(stat_prop, df, lower.tail = FALSE))
}

# Prediction at a new covariate value x0 (e.g. a normal storage temperature):
# S_x0(IT_j) = S_1(IT_j)^exp{(x0 - x_1) beta}, delta method on the log-log scale
predict_ph <- function(fit, x0, level = 0.95) {
  stopifnot(fit$link == "log")
  st <- fit$st; J <- st$J; V <- summary_ph(fit, level)$V
  d0 <- x0 - st$xk[1]; lam0 <- exp(d0 * fit$beta)
  S0 <- fit$S1^lam0
  # eta = log(-log S0) = log lam0 + log(-log S1)
  se_eta <- sapply(1:J, function(j) {
    if (!is.finite(V[j, j])) return(NA_real_)
    g <- numeric(J + 1)
    g[j] <- 1 / (fit$S1[j] * log(fit$S1[j]))
    g[J + 1] <- d0
    ok <- g != 0
    sqrt(drop(t(g[ok]) %*% V[ok, ok, drop = FALSE] %*% g[ok]))
  })
  z <- qnorm(1 - (1 - level) / 2)
  eta <- log(-log(S0))
  data.frame(IT = st$grid, S = S0, lo = exp(-exp(eta + z * se_eta)),
             hi = exp(-exp(eta - z * se_eta)))
}

# =============================================================================
# Parametric competitors used in the Introduction of the paper (Tables 4-5):
# exponential and Weibull lifetimes, scale = exp(alpha + b x)
# =============================================================================
fit_param <- function(d, family = c("exponential", "weibull")) {
  family <- match.arg(family)
  Sfun <- function(th) {
    sc <- exp(th[1] + th[2] * d$x)
    if (family == "exponential") exp(-d$IT / sc) else exp(-(d$IT / sc)^exp(th[3]))
  }
  nll <- function(th) { S <- Sfun(th); -sum(xlogy(d$r, 1 - S) + xlogy(d$N - d$r, S)) }
  st <- if (family == "exponential") c(log(mean(d$IT)) + 1, 0) else c(log(mean(d$IT)) + 1, 0, 0)
  best <- NULL
  for (s3 in if (family == "weibull") c(0, -1, -2, -3) else 0) {
    s <- st; if (family == "weibull") s[3] <- s3
    f <- optim(s, nll, method = "BFGS", control = list(reltol = 1e-16, maxit = 1e5))
    f <- optim(f$par, nll, method = "Nelder-Mead", control = list(reltol = 1e-16, maxit = 1e5))
    if (is.null(best) || f$value < best$value) best <- f
  }
  S <- Sfun(best$par)
  e <- d$N * (1 - S)
  list(par = best$par, loglik = -best$value, S = S,
       stat = sum((e - d$r)^2 / (d$N * S * (1 - S))),
       stat_prop = sum(((e - d$r) / d$N)^2 / (S * (1 - S))),
       df = nrow(d) - length(best$par))
}

# =============================================================================
# Output helpers
# =============================================================================
f4 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 4))

# One emitter for every LaTeX table: tabular* flush to both margins
latex_table <- function(df, caption, label, align = NULL) {
  nc <- ncol(df)
  if (is.null(align)) align <- paste0("l", strrep("r", nc - 1))
  cells <- as.matrix(df)
  lines <- c("\\begin{table}[htbp]", "\\centering", "\\footnotesize", "\\setlength{\\tabcolsep}{2pt}",
             sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
             sprintf("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}%s@{}}", align),
             "\\hline",
             paste0(paste(colnames(df), collapse = " & "), " \\\\"), "\\hline",
             apply(cells, 1, function(r) paste0(paste(r, collapse = " & "), " \\\\")),
             "\\hline", "\\end{tabular*}", "\\end{table}")
  paste(lines, collapse = "\n")
}

# Figures: written as .eps and .pdf through one helper, rendered to a local
# temporary file and copied into figures/
save_figure <- function(name, plot_fun, width = 6, height = 4.5) {
  dir.create("figures", showWarnings = FALSE)
  for (ext in c("eps", "pdf")) {
    tmp <- tempfile(fileext = paste0(".", ext))
    # cairo devices embed the fonts, as journals require for vector artwork;
    # the standard devices are used only when cairo is not available
    cairo <- isTRUE(capabilities("cairo"))
    if (ext == "eps") {
      if (cairo) cairo_ps(tmp, width = width, height = height, fallback_resolution = 600)
      else postscript(tmp, width = width, height = height,
                      horizontal = FALSE, onefile = FALSE, paper = "special")
    } else {
      if (cairo) cairo_pdf(tmp, width = width, height = height) else pdf(tmp, width = width, height = height)
    }
    par(mar = c(4.2, 4.2, 1, 1), cex = 1)
    plot_fun()
    dev.off()
    file.copy(tmp, file.path("figures", paste0(name, ".", ext)), overwrite = TRUE)
  }
  invisible(NULL)
}

# =============================================================================
# Weibull hazard method (Gacula, 1975; Cardelli and Labuza, 2001) as applied by
# Keklik, Isikli and Sur (2017): the unacceptable responses are ranked in
# reverse order, hazard = 100 / rank, cumulative hazard CH summed over the
# ranked responses; log10(t) is regressed on log10(CH) using the points with
# CH <= CH_max (100 in the source article) and the shelf life is the time at
# which CH = 69.3, with the 95% confidence limits of the fitted mean.
# =============================================================================
gacula_fit <- function(IT, r, CH_max = 100, level = 0.95, N = NULL) {
  # N = NULL: the version of Keklik et al. (2017), in which only the
  #   unacceptable responses are ranked.
  # N given: the original formulation of Gacula (1975): all responses are ranked in reverse order, the acceptable ones
  #   acting as censored observations (failures listed first at a tie), and
  #   the hazard 100 / reverse rank is accumulated over the failures only.
  if (is.null(N)) {
    t <- rep(IT, r); n <- length(t)
    if (n < 3) return(NULL)
    CH <- cumsum(100 / (n:1))
  } else {
    tt <- unlist(mapply(function(d, n, rr) rep(d, n), IT, N, r))
    fail <- unlist(mapply(function(n, rr) c(rep(TRUE, rr), rep(FALSE, n - rr)), N, r))
    K <- length(tt):1
    CH <- cumsum(ifelse(fail, 100 / K, 0))[fail]; t <- tt[fail]
    if (length(t) < 3) return(NULL)
  }
  keep <- CH <= CH_max + 1e-9
  if (sum(keep) < 3 || length(unique(t[keep])) < 2) return(NULL)
  dd <- data.frame(x = log10(CH[keep]), y = log10(t[keep]))
  f <- lm(y ~ x, data = dd)
  p <- predict(f, data.frame(x = log10(69.3)), interval = "confidence", level = level)
  a <- unname(coef(f)[1]); b <- unname(coef(f)[2])
  list(beta = 1 / b, t50 = 10^p[1], lo = 10^p[2], hi = 10^p[3],
       r2 = summary(f)$r.squared, n_used = sum(keep), t = t, CH = CH,
       S = function(tt) exp(-10^((log10(tt) - a) / b) / 100))
}
