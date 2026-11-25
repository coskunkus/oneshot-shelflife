# =============================================================================
# check_formulas.R
# Independent numerical verification of the formulas implemented in
# core_oneshot.R. Layer 1: exact checks. Layer 2: Monte Carlo checks with a
# fixed seed, reported as z-scores against the Monte Carlo standard error
# (|z| < 4 passes). Writes results/check_formulas_report.txt.
# =============================================================================
locate_root <- function(core = "core_oneshot.R") {
  cand <- character(0)
  for (i in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(f)) cand <- c(cand, dirname(normalizePath(f)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (is.character(p) && length(p) == 1 && nzchar(p)) cand <- c(cand, dirname(normalizePath(p)))
  }
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) cand <- c(cand, dirname(normalizePath(sub("^--file=", "", fa[1]))))
  cand <- c(cand, getwd())
  for (cd in cand) for (up in c(cd, dirname(cd))) {
    if (file.exists(file.path(up, "code", core))) return(up)
    if (file.exists(file.path(up, core)) && basename(up) == "code") return(dirname(up))
  }
  stop("project root not found")
}
setwd(locate_root())
source(file.path("code", "core_oneshot.R"))
dir.create("results", showWarnings = FALSE)

res <- data.frame(check = character(0), value = numeric(0), tol = numeric(0),
                  pass = logical(0))
chk <- function(name, value, tol) {
  res <<- rbind(res, data.frame(check = name, value = value, tol = tol,
                                pass = is.finite(value) && abs(value) < tol))
}
set.seed(20261122)

temp <- c(35, 45, 55, 65)
d1 <- os_data(env = rep(1:4, 3), x = rep(temp, 3), IT = rep(c(10, 20, 30), each = 4),
              N = 10, r = c(2, 5, 5, 4, 0, 3, 3, 4, 1, 2, 6, 8))
d2 <- os_data(env = rep(1:2, each = 3), x = rep(c(60, 400), each = 3),
              IT = rep(c(9.37, 14.07, 18.70), 2),
              N = c(72, 49, 46, 25, 35, 11), r = c(2, 1, 11, 3, 21, 10))

# ---------------------------------------------------------- exact layer ---
# E1. EM of Section 2.1 reaches the NPMLE (isotonic regression) on random data
for (b in 1:20) {
  J <- 5; N <- sample(5:30, J, TRUE); r <- rbinom(J, N, sort(runif(J, 0.05, 0.9)))
  e <- em_np(1:J, N, r); p <- pava_np(N, r)
  chk(sprintf("E1.%02d EM vs PAVA, logL difference", b),
      loglik_np(p, N, r) - e$loglik, 1e-7)
}
# E2. Lemma 1: with monotone proportions the MLE is 1 - r/N
N <- c(20, 20, 20); r <- c(2, 7, 15)
chk("E2 Lemma 1 closed form", max(abs(pava_np(N, r) - (1 - r / N))), 1e-14)
# E3. Paper Table 12 intervals reproduced to 4 decimals
ci <- ci_S(0.9, se_np(0.9, 10))
chk("E3 Table 12 plain lower", round(ci$plain[1], 4) - 0.7141, 1e-12)
chk("E3 Table 12 log-log lower", round(ci$loglog[1], 4) - 0.4730, 1e-12)
chk("E3 Table 12 log-log upper", round(ci$loglog[2], 4) - 0.9853, 1e-12)
chk("E3 Table 12 logit lower", round(ci$logit[1], 4) - 0.5328, 1e-12)
chk("E3 Table 12 logit upper", round(ci$logit[2], 4) - 0.9861, 1e-12)
# E4. PH EM equals direct maximisation (log-likelihood) on the two examples
for (cs in list(list("Ex1 none", d1, "none"), list("Ex1 log", d1, "log"),
                list("Ex2 log", d2, "log"))) {
  f <- em_ph(cs[[2]], cs[[3]]); g <- direct_ph(cs[[2]], cs[[3]])
  # one-sided on purpose: the EM must reach a log-likelihood at least as high as
  # the direct maximization (a lower value of the direct route is not an error of the EM)
  chk(sprintf("E4 %s: logL(direct) - logL(EM)", cs[[1]]), max(0, g$loglik - f$loglik), 1e-8)
}
# E5. PH EM with K = 1 reduces to the NPMLE of Section 2
s35 <- d1[d1$env == 1, ]; s35$r <- c(1, 4, 6)
f1 <- em_ph(s35, "none")
chk("E5 PH(K=1) vs PAVA", max(abs(f1$S1 - pava_np(s35$N, s35$r))), 1e-6)
# E6-E7. Expected information (Appendix) against exact enumeration of all
# outcomes of a tiny design: E[-Hessian] = I and E[score score^T] = I
dt <- os_data(env = c(1, 1, 2, 2), x = c(0, 0, 1, 1), IT = c(1, 2, 1, 2),
              N = c(3, 3, 3, 3), r = c(0, 0, 0, 0))
for (lk in c("none", "log")) {
  fit <- list(st = ph_setup(dt), link = lk, S1 = c(0.8, 0.55),
              lam = c(1, 1.7), beta = log(1.7))
  I <- info_ph(fit)
  grid <- as.matrix(expand.grid(0:3, 0:3, 0:3, 0:3))
  EH <- matrix(0, nrow(I), ncol(I)); ES <- EH
  pcell <- 1 - fit$S1[fit$st$g]^fit$lam[dt$env]
  for (i in seq_len(nrow(grid))) {
    rr <- grid[i, ]; pr <- prod(dbinom(rr, 3, pcell))
    ff <- fit; ff$st$d$r <- rr
    EH <- EH + pr * obs_info_ph(ff, h = 1e-4)
    s <- score_ph(ff); ES <- ES + pr * outer(s, s)
  }
  chk(sprintf("E6 [%s] max rel |E(-H) - I|", lk), max(abs(EH - I)) / max(abs(I)), 1e-5)
  chk(sprintf("E7 [%s] max rel |E(ss') - I|", lk), max(abs(ES - I)) / max(abs(I)), 1e-6)
}
# E8. Delta-method gradient in predict_ph against a numerical gradient
f14 <- em_ph(d1, "log"); x0 <- 20; V <- solve(info_ph(f14)); J <- 3
eta <- function(th) (x0 - 35) * th[J + 1] + log(-log(th[1:J]))
th <- par_vec(f14)
G <- sapply(1:(J + 1), function(a) { e <- numeric(J + 1); e[a] <- 1e-6
  (eta(th + e) - eta(th - e)) / 2e-6 })
se_num <- sqrt(diag(G %*% V %*% t(G)))
pr <- predict_ph(f14, x0); z <- qnorm(0.975)
se_core <- (log(-log(pr$lo)) - log(-log(pr$S))) / z
chk("E8 delta-method se, numerical vs core", max(abs(se_num - se_core) / se_num), 1e-6)
# E9. Pearson statistic reproduces the p-value of Table 14 to 3 decimals
chk("E9 Table 14 GOF p-value (|diff| < 0.001)", gof_ph(f14)$p - 0.1889, 1e-3)

# E10. The two food analyses of the article: the EM route of fit_ph and a direct
# maximization started from a different point (equal increments of the
# cumulative baseline hazard, lambda = 1, beta = 0) reach the same maximum of
# the log-likelihood and the same MLEs (two-sided checks).
pz <- os_data(env = c(rep(1, 10), rep(2, 9)), x = c(rep(7, 10), rep(9, 9)),
              IT = c(1, 20, 40, 60, 80, 100, 120, 140, 160, 170, 1, 20, 40, 60, 80, 100, 120, 140, 150),
              N = c(3, 4, 5, 6, 7, 8, 9, 10, 11, 19, 3, 4, 5, 6, 7, 8, 9, 10, 18),
              r = c(0, 0, 0, 2, 1, 1, 2, 3, 7, 12, 0, 0, 1, 2, 1, 2, 4, 7, 13))
gt_IT <- list(seq(5, 70, 5), seq(5, 55, 5), seq(5, 45, 5))
gt_N  <- list(c(2:10, 10, 10, 10, 10, 10), c(2:8, 9, 10, 10, 10), c(2:5, 6, 9, 10, 10, 10))
gt_r  <- list(c(0,0,0,0,0,0,0,0,0,2,0,1,3,6), c(0,0,0,0,0,0,0,2,2,5,8), c(0,0,0,0,2,2,2,3,8))
gt <- do.call(rbind, lapply(1:3, function(k)
  os_data(env = k, x = c(20, 30, 40)[k], IT = gt_IT[[k]], N = gt_N[[k]], r = gt_r[[k]])))
for (cs in list(list("pezik, free lambda", pz, "none"), list("green tea, free lambda", gt, "none"),
                list("green tea, log link", gt, "log"))) {
  f <- fit_ph(cs[[2]], cs[[3]]); g <- direct_ph(cs[[2]], cs[[3]])
  chk(sprintf("E10 %s: |logL(fit_ph) - logL(direct)|", cs[[1]]), f$loglik - g$loglik, 1e-6)
  pf <- c(f$S1, if (cs[[3]] == "log") f$beta else f$lam[-1])
  pg <- c(g$S1, if (cs[[3]] == "log") g$beta else g$lam[-1])
  chk(sprintf("E10 %s: max |MLE diff| / max(1, |MLE|)", cs[[1]]), max(abs(pf - pg) / pmax(1, abs(pf))), 1e-3)
}

# ----------------------------------------------------- Monte Carlo layer ---
# Setting of Section 5 (Tables 10-11): exponential lifetimes with mean
# exp(6 - 0.05 x), x = 35, 45, 55, 65, IT = 10, 20, 30, N_kj = 80.
B <- 400; Nkj <- 80
S1_true <- exp(-c(10, 20, 30) / exp(6 - 0.05 * 35)); beta_true <- 0.05
cover <- matrix(NA, B, 4); est <- matrix(NA, B, 4)
for (b in 1:B) {
  db <- d1; db$N <- Nkj
  db$r <- rbinom(nrow(db), Nkj, 1 - S1_true[match(db$IT, c(10, 20, 30))]^exp(beta_true * (db$x - 35)))
  f <- em_ph(db, "log"); s <- summary_ph(f)$tab
  est[b, ] <- s$mle
  cover[b, 1:3] <- s$logit_lo[1:3] <= S1_true & S1_true <= s$logit_hi[1:3]
  cover[b, 4] <- s$plain_lo[4] <= beta_true & beta_true <= s$plain_hi[4]
}
truth <- c(S1_true, beta_true); nm <- c("S1(10)", "S1(20)", "S1(30)", "beta")
for (a in 1:4) {
  zc <- (mean(cover[, a]) - 0.95) / sqrt(0.95 * 0.05 / B)
  chk(sprintf("M1 coverage %s = %.4f, z", nm[a], mean(cover[, a])), zc, 4)
  zb <- (mean(est[, a]) - truth[a]) / (sd(est[, a]) / sqrt(B))
  chk(sprintf("M2 mean %s = %.4f (true %.4f), z", nm[a], mean(est[, a]), truth[a]), zb, 4)
}

# ----------------------------------------------------------------- report ---
res$status <- ifelse(res$pass, "PASS", "FAIL")
lines <- c(sprintf("check_formulas_report.txt  --  %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
           sprintf("Monte Carlo: B = %d, seed 20261122", B), "",
           sprintf("%-55s %14s %10s  %s", "check", "value", "tol", "status"),
           sprintf("%-55s %14.3e %10.1e  %s", res$check, res$value, res$tol, res$status),
           "", sprintf("SUMMARY: %d PASS, %d FAIL", sum(res$pass), sum(!res$pass)))
writeLines(lines, file.path("results", "check_formulas_report.txt"))
cat(tail(lines, 1), "\n")
