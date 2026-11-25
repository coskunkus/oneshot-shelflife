# =============================================================================
# run_pezik.R
# Zhu and Balakrishnan (2022) applied to the staggered-sampling (Gacula, 1975)
# sensory data of Keklik, Isikli and Sur (2017), Food Sci. Technol. (Campinas)
# 37(suppl 1), doi:10.1590/1678-457X.33216, Tables 1 and 2 (CC BY 4.0).
# At each test day N_j panellists label a sample acceptable (+) or not (-):
# every evaluation is a one-shot observation. Groups: blanching 7 and 9 min.
# Writes results/run_pezik_report.txt and figures/pezik_survival.{eps,pdf}.
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
rep_lines <- character(0)
out <- function(...) rep_lines <<- c(rep_lines, sprintf(...))
blk <- function(x) rep_lines <<- c(rep_lines, capture.output(x))
hdr <- function(s) out("\n%s\n%s", s, strrep("=", nchar(s)))

# Tables 1 and 2 of the source article (counts of + and - per test day)
day1 <- c(1, 20, 40, 60, 80, 100, 120, 140, 160, 170)
N1   <- c(3, 4, 5, 6, 7, 8, 9, 10, 11, 19)
r1   <- c(0, 0, 0, 2, 1, 1, 2, 3, 7, 12)            # 28 unacceptable in total
day2 <- c(1, 20, 40, 60, 80, 100, 120, 140, 150)
N2   <- c(3, 4, 5, 6, 7, 8, 9, 10, 18)
r2   <- c(0, 0, 1, 2, 1, 2, 4, 7, 13)               # 30 unacceptable in total
d <- os_data(env = c(rep(1, 10), rep(2, 9)), x = c(rep(7, 10), rep(9, 9)),
             IT = c(day1, day2), N = c(N1, N2), r = c(r1, r2))

out("run_pezik_report.txt  --  %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
out("Kaynak veri: Keklik, Isikli, Sur (2017), Tablo 1-2; toplam '-' sayilari 28 ve 30 (makaleyle ayni).")
hdr("Veri (env 1 = 7 dk haslama, env 2 = 9 dk haslama)"); blk(d)

hdr("Bolum 2: her grup icin parametrik olmayan MLE ve %95 GA (logit)")
for (k in 1:2) {
  s <- d[d$env == k, ]; p <- pava_np(s$N, s$r); e <- em_np(s$IT, s$N, s$r)
  ci <- ci_S(p, se_np(p, s$N))
  out("Grup %d (%d dk): EM %d iterasyon, EM-PAVA logL farki %.2e", k, c(7, 9)[k], e$iter,
      loglik_np(p, s$N, s$r) - e$loglik)
  blk(data.frame(day = s$IT, N = s$N, r = s$r, raw = round(1 - s$r / s$N, 4), S = round(p, 4),
                 lo = round(ci$logit[, 1], 4), hi = round(ci$logit[, 2], 4)))
}

hdr("Bolum 4: PH modeli, S_9dk(t) = S_7dk(t)^lambda")
f <- fit_ph(d, "none"); g <- direct_ph(d, "none"); s <- summary_ph(f)
blk(format(s$tab, digits = 4))
out("Yol: %s; logL = %.6f; dogrudan maksimizasyon logL = %.6f", f$route, f$loglik, g$loglik)
lam <- f$lam[2]; se_lam <- s$tab$se[nrow(s$tab)]
out("lambda (9 dk / 7 dk tehlike orani) = %.4f; log-olcekte %%95 GA = (%.4f, %.4f)", lam,
    exp(log(lam) - z975 * se_lam / lam), exp(log(lam) + z975 * se_lam / lam))
gf <- gof_ph(f)
out("Uyum iyiligi: Pearson = %.4f, sd = %d, p = %.4f", gf$stat, gf$df, gf$p)
f_null <- fit_ph(transform(d, env = 1L), "none")
out("Sayisal ayrinti: baslangic S1(IT_j) = 1 - j/(J+1) (J = %d), lambda = 1; EM %d iterasyon (%s); logL = %.6f; lambda = 1 altinda logL = %.6f",
    f$st$J, f$iter, f$route, f$loglik, f_null$loglik)
out("Sayisal ayrinti: degerlendirme sayisi %d (7 dk) ve %d (9 dk); ret sayisi %d ve %d",
    sum(d$N[d$env == 1]), sum(d$N[d$env == 2]), sum(d$r[d$env == 1]), sum(d$r[d$env == 2]))
# reduced model: no blanching effect (lambda = 1)
fr <- fit_ph(transform(d, env = 1L), "none")
lr <- 2 * (f$loglik - fr$loglik)
out("Haslama etkisi (LR testi, lambda = 1): %.4f, chi2(1) p = %.4f", lr, pchisq(lr, 1, lower.tail = FALSE))

hdr("Raf omru: S(t) = 0.5 ve 0.9 icin ilk muayene gunu (NPMLE basamaklari)")
for (k in 1:2) {
  Sk <- f$S1^f$lam[k]; tt <- f$st$grid
  out("Grup %d: S <= 0.9 ilk gun = %s; S <= 0.5 ilk gun = %s", k,
      tt[which(Sk <= 0.9)[1]], tt[which(Sk <= 0.5)[1]])
}
out("Karsilastirma: makale (Weibull hazard, Gacula 1975) %%50 reddi 167 gun (7 dk) ve 148 gun (9 dk) verir.")

save_figure("pezik_survival", function() {
  plot(NA, xlim = c(0, 175), ylim = c(0, 1), xlab = "Storage time (days)",
       ylab = "Probability of acceptance")
  for (k in 1:2) {
    lines(stepfun(f$st$grid, c(1, f$S1^f$lam[k])), do.points = FALSE, lty = k)
    s <- d[d$env == k, ]; points(s$IT, 1 - s$r / s$N, pch = k)
  }
  legend("bottomleft", legend = c("Blanched 7 min", "Blanched 9 min"), lty = 1:2, pch = 1:2, bty = "n")
})
writeLines(rep_lines, file.path("results", "run_pezik_report.txt"))
cat("run_pezik.R finished\n")

# ------------------------------------------------ Weibull hazard method ---
rep_lines <- readLines(file.path("results", "run_pezik_report.txt"))
hdr("Kaynak makalenin yontemi (Weibull hazard, Gacula 1975) yeniden uretimi")
for (k in 1:2) {
  s <- d[d$env == k, ]
  for (chm in c(Inf, 100)) {
    g <- gacula_fit(s$IT, s$r, CH_max = chm)
    out("Grup %d (%d dk), CH <= %s: %d nokta, beta = %.2f, R2 = %.2f, t50 = %.1f gun, %%95 GA (%.1f, %.1f)",
        k, c(7, 9)[k], format(chm), g$n_used, g$beta, g$r2, g$t50, g$lo, g$hi)
  }
}
hdr("Gacula (1975) ozgun hali: kabul edilenler dahil tum numuneler ters sirada siralanir")
for (k in 1:2) {
  s <- d[d$env == k, ]
  for (chm in c(Inf, 100)) {
    g <- gacula_fit(s$IT, s$r, CH_max = chm, N = s$N)
    out("Grup %d (%d dk), CH <= %s: %d nokta, beta = %.2f, R2 = %.2f, t50 = %.1f gun, %%95 GA (%.1f, %.1f); S(%d) = %.4f",
        k, c(7, 9)[k], format(chm), g$n_used, g$beta, g$r2, g$t50, g$lo, g$hi,
        c(160, 140)[k], g$S(c(160, 140)[k]))
  }
}
out("Makale: tum noktalarla beta = 4.27 (R2 = 0.78) ve 3.90 (R2 = 0.79); CH <= 100 ile beta = 2.84 ve 2.58,")
out("t50 = 167 (156-179) ve 148 (138-159) gun. Degerler yeniden elde edildi (GA siniri en fazla 1 gun farkli,")
out("makaledeki H degerlerinin bir ondaliga yuvarlanmasindan).")
hdr("Ayni veride iki yaklasimin S(t) tahminleri")
for (k in 1:2) {
  s <- d[d$env == k, ]; g <- gacula_fit(s$IT, s$r); p <- pava_np(s$N, s$r)
  blk(data.frame(group = c(7, 9)[k], day = s$IT, raw = round(1 - s$r / s$N, 4),
                 NPMLE = round(p, 4), Weibull_hazard = round(g$S(s$IT), 4)))
}
writeLines(rep_lines, file.path("results", "run_pezik_report.txt"))

# ------------------------------------------------------------ LaTeX tables ---
rep_lines <- readLines(file.path("results", "run_pezik_report.txt"))
hdr("LaTeX tablolari")
grid <- sort(unique(d$IT))
cell <- function(k, v) sapply(grid, function(t) { i <- which(d$env == k & d$IT == t)
  if (length(i)) format(d[[v]][i]) else "--" })
tab_data <- data.frame(Day = grid, `$N$ (7 min)` = cell(1, "N"), `Unacceptable (7 min)` = cell(1, "r"),
                       `$N$ (9 min)` = cell(2, "N"), `Unacceptable (9 min)` = cell(2, "r"), check.names = FALSE)
out("%s\n", latex_table(tab_data,
  "Staggered sensory data of pezik pickles: number of evaluations and number of unacceptable responses on each test day.",
  "tab:data", align = "rrrrr"))
est_rows <- list()
for (k in 1:2) {
  s <- d[d$env == k, ]; p <- pava_np(s$N, s$r); ci <- ci_S(p, se_np(p, s$N)); g <- gacula_fit(s$IT, s$r)
  est_rows[[k]] <- data.frame(Group = c("7 min", "9 min")[k], Day = s$IT, Observed = f4(1 - s$r / s$N),
    NPMLE = f4(p), `95\\% CI (logit)` = ifelse(is.finite(ci$logit[, 1]),
      sprintf("(%s, %s)", f4(ci$logit[, 1]), f4(ci$logit[, 2])), "--"),
    `Weibull hazard` = f4(g$S(s$IT)), check.names = FALSE)
}
out("%s\n", latex_table(do.call(rbind, est_rows),
  "Estimated probability of acceptance on each test day: observed proportion, non-parametric MLE with its 95\\% confidence interval, and the fit of the Weibull hazard method.",
  "tab:npmle", align = "lrrrrr"))
tb <- s_ph <- summary_ph(f)$tab
ph_rows <- data.frame(Parameter = c(sprintf("$S_1(%g)$", f$st$grid), "$\\lambda$"),
  MLE = f4(tb$mle), SE = f4(tb$se),
  `95\\% CI` = c(ifelse(is.finite(tb$logit_lo[1:f$st$J]),
                        sprintf("(%s, %s)", f4(tb$logit_lo[1:f$st$J]), f4(tb$logit_hi[1:f$st$J])), "--"),
                 sprintf("(%s, %s)", f4(exp(log(lam) - z975 * se_lam / lam)), f4(exp(log(lam) + z975 * se_lam / lam)))),
  check.names = FALSE)
out("%s\n", latex_table(ph_rows,
  "Proportional hazards model for the two blanching times: MLEs, standard errors (SE) and 95\\% confidence intervals (logit transformation for $S_1$, log transformation for $\\lambda$).",
  "tab:ph", align = "lrrr"))
out("Grup 2'nin ayri 140. gun tahmini, temel olcege donusturulmus: 0.3000^(1/lambda) = %.4f; PH S1(140) GA genisligi %.4f, grup 1 tek basina GA genisligi %.4f",
    0.3^(1 / lam), diff(unlist(tb[tb$par == "S1(140)", c("logit_lo", "logit_hi")])),
    diff(as.numeric(ci_S(0.7, se_np(0.7, 10))$logit)))
fs <- fit_ph(d, "none", it_em = 100)
out("Shiny yolu (it_em = 100, %s): Pearson %.4f, df %d, p %.4f", fs$route, gof_ph(fs)$stat, gof_ph(fs)$df, gof_ph(fs)$p)
writeLines(rep_lines, file.path("results", "run_pezik_report.txt"))

# ----------------------------------------------------------------- figures ---
# Figure 1: non-parametric MLE with 95% logit intervals, one panel per group.
# Figure 2: proportional hazards fit for both groups.
# Figure 3: Weibull hazard method against the one-shot analysis (7 min group).
col_g <- c("#2a78d6", "#eb6834")       # group 1 (7 min), group 2 (9 min)
ink <- "#3d3d3a"
grp_lab <- c("Blanched 7 min", "Blanched 9 min")
fig_np <- function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 0.8))
  for (k in 1:2) {
    s <- d[d$env == k, ]; p <- pava_np(s$N, s$r); ci <- ci_S(p, se_np(p, s$N))
    plot(NA, xlim = c(0, 175), ylim = c(0, 1), xlab = "Storage time (days)",
         ylab = "Probability of acceptance", las = 1)
    title(main = sprintf("(%s) %s", letters[k], grp_lab[k]), font.main = 1, cex.main = 1)
    ok <- is.finite(ci$logit[, 1])
    segments(s$IT[ok], ci$logit[ok, 1], s$IT[ok], ci$logit[ok, 2], col = col_g[k], lwd = 1.5)
    segments(s$IT[ok] - 2.5, ci$logit[ok, 1], s$IT[ok] + 2.5, ci$logit[ok, 1], col = col_g[k], lwd = 1.5)
    segments(s$IT[ok] - 2.5, ci$logit[ok, 2], s$IT[ok] + 2.5, ci$logit[ok, 2], col = col_g[k], lwd = 1.5)
    lines(stepfun(s$IT, c(1, p)), do.points = FALSE, col = col_g[k], lwd = 2, verticals = TRUE)
    points(s$IT, 1 - s$r / s$N, pch = 21, bg = "white", col = ink, cex = 0.6 + 0.12 * sqrt(s$N) * 2)
    legend("bottomleft", c("Non-parametric MLE", "95% CI (logit)", "Observed proportion"),
           lty = c(1, 1, NA), lwd = c(2, 1.5, NA), pch = c(NA, NA, 21), pt.bg = "white",
           col = c(col_g[k], col_g[k], ink), bty = "n", cex = 0.85)
  }
}
save_figure("pezik_fig1_npmle", fig_np, width = 9, height = 4.3)

fig_ph <- function() {
  par(mar = c(4.2, 4.2, 1, 1))
  plot(NA, xlim = c(0, 175), ylim = c(0, 1), xlab = "Storage time (days)",
       ylab = "Probability of acceptance", las = 1)
  for (k in 1:2) {
    lines(stepfun(f$st$grid, c(1, f$S1^f$lam[k])), do.points = FALSE, col = col_g[k],
          lwd = 2, lty = k, verticals = TRUE)
    s <- d[d$env == k, ]
    points(s$IT, 1 - s$r / s$N, pch = c(21, 24)[k], bg = "white", col = col_g[k], cex = 1.2, lwd = 1.5)
  }
  legend("bottomleft", c(sprintf("%s: fitted", grp_lab), sprintf("%s: observed", grp_lab)),
         lty = c(1, 2, NA, NA), lwd = c(2, 2, 1.5, 1.5), pch = c(NA, NA, 21, 24), pt.bg = "white",
         col = rep(col_g, 2), bty = "n", cex = 0.9)
}
save_figure("pezik_fig2_ph", fig_ph, width = 6.5, height = 4.6)

fig_wh <- function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 0.8))
  s <- d[d$env == 1, ]; tt <- rep(s$IT, s$r); n <- length(tt); CH <- cumsum(100 / (n:1))
  gA <- gacula_fit(s$IT, s$r, CH_max = Inf); gT <- gacula_fit(s$IT, s$r, CH_max = 100)
  plot(CH, tt, log = "xy", xlim = c(3, 450), ylim = c(55, 230), las = 1,
       xlab = "Cumulative hazard (%)", ylab = "Storage time (days)",
       pch = 21, bg = ifelse(CH <= 100, ink, "white"), col = ink)
  title(main = "(a) Weibull hazard plot, 7 min", font.main = 1, cex.main = 1)
  xs <- exp(seq(log(3), log(450), length.out = 100))
  for (g in list(list(gA, 2), list(gT, 1))) {
    b <- 1 / g[[1]]$beta; a <- log10(g[[1]]$t50) - b * log10(69.3)
    lines(xs, 10^(a + b * log10(xs)), col = col_g[2], lwd = 2, lty = g[[2]])
  }
  abline(v = 69.3, col = "grey60", lty = 3)
  points(c(69.3, 69.3), c(gT$t50, gA$t50), pch = 4, col = col_g[2], cex = 1.3, lwd = 2)
  segments(69.3, gT$t50, 115, 118, col = "grey50"); text(118, 118, sprintf("%.0f days", gT$t50), adj = 0, cex = 0.85, col = ink)
  segments(69.3, gA$t50, 115, 100, col = "grey50"); text(118, 100, sprintf("%.0f days", gA$t50), adj = 0, cex = 0.85, col = ink)
  legend("topleft", c("Fit to points up to 100%", "Fit to all points", "Point used in the 100% fit", "Point not used"),
         lty = c(1, 2, NA, NA), lwd = c(2, 2, NA, NA), pch = c(NA, NA, 21, 21), pt.bg = c(NA, NA, ink, "white"),
         col = c(col_g[2], col_g[2], ink, ink), bty = "n", cex = 0.8)
  # (b) probability of acceptance: data, one-shot MLE and Weibull hazard fit
  p <- pava_np(s$N, s$r)
  plot(NA, xlim = c(0, 175), ylim = c(0, 1), xlab = "Storage time (days)",
       ylab = "Probability of acceptance", las = 1)
  title(main = "(b) Probability of acceptance, 7 min", font.main = 1, cex.main = 1)
  curve(gT$S(x), from = 1, to = 175, add = TRUE, col = col_g[2], lwd = 2, lty = 2)
  lines(stepfun(s$IT, c(1, p)), do.points = FALSE, col = col_g[1], lwd = 2, verticals = TRUE)
  points(s$IT, 1 - s$r / s$N, pch = 21, bg = "white", col = ink, cex = 1.1)
  abline(h = 0.5, col = "grey60", lty = 3)
  legend("bottomleft", c("Non-parametric MLE", "Weibull hazard fit", "Observed proportion"),
         lty = c(1, 2, NA), lwd = c(2, 2, NA), pch = c(NA, NA, 21), pt.bg = "white",
         col = c(col_g[1], col_g[2], ink), bty = "n", cex = 0.85)
}
save_figure("pezik_fig3_weibull", fig_wh, width = 9, height = 4.3)
cat("figures written\n")
