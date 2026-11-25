# =============================================================================
# run_greentea.R
# Zhu and Balakrishnan (2022) applied to the staggered acceptability data of
# Xiao et al. (2022), J. Chem. 2022:2955839, doi:10.1155/2022/2955839, Table 4
# (CC BY 4.0): green tea stored at 20, 30 and 40 C, tested every 5 days.
# Writes results/run_greentea_report.txt and figures/greentea_fig*.{eps,pdf}.
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

# Table 4 of the source article: accepted (ok) and rejected (r) marks per test
temp <- c(20, 30, 40)
days <- list(seq(5, 70, 5), seq(5, 55, 5), seq(5, 45, 5))
r_l  <- list(c(0,0,0,0,0,0,0,0,0,2,0,1,3,6), c(0,0,0,0,0,0,0,2,2,5,8), c(0,0,0,0,2,2,2,3,8))
N_l  <- list(c(2:10, 10,10,10,10,10), c(2:8, 9,10,10,10), c(2:5, 6, 9,10,10,10))
d <- do.call(rbind, lapply(1:3, function(k)
  os_data(env = k, x = temp[k], IT = days[[k]], N = N_l[[k]], r = r_l[[k]])))
rownames(d) <- NULL

out("run_greentea_report.txt  --  %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
out("Kaynak veri: Xiao ve ark. (2022), Tablo 4. Toplam ret sayilari: %s.",
    paste(tapply(d$r, d$env, sum), collapse = ", "))
out("Not: kaynak makalenin Tablo 1'i 30 C icin 18 ret listeler (50. gunde 6), Tablo 4 ise 17 ret");
out("(50. gunde 5) gosterir; 20 C icin Tablo 1 son 6 reddi 70. gun yerine 67. gune yazar. Bu analiz Tablo 4'u kullanir.")
hdr("Veri"); blk(d)

hdr("Kaynak makalenin Weibull hazard analizi (Tablo 1 yaslari ile)")
t1 <- list(c(50,50,60,65,65,65,rep(67,6)), c(40,40,45,45,rep(50,6),rep(55,8)), c(25,25,30,30,35,35,40,40,40,rep(45,8)))
for (k in 1:3) {
  tab <- table(t1[[k]]); g <- gacula_fit(as.numeric(names(tab)), as.integer(tab), CH_max = Inf)
  out("%d C: beta = %.2f, t50 = %.1f gun (%%95 GA %.1f-%.1f)", temp[k], g$beta, g$t50, g$lo, g$hi)
}
out("Makale: beta = 3.37, 3.13, 3.42; FLGT = 67, 57, 42 gun.")
hdr("Ayni yontem Tablo 4 verisiyle")
gac4 <- list()
for (k in 1:3) {
  s <- d[d$env == k, ]; g <- gacula_fit(s$IT, s$r, CH_max = Inf); gac4[[k]] <- g
  gv <- gacula_fit(s$IT, s$r, CH_max = Inf, N = s$N)
  out("%d C: yalniz ret siralamasi beta = %.2f, t50 = %.1f gun; tum yanitlar siralaninca t50 = %.1f gun",
      temp[k], g$beta, g$t50, gv$t50)
}

hdr("Her sicaklik icin parametrik olmayan MLE ve logit %95 GA")
np <- list()
for (k in 1:3) {
  s <- d[d$env == k, ]; p <- pava_np(s$N, s$r); ci <- ci_S(p, se_np(p, s$N))
  np[[k]] <- data.frame(temp = temp[k], day = s$IT, N = s$N, r = s$r, raw = 1 - s$r / s$N, S = p,
                        lo = ci$logit[, 1], hi = ci$logit[, 2], WH = gac4[[k]]$S(s$IT))
  blk(format(np[[k]], digits = 4))
}

hdr("PH modelleri")
f0 <- fit_ph(d, "none"); f1 <- fit_ph(d, "log")
fr <- fit_ph(transform(d, env = 1L), "none")
s0 <- summary_ph(f0)$tab; s1 <- summary_ph(f1)$tab
out("Baglantisiz: lambda_30 = %.4f, lambda_40 = %.4f; logL = %.4f; yol %s", f0$lam[2], f0$lam[3], f0$loglik, f0$route)
out("Log-baglanti lambda_k = exp{beta (x_k - 20)}: beta = %.4f (SE %.4f); logL = %.4f; yol %s",
    f1$beta, s1$se[nrow(s1)], f1$loglik, f1$route)
be <- f1$beta; seb <- s1$se[nrow(s1)]
out("beta %%95 GA = (%.4f, %.4f); 10 C artista tehlike orani exp(10 beta) = %.4f (%.4f, %.4f)", be - z975 * seb, be + z975 * seb,
    exp(10 * be), exp(10 * (be - z975 * seb)), exp(10 * (be + z975 * seb)))
lr_t <- 2 * (f1$loglik - fr$loglik); lr_l <- 2 * (f0$loglik - f1$loglik)
out("Sayisal ayrinti: baslangic S1(IT_j) = 1 - j/(J+1) (J = %d), lambda = 1, beta = 0; baglantisiz EM %d iterasyon (%s), log-baglanti EM %d iterasyon (%s); beta = 0 altinda logL = %.4f",
    f1$st$J, f0$iter, f0$route, f1$iter, f1$route, fr$loglik)
out("Sayisal ayrinti: log-baglanti altinda lambda_30 = exp(10 beta) = %.4f, lambda_40 = exp(20 beta) = %.4f", exp(10 * f1$beta), exp(20 * f1$beta))
out("Sayisal ayrinti: degerlendirme sayisi %s; ret sayisi %s (20, 30, 40 C)",
    paste(tapply(d$N, d$env, sum), collapse = ", "), paste(tapply(d$r, d$env, sum), collapse = ", "))
out("Sicaklik etkisi (beta = 0) LR = %.4f, sd 1, p = %.6f", lr_t, pchisq(lr_t, 1, lower.tail = FALSE))
out("Log-dogrusal baglanti yeterliligi (baglantisiz modele karsi) LR = %.4f, sd 1, p = %.4f", lr_l, pchisq(lr_l, 1, lower.tail = FALSE))
g0 <- gof_ph(f0); g1 <- gof_ph(f1)
out("Pearson: baglantisiz %.4f (sd %d, p = %.4f); log-baglanti %.4f (sd %d, p = %.4f)", g0$stat, g0$df, g0$p, g1$stat, g1$df, g1$p)
blk(format(s1, digits = 4))
hdr("Log-baglanti modelinden S_k(t) ve 4 C tahmini")
pr <- lapply(c(20, 30, 40, 4), function(x0) predict_ph(f1, x0))
names(pr) <- c("20", "30", "40", "4")
blk(data.frame(day = f1$st$grid, S20 = round(pr[["20"]]$S, 4), S30 = round(pr[["30"]]$S, 4),
               S40 = round(pr[["40"]]$S, 4), S4 = round(pr[["4"]]$S, 4),
               S4_lo = round(pr[["4"]]$lo, 4), S4_hi = round(pr[["4"]]$hi, 4)))
for (x0 in c(20, 30, 40)) {
  Sx <- pr[[as.character(x0)]]$S; out("%d C: S <= 0.5 ilk gun = %s", x0, f1$st$grid[which(Sx <= 0.5)[1]])
}

# ------------------------------------------------------------ LaTeX tables ---
hdr("LaTeX tablolari")
grid <- f1$st$grid
cellv <- function(k, v) sapply(grid, function(t) { i <- which(d$env == k & d$IT == t)
  if (length(i)) format(d[[v]][i]) else "--" })
tab_data <- data.frame(Day = grid)
for (k in 1:3) { tab_data[[sprintf("$N$ (%d \\textdegree C)", temp[k])]] <- cellv(k, "N")
                 tab_data[[sprintf("Unacc. (%d \\textdegree C)", temp[k])]] <- cellv(k, "r") }
out("%s\n", latex_table(tab_data,
  "Staggered acceptability data of green tea stored at three temperatures: number of evaluations and number of unacceptable responses on each test day.",
  "tab:gtdata", align = "rrrrrrr"))
rows <- do.call(rbind, lapply(1:3, function(k) {
  x <- np[[k]]; x <- x[x$day >= 20, ]
  Sk <- pr[[as.character(temp[k])]]; j <- match(x$day, Sk$IT)
  data.frame(`Temp.` = sprintf("%d \\textdegree C", temp[k]), Day = x$day, Observed = f4(x$raw), NPMLE = f4(x$S),
             `PH model` = f4(Sk$S[j]),
             `95\\% CI (PH)` = ifelse(is.finite(Sk$lo[j]), sprintf("(%s, %s)", f4(Sk$lo[j]), f4(Sk$hi[j])), "--"),
             check.names = FALSE)
}))
out("%s\n", latex_table(rows,
  "Green tea: estimated probability of acceptance from day 20 onwards. Observed proportion, non-parametric MLE for each temperature, and the proportional hazards model with log-linear link with its 95\\% confidence interval.",
  "tab:gtest", align = "lrrrrr"))
# Table in the layout of Tables 13 and 14 of Zhu and Balakrishnan (2022):
# S1(IT_j) and the hazard ratios under free hazard ratios and under the log-linear link
ss0 <- summary_ph(f0)$tab; ss1 <- summary_ph(f1)$tab; J <- f1$st$J
cS <- function(tab, j) if (is.finite(tab$se[j])) sprintf("(%s, %s)", f4(tab$logit_lo[j]), f4(tab$logit_hi[j])) else "--"
jj <- which(grid >= 20)
rows_ph <- data.frame(Parameter = sprintf("$S_1(%d)$", grid[jj]),
  `MLE (free)` = f4(ss0$mle[jj]), `SE (free)` = f4(ss0$se[jj]), `95\\% CI (free)` = sapply(jj, function(j) cS(ss0, j)),
  `MLE (link)` = f4(ss1$mle[jj]), `SE (link)` = f4(ss1$se[jj]), `95\\% CI (link)` = sapply(jj, function(j) cS(ss1, j)),
  check.names = FALSE)
lamrow <- function(k, mult) {
  l0 <- f0$lam[k]; s0l <- ss0$se[J + k - 1]
  l1 <- exp(mult * be); s1l <- mult * l1 * seb
  data.frame(Parameter = sprintf("$\\lambda_%d$", k), `MLE (free)` = f4(l0), `SE (free)` = f4(s0l),
    `95\\% CI (free)` = sprintf("(%s, %s)", f4(exp(log(l0) - z975 * s0l / l0)), f4(exp(log(l0) + z975 * s0l / l0))),
    `MLE (link)` = f4(l1), `SE (link)` = f4(s1l),
    `95\\% CI (link)` = sprintf("(%s, %s)", f4(exp(mult * (be - z975 * seb))), f4(exp(mult * (be + z975 * seb)))), check.names = FALSE)
}
rows_ph <- rbind(rows_ph, lamrow(2, 10), lamrow(3, 20),
  data.frame(Parameter = "$\\beta$", `MLE (free)` = "--", `SE (free)` = "--", `95\\% CI (free)` = "--",
    `MLE (link)` = f4(be), `SE (link)` = f4(seb), `95\\% CI (link)` = sprintf("(%s, %s)", f4(be - z975 * seb), f4(be + z975 * seb)),
    check.names = FALSE))
out("%s\n", latex_table(rows_ph,
  "Green tea: MLEs, standard errors (SE) and 95\\% confidence intervals under the proportional hazards model with free hazard ratios and with the log-linear link (logit transformation for $S_1$, log transformation for $\\lambda_2$ and $\\lambda_3$). Under the link, $\\lambda_2 = \\exp(10\\beta)$ and $\\lambda_3 = \\exp(20\\beta)$.",
  "tab:gtph", align = "lrrrrrr"))

# ----------------------------------------------------------------- figures ---
col_t <- c("#2a78d6", "#1baf7a", "#eb6834"); ink <- "#3d3d3a"
step_line <- function(t, S, tmax, ...) {   # right-continuous step from (0, 1) to tmax
  xs <- c(0, rep(t, each = 2), tmax); ys <- c(rep(c(1, S), each = 2))
  lines(xs, ys, ...)
}
fig <- function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 0.8))
  plot(NA, xlim = c(0, 72), ylim = c(0, 1), xlab = "Storage time (days)",
       ylab = "Probability of acceptance", las = 1)
  title(main = "(a) Non-parametric MLE", font.main = 1, cex.main = 1)
  for (k in 1:3) {
    x <- np[[k]]
    step_line(x$day, x$S, max(x$day) + 2, col = col_t[k], lwd = 2, lty = k)
    points(x$day, x$raw, pch = c(21, 22, 24)[k], bg = "white", col = col_t[k], cex = 1.1, lwd = 1.5)
  }
  legend("bottomleft", sprintf("%d \u00b0C", temp), lty = 1:3, lwd = 2, pch = c(21, 22, 24),
         pt.bg = "white", col = col_t, bty = "n", cex = 0.85)
  plot(NA, xlim = c(0, 72), ylim = c(0, 1), xlab = "Storage time (days)",
       ylab = "Probability of acceptance", las = 1)
  title(main = "(b) Proportional hazards model, log-linear link", font.main = 1, cex.main = 1)
  for (k in 1:3) {
    Sk <- pr[[as.character(temp[k])]]; keep <- Sk$IT <= max(days[[k]])
    step_line(Sk$IT[keep], Sk$S[keep], max(days[[k]]) + 2, col = col_t[k], lwd = 2, lty = k)
    s <- d[d$env == k, ]
    points(s$IT, 1 - s$r / s$N, pch = c(21, 22, 24)[k], bg = "white", col = col_t[k], cex = 1.1, lwd = 1.5)
  }
  S4 <- pr[["4"]]
  step_line(S4$IT, S4$S, 72, col = ink, lwd = 2, lty = 4)
  legend("bottomleft", c(sprintf("%d \u00b0C", temp), "4 \u00b0C (predicted)"), lty = c(1, 2, 3, 4), lwd = 2,
         pch = c(21, 22, 24, NA), pt.bg = "white", col = c(col_t, ink), bty = "n", cex = 0.85)
}
save_figure("greentea_fig4", fig, width = 9, height = 4.3)
writeLines(rep_lines, file.path("results", "run_greentea_report.txt"))
cat("run_greentea.R finished\n")
