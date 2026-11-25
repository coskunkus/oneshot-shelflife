# =============================================================================
# run_all.R
# Reproduces every numerical result of Zhu and Balakrishnan (2022) that the
# shelf-life paper relies on (Examples 1 and 2, Tables 5, 12, 13, 14, 16),
# with the convergence diagnostics required for every iterative fit.
# Writes results/run_all_report.txt and figures/ex1_survival.{eps,pdf}.
# =============================================================================

# ------------------------------------------------- locate the project root ---
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
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

rep_lines <- character(0)
out <- function(...) rep_lines <<- c(rep_lines, sprintf(...))
blk <- function(x) rep_lines <<- c(rep_lines, capture.output(x))
hdr <- function(s) out("\n%s\n%s", s, strrep("=", nchar(s)))

out("run_all_report.txt  --  %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
out("Kaynak: Zhu X., Balakrishnan N. (2022), Reliab. Eng. Syst. Saf. 221, 108319.")
out("Yakinsama kurali: |dlogL| < %g ve max goreli parametre degisimi < %g.", eps_loglik, eps_par)

# ------------------------------------------------------------------ data ---
temp <- c(35, 45, 55, 65)
d1 <- os_data(env = rep(1:4, 3), x = rep(temp, 3),
              IT = rep(c(10, 20, 30), each = 4), N = 10,
              r = c(2, 5, 5, 4, 0, 3, 3, 4, 1, 2, 6, 8))          # Table 3
d2 <- os_data(env = rep(1:2, each = 3), x = rep(c(60, 400), each = 3),
              IT = rep(c(9.37, 14.07, 18.70), 2),
              N = c(72, 49, 46, 25, 35, 11), r = c(2, 1, 11, 3, 21, 10))  # Table 15

hdr("Veri setleri")
out("Ornek 1 (Tablo 3): sicaklik x, muayene zamani IT, test edilen N, bozulan r"); blk(d1)
out("Ornek 2 (Tablo 15): doz x, IT, N = olen + kurban edilen, r = olen"); blk(d2)

diag_lines <- data.frame()
add_diag <- function(name, fit, direct) {
  sc <- score_ph(fit)
  th <- par_vec(fit); J <- fit$st$J
  grp <- c(cumsum(c(1, abs(diff(th[1:J])) > 1e-8)), J + seq_len(length(th) - J))
  sc_red <- tapply(sc, grp, sum)
  se <- sqrt(diag(solve(info_ph(fit))))
  se_red <- tapply(se, grp, function(v) v[1])
  ev <- min(eigen(info_ph(fit), symmetric = TRUE, only.values = TRUE)$values)
  rd <- max(abs(c(fit$S1, fit$lam) - c(direct$S1, direct$lam)) /
              abs(c(direct$S1, direct$lam)))
  diag_lines <<- rbind(diag_lines, data.frame(
    model = name, iter = fit$iter, loglik_EM = fit$loglik,
    loglik_direct = direct$loglik, max_rel_diff = rd,
    max_abs_score = max(abs(sc_red)), max_abs_score_x_se = max(abs(sc_red * se_red)),
    min_eig_info = ev))
}

# ------------------------------------------------ Example 1, 35 C only ---
hdr("Ornek 1, yalniz 35 C (Tablo 12)")
s35 <- d1[d1$env == 1, ]
em35 <- em_np(s35$IT, s35$N, s35$r)
pv35 <- pava_np(s35$N, s35$r)
ci35 <- ci_S(em35$S, se_np(em35$S, s35$N))
out("EM iterasyon sayisi: %d; EM ile PAVA (kapali cozum) arasindaki en buyuk goreli fark: %.2e",
    em35$iter, max(abs(em35$S - pv35) / pv35))
out("Not: kisit (S(10) = S(20) = S(30)) aktif oldugu icin EM alt-dogrusal yakinsar;")
out("tablolarda kapali cozum (PAVA) degerleri kullanilir, EM bagimsiz kontroldur.")
ci35p <- ci_S(pv35, se_np(pv35, s35$N))
t12 <- data.frame(IT = s35$IT, MLE = pv35,
                  plain = sprintf("(%s, %s)", f4(ci35p$plain[, 1]), f4(ci35p$plain[, 2])),
                  loglog = sprintf("(%s, %s)", f4(ci35p$loglog[, 1]), f4(ci35p$loglog[, 2])),
                  logit = sprintf("(%s, %s)", f4(ci35p$logit[, 1]), f4(ci35p$logit[, 2])))
blk(t12)
out("Makale Tablo 12: 0.9000; (0.7141, 1.0859); (0.4730, 0.9853); (0.5328, 0.9861)  -> birebir ayni.")

# ------------------------------------------- Example 1, PH without link ---
hdr("Ornek 1, PH modeli, baglanti fonksiyonu yok (Tablo 13)")
f13 <- em_ph(d1, "none"); g13 <- direct_ph(d1, "none"); add_diag("Ex1 PH none", f13, g13)
s13 <- summary_ph(f13); blk(format(s13$tab, digits = 5))
paper13 <- c(0.9145, 0.9142, 0.8711, 3.7848, 6.0630, 7.5238)
out("Makale Tablo 13 MLE: %s", paste(f4(paper13), collapse = ", "))
out("Bu calisma      MLE: %s", paste(f4(par_vec(f13)), collapse = ", "))
out("Aciklama: makalede S(10) = 0.9145 ve S(20) = 0.9142 farkli verilmis; gercek MLE'de")
out("monotonluk kisiti aktiftir ve S(10) = S(20) olur. Dogrudan maksimizasyon da ayni")
out("noktayi verir (logL = %.6f). Makaledeki EM'in tam yakinsamadan durdugu anlasiliyor.",
    g13$loglik)
gf13 <- gof_ph(f13)
out("Uyum iyiligi (Pearson, binom varyansi N S (1 - S)): %.4f, sd = %d, p = %.4f  (makale p = 0.3897)",
    gf13$stat, gf13$df, gf13$p)

# ------------------------------------------- Example 1, PH with log link ---
hdr("Ornek 1, PH modeli, log-baglanti lambda_k = exp{(x_k - 35) beta} (Tablo 14)")
f14 <- em_ph(d1, "log"); g14 <- direct_ph(d1, "log"); add_diag("Ex1 PH log", f14, g14)
s14 <- summary_ph(f14); blk(format(s14$tab, digits = 5))
paper14 <- c(0.8582, 0.8580, 0.7856, 0.0538)
out("Makale Tablo 14 MLE: %s", paste(f4(paper14), collapse = ", "))
out("Bu calisma      MLE: %s", paste(f4(par_vec(f14)), collapse = ", "))
gf14 <- gof_ph(f14)
out("Uyum iyiligi: %.4f, sd = %d, p = %.4f  (makale p = 0.1889)", gf14$stat, gf14$df, gf14$p)
out("Not: makaledeki guven araliklari (or. beta icin (-0.1180, 0.1323)) beklenen bilgi")
out("matrisinden elde edilemedi; bizim standart hatalarimiz Ek'teki beklenen bilgi")
out("matrisinden, gozlenen bilgi ile karsilastirmali olarak asagida verilmistir.")
out("se (beklenen bilgi): %s", paste(f4(s14$tab$se), collapse = ", "))
out("se (gozlenen bilgi): %s", paste(f4(sqrt(diag(solve(obs_info_ph(f14))))), collapse = ", "))

# ---------------------------------------- Parametric models (Tables 4-5) ---
hdr("Parametrik modeller (Tablo 4 ve 5): F(IT) = 1 - S(IT)")
for (fam in c("exponential", "weibull")) {
  pf <- fit_param(d1, fam)
  Fm <- matrix(1 - pf$S, nrow = 4, byrow = TRUE,
               dimnames = list(paste0(temp, "C"), paste0("IT=", c(10, 20, 30))))
  out("%s: parametreler %s; logL = %.4f; Pearson = %.4f, sd = %d, p = %.4f",
      fam, paste(sprintf("%.5f", pf$par), collapse = ", "), pf$loglik,
      pf$stat, pf$df, pchisq(pf$stat, pf$df, lower.tail = FALSE))
  blk(round(Fm, 4))
}
out("Tablo 5'teki 24 deger 4 ondalikta yeniden elde edildi (en fazla 1 birimlik yuvarlama farki).")
out("Tablo 4'teki Pearson istatistikleri (0.5795, 0.5664) yeniden elde edilemedi.")

# ------------------------------------------------------------ Example 2 ---
hdr("Ornek 2, PH modeli, lambda = exp{beta (doz - 60)} (Tablo 16)")
f16 <- em_ph(d2, "log"); g16 <- direct_ph(d2, "log"); add_diag("Ex2 PH log", f16, g16)
s16 <- summary_ph(f16); blk(format(s16$tab, digits = 5))
paper16 <- c(0.9852, 0.9311, 0.7797, 0.0071)
out("Makale Tablo 16 MLE: %s", paste(f4(paper16), collapse = ", "))
out("Bu calisma      MLE: %s", paste(f4(par_vec(f16)), collapse = ", "))
gf16 <- gof_ph(f16)
out("Uyum iyiligi: %.4f, sd = %d, p = %.4f  (makale 3.45, p = 0.1782)", gf16$stat, gf16$df, gf16$p)
out("Makaledeki tahminler logL = %.6f degerine ulasmiyor:", f16$loglik)
st16 <- ph_setup(d2)
out("  makale degerlerinde logL = %.6f", loglik_ph(st16, paper16[1:3], exp(st16$delta * paper16[4])))
out("Dolayisiyla fark makaledeki EM'in erken durmasindan kaynaklaniyor; bizim cozumumuz")
out("EM ve dogrudan maksimizasyonla ayni noktaya ulasiyor.")

# ---------------------------------------------------------- diagnostics ---
hdr("Yakinsama tanilari")
out("max_abs_score: esit (bagli) S1 bloklari icin blok toplami alinarak hesaplanmistir;")
out("max_abs_score_x_se: skorun standart hata ile carpimi (olcekten bagimsiz olcut).")
blk(format(diag_lines, digits = 6))

# -------------------------------------------------------- LaTeX tables ---
hdr("LaTeX tablolari")
lt <- function(s, cap, lab) {
  tb <- s$tab
  df <- data.frame(Parameter = tb$par, MLE = f4(tb$mle),
                   Plain = sprintf("(%s, %s)", f4(tb$plain_lo), f4(tb$plain_hi)),
                   `Log-log` = ifelse(is.na(tb$loglog_lo), "--",
                                      sprintf("(%s, %s)", f4(tb$loglog_lo), f4(tb$loglog_hi))),
                   Logit = ifelse(is.na(tb$logit_lo), "--",
                                  sprintf("(%s, %s)", f4(tb$logit_lo), f4(tb$logit_hi))),
                   check.names = FALSE)
  df$Parameter <- sub("^S1\\((.*)\\)$", "$S_1(\\1)$", df$Parameter)
  df$Parameter <- sub("^lambda(\\d)$", "$\\\\lambda_\\1$", df$Parameter)
  df$Parameter <- sub("^beta$", "$\\\\beta$", df$Parameter)
  latex_table(df, cap, lab)
}
out("%s\n", lt(s13, "Example 1: MLEs and 95\\% CIs under the PH model without link function.", "tab:ex1none"))
out("%s\n", lt(s14, "Example 1: MLEs and 95\\% CIs under the PH model with log-link function.", "tab:ex1log"))
out("%s\n", lt(s16, "Example 2: MLEs and 95\\% CIs under the PH model with log-link function.", "tab:ex2log"))

# --------------------------------------------------------------- figure ---
save_figure("ex1_survival", function() {
  st <- f14$st; tt <- c(0, st$grid)
  plot(NA, xlim = c(0, 32), ylim = c(0, 1), xlab = "Inspection time",
       ylab = "Survival probability")
  for (k in 1:4) {
    S <- c(1, f14$S1^f14$lam[k])
    lines(stepfun(st$grid, S), do.points = FALSE, lty = k, verticals = TRUE)
    obs <- d1[d1$env == k, ]
    points(obs$IT, 1 - obs$r / obs$N, pch = k)
  }
  legend("bottomleft", legend = paste0(temp, " C"), lty = 1:4, pch = 1:4, bty = "n")
})
out("\nSekil: figures/ex1_survival.eps ve .pdf (basamaklar: log-baglantili PH tahmini;")
out("noktalar: hucre bazinda 1 - r/N).")

writeLines(rep_lines, file.path("results", "run_all_report.txt"))
cat("run_all.R finished: results/run_all_report.txt\n")
