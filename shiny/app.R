# =============================================================================
# app.R  --  One-shot shelf-life analysis of staggered sensory data
# Shiny application accompanying the P96_FOOD manuscript. It uses the same
# core (core_oneshot.R) as the paper, so its numbers agree with the tables.
# Run locally: open this file in RStudio and press "Run App".
# =============================================================================
library(shiny)

for (p in c("core_oneshot.R", file.path("..", "code", "core_oneshot.R"))) {
  if (file.exists(p)) { source(p, local = FALSE); break }
}
if (!exists("fit_ph")) stop("core_oneshot.R not found next to app.R or in ../code/")

# ------------------------------------------------------------ example data ---
ex_pezik <- data.frame(
  group = c(rep("Blanched 7 min", 10), rep("Blanched 9 min", 9)),
  x = c(rep(7, 10), rep(9, 9)),
  day = c(1, 20, 40, 60, 80, 100, 120, 140, 160, 170, 1, 20, 40, 60, 80, 100, 120, 140, 150),
  N = c(3:11, 19, 3:10, 18),
  r = c(0, 0, 0, 2, 1, 1, 2, 3, 7, 12, 0, 0, 1, 2, 1, 2, 4, 7, 13))
ex_tea <- data.frame(
  group = c(rep("20 C", 14), rep("30 C", 11), rep("40 C", 9)),
  x = c(rep(20, 14), rep(30, 11), rep(40, 9)),
  day = c(seq(5, 70, 5), seq(5, 55, 5), seq(5, 45, 5)),
  N = c(2:10, 10, 10, 10, 10, 10, 2:8, 9, 10, 10, 10, 2:5, 6, 9, 10, 10, 10),
  r = c(0,0,0,0,0,0,0,0,0,2,0,1,3,6, 0,0,0,0,0,0,0,2,2,5,8, 0,0,0,0,2,2,2,3,8))

pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#7a5cc9", "#5f5e5a")
ink <- "#3d3d3a"

step_line <- function(t, S, tmax, ...) {
  lines(c(0, rep(t, each = 2), tmax), rep(c(1, S), each = 2), ...)
}
fmt <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 4))

includeMarkdownText <- function() {
  div(
    p("Each evaluation in a staggered sensory design shows only whether the product had already become unacceptable on the day of the test.",
      "Such responses are one-shot device test data: an unacceptable response is a left-censored observation and an acceptable one is right-censored."),
    p("The application estimates the probability of acceptance on every test day by the non-parametric maximum likelihood estimator",
      "(EM algorithm), compares groups through the semi-parametric proportional hazards model, optionally with a log-linear link in",
      "a numeric covariate such as storage temperature, and contrasts the results with the Weibull hazard method."),
    p("Confidence intervals: logit scale for single-group estimates, log scale for hazard ratios, log(-log) scale for predictions.",
      "Groups must share test days for the proportional hazards model to be identifiable."),
    p(class = "note", "Code: Coskun Kus, Selcuk University. All computations use base R; the code is the one used for the tables of the accompanying paper.")
  )
}

# ---------------------------------------------------------------------- UI ---
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Segoe UI', Arial, sans-serif; }
    .well { background: #f6f6f4; }
    h2 { font-size: 22px; }
    .note { color: #555; font-size: 13px; }
    pre { background: #f6f6f4; border: none; }"))),
  titlePanel("Shelf life from staggered sensory tests: one-shot device analysis"),
  sidebarLayout(
    sidebarPanel(width = 3,
      radioButtons("src", "Data", c("Pezik pickles (Keklik et al., 2017)" = "pezik",
                                    "Green tea (Xiao et al., 2022)" = "tea",
                                    "Upload a CSV file" = "upload")),
      conditionalPanel("input.src == 'upload'",
        fileInput("file", NULL, accept = c(".csv", ".txt")),
        p(class = "note", "Columns: group, x (numeric, e.g. temperature; optional), day, N, r.",
          "N = number of evaluations on that day, r = number found unacceptable.")),
      uiOutput("ref_ui"),
      radioButtons("link", "Proportional hazards model",
                   c("Free hazard ratio for each group" = "none",
                     "Log-linear link in x" = "log")),
      conditionalPanel("input.link == 'log'",
        numericInput("x0", "Predict at x =", value = 4)),
      hr(),
      p(class = "note", "Methods: Zhu X, Balakrishnan N (2022), Reliability Engineering and System Safety 221:108319.")
    ),
    mainPanel(width = 9,
      tabsetPanel(
        tabPanel("Data", br(), tableOutput("data_tab"), uiOutput("data_msg")),
        tabPanel("Non-parametric MLE", br(),
                 plotOutput("np_plot", height = "420px"),
                 tableOutput("np_tab"), downloadButton("np_dl", "Download table")),
        tabPanel("Proportional hazards", br(),
                 verbatimTextOutput("ph_txt"),
                 plotOutput("ph_plot", height = "420px"),
                 tableOutput("ph_tab"), downloadButton("ph_dl", "Download table")),
        tabPanel("Prediction", br(), uiOutput("pred_msg"),
                 plotOutput("pred_plot", height = "380px"), tableOutput("pred_tab")),
        tabPanel("Weibull hazard method", br(),
                 p(class = "note", "Modified version (Cardelli and Labuza, 2001): only the unacceptable responses are ranked.",
                   "Original version (Gacula, 1975): all samples are ranked. Shelf life = time at which the cumulative hazard reaches 69.3%."),
                 checkboxInput("trunc", "Use only the points with cumulative hazard up to 100% (Cardelli and Labuza rule)", TRUE),
                 tableOutput("wh_tab"), uiOutput("wh_grp_ui"), plotOutput("wh_plot", height = "400px")),
        tabPanel("About", br(), includeMarkdownText())
      ))
  )
)

# ------------------------------------------------------------------ server ---
server <- function(input, output, session) {

  raw <- reactive({
    if (input$src == "pezik") return(ex_pezik)
    if (input$src == "tea") return(ex_tea)
    req(input$file)
    df <- read.csv(input$file$datapath, stringsAsFactors = FALSE)
    names(df) <- tolower(trimws(names(df)))
    validate(need(all(c("group", "day", "n", "r") %in% names(df)),
                  "The file needs the columns group, day, N and r (x is optional)."))
    names(df)[names(df) == "n"] <- "N"
    if (!"x" %in% names(df)) df$x <- NA_real_
    df <- df[, c("group", "x", "day", "N", "r")]
    validate(need(all(df$r >= 0 & df$r <= df$N), "Every r must lie between 0 and N."))
    df
  })

  output$ref_ui <- renderUI({
    g <- unique(raw()$group)
    selectInput("ref", "Reference group", g, selected = g[1])
  })

  dat <- reactive({
    df <- raw(); req(input$ref, input$ref %in% df$group)
    lev <- c(input$ref, setdiff(unique(df$group), input$ref))
    env <- match(df$group, lev)
    xx <- if (all(is.finite(df$x))) df$x else env
    d <- os_data(env = env, x = xx, IT = df$day, N = df$N, r = df$r)
    list(d = d, lev = lev, has_x = all(is.finite(df$x)))
  })

  output$data_tab <- renderTable(raw(), digits = 0)
  output$data_msg <- renderUI({
    d <- dat()$d; K <- max(d$env)
    shared <- sapply(1:K, function(k) any(d$IT[d$env == k] %in% d$IT[d$env != k]))
    if (K > 1 && !all(shared))
      p(style = "color:#b3401a", "Warning: some groups share no test day with the others; the proportional hazards model is not identifiable for them.")
  })

  # ----------------------------------------------------- non-parametric ---
  np <- reactive({
    d <- dat()$d; lev <- dat()$lev
    do.call(rbind, lapply(seq_along(lev), function(k) {
      s <- d[d$env == k, ]; p <- pava_np(s$N, s$r); ci <- ci_S(p, se_np(p, s$N))
      data.frame(group = lev[k], day = s$IT, N = s$N, r = s$r, observed = 1 - s$r / s$N,
                 NPMLE = p, lower = ifelse(is.finite(ci$logit[, 1]), ci$logit[, 1], NA),
                 upper = ifelse(is.finite(ci$logit[, 2]), ci$logit[, 2], NA))
    }))
  })
  output$np_tab <- renderTable({
    x <- np(); x[, 5:8] <- lapply(x[, 5:8], fmt); x }, digits = 0)
  output$np_dl <- downloadHandler("npmle.csv", function(f) write.csv(np(), f, row.names = FALSE))
  output$np_plot <- renderPlot({
    x <- np(); lev <- dat()$lev; tmax <- max(x$day) * 1.03
    par(mar = c(4.2, 4.2, 1, 1))
    plot(NA, xlim = c(0, tmax), ylim = c(0, 1), xlab = "Storage time (days)",
         ylab = "Probability of acceptance", las = 1)
    for (k in seq_along(lev)) {
      s <- x[x$group == lev[k], ]
      step_line(s$day, s$NPMLE, max(s$day) + 0.02 * tmax, col = pal[k], lwd = 2, lty = k)
      ok <- is.finite(s$lower); off <- (k - (length(lev) + 1) / 2) * 0.008 * tmax
      segments(s$day[ok] + off, s$lower[ok], s$day[ok] + off, s$upper[ok], col = pal[k], lwd = 1.5)
      points(s$day + off, s$observed, pch = 20 + k %% 5, bg = "white", col = pal[k], cex = 1.1)
    }
    legend("bottomleft", lev, col = pal[seq_along(lev)], lty = seq_along(lev), lwd = 2, bty = "n")
  })

  # ---------------------------------------------- proportional hazards ---
  ph <- reactive({
    d <- dat()$d; K <- max(d$env)
    validate(need(K >= 2, "At least two groups are needed."))
    link <- input$link
    if (link == "log") validate(need(dat()$has_x && length(unique(d$x)) == K,
                                     "The log-linear link needs a numeric x that differs between groups."))
    f <- fit_ph(d, link, it_em = 100)
    fr <- fit_ph(transform(d, env = 1L), "none", it_em = 100)
    out <- list(f = f, s = summary_ph(f), gof = gof_ph(f),
                lr = 2 * (f$loglik - fr$loglik), df_lr = if (link == "log") 1 else K - 1)
    if (link == "log") {
      f0 <- fit_ph(d, "none", it_em = 100)
      out$lr_link <- 2 * (f0$loglik - f$loglik); out$df_link <- K - 2
    }
    out
  })
  output$ph_txt <- renderPrint({
    o <- ph(); f <- o$f; tb <- o$s$tab; lev <- dat()$lev; J <- f$st$J
    if (f$link == "none") {
      cat("Hazard ratios relative to", lev[1], "(95% CI on the log scale)\n")
      for (k in 2:length(lev)) {
        l <- f$lam[k]; se <- tb$se[J + k - 1]
        cat(sprintf("  %-20s %.4f  (%.4f, %.4f)\n", lev[k], l,
                    exp(log(l) - z975 * se / l), exp(log(l) + z975 * se / l)))
      }
    } else {
      b <- f$beta; se <- tb$se[J + 1]
      cat(sprintf("beta = %.4f  (SE %.4f), 95%% CI (%.4f, %.4f)\n", b, se, b - z975 * se, b + z975 * se))
      cat(sprintf("Hazard ratio for a 10-unit increase in x: %.4f  (%.4f, %.4f)\n",
                  exp(10 * b), exp(10 * (b - z975 * se)), exp(10 * (b + z975 * se))))
      if (!is.null(o$lr_link) && o$df_link > 0)
        cat(sprintf("Link adequacy (LR against free hazard ratios): %.4f, df %d, p = %.4f\n",
                    o$lr_link, o$df_link, pchisq(o$lr_link, o$df_link, lower.tail = FALSE)))
    }
    cat(sprintf("No group effect (LR test): %.4f, df %d, p = %.4f\n", o$lr, o$df_lr,
                pchisq(o$lr, o$df_lr, lower.tail = FALSE)))
    cat(sprintf("Goodness of fit (Pearson): %.4f, df %d, p = %.4f\n", o$gof$stat, o$gof$df, o$gof$p))
  })
  ph_fit_tab <- reactive({
    o <- ph(); f <- o$f; d <- dat()$d; lev <- dat()$lev
    S <- f$S1[f$st$g]^f$lam[d$env]
    data.frame(group = lev[d$env], day = d$IT, N = d$N, r = d$r,
               observed = 1 - d$r / d$N, fitted = S)
  })
  output$ph_tab <- renderTable({ x <- ph_fit_tab(); x[, 5:6] <- lapply(x[, 5:6], fmt); x }, digits = 0)
  output$ph_dl <- downloadHandler("ph_model.csv", function(f) write.csv(ph_fit_tab(), f, row.names = FALSE))
  output$ph_plot <- renderPlot({
    o <- ph(); f <- o$f; d <- dat()$d; lev <- dat()$lev; tmax <- max(d$IT) * 1.03
    par(mar = c(4.2, 4.2, 1, 1))
    plot(NA, xlim = c(0, tmax), ylim = c(0, 1), xlab = "Storage time (days)",
         ylab = "Probability of acceptance", las = 1)
    for (k in seq_along(lev)) {
      keep <- f$st$grid <= max(d$IT[d$env == k])
      step_line(f$st$grid[keep], (f$S1^f$lam[k])[keep], max(d$IT[d$env == k]) + 0.02 * tmax,
                col = pal[k], lwd = 2, lty = k)
      s <- d[d$env == k, ]
      points(s$IT, 1 - s$r / s$N, pch = 20 + k %% 5, bg = "white", col = pal[k], cex = 1.1)
    }
    legend("bottomleft", lev, col = pal[seq_along(lev)], lty = seq_along(lev), lwd = 2, bty = "n")
  })

  # ----------------------------------------------------------- prediction ---
  pred <- reactive({
    req(input$link == "log"); o <- ph(); req(is.finite(input$x0))
    predict_ph(o$f, input$x0)
  })
  output$pred_msg <- renderUI({
    if (input$link != "log")
      p("Choose the log-linear link in the sidebar to predict at an untested value of x (for example, a storage temperature).")
    else p(class = "note", sprintf(
      "Probability of acceptance at x = %g on the test days of the data, with 95%% CI on the log(-log) scale. Values far outside the tested range of x are extrapolations.",
      input$x0))
  })
  output$pred_tab <- renderTable({
    x <- pred(); data.frame(day = x$IT, S = fmt(x$S), lower = fmt(x$lo), upper = fmt(x$hi)) }, digits = 0)
  output$pred_plot <- renderPlot({
    x <- pred(); tmax <- max(x$IT) * 1.03
    par(mar = c(4.2, 4.2, 1, 1))
    plot(NA, xlim = c(0, tmax), ylim = c(0, 1), xlab = "Storage time (days)",
         ylab = "Probability of acceptance", las = 1)
    ok <- is.finite(x$lo)
    segments(x$IT[ok], x$lo[ok], x$IT[ok], x$hi[ok], col = ink)
    step_line(x$IT, x$S, max(x$IT) + 0.02 * tmax, col = ink, lwd = 2)
    legend("bottomleft", sprintf("Predicted at x = %g", input$x0), col = ink, lwd = 2, bty = "n")
  })

  # ------------------------------------------------------- Weibull hazard ---
  wh <- reactive({
    d <- dat()$d; lev <- dat()$lev
    do.call(rbind, lapply(seq_along(lev), function(k) {
      s <- d[d$env == k, ]
      chm <- if (isTRUE(input$trunc)) 100 else Inf
      gm <- gacula_fit(s$IT, s$r, CH_max = chm); go <- gacula_fit(s$IT, s$r, CH_max = chm, N = s$N)
      row <- function(g, v) if (is.null(g)) data.frame(group = lev[k], version = v, beta = NA, t50 = NA, lower = NA, upper = NA)
        else data.frame(group = lev[k], version = v, beta = g$beta, t50 = g$t50, lower = g$lo, upper = g$hi)
      rbind(row(gm, "Modified (failures only)"), row(go, "Original (all samples)"))
    }))
  })
  output$wh_tab <- renderTable({
    x <- wh(); x$beta <- formatC(x$beta, format = "f", digits = 2)
    x[, 4:6] <- lapply(x[, 4:6], function(v) formatC(v, format = "f", digits = 1))
    names(x) <- c("group", "version", "shape", "shelf life (days)", "lower 95%", "upper 95%"); x })
  output$wh_grp_ui <- renderUI(selectInput("wh_grp", "Group for the plot", dat()$lev))
  output$wh_plot <- renderPlot({
    req(input$wh_grp); d <- dat()$d; k <- match(input$wh_grp, dat()$lev); req(!is.na(k)); s <- d[d$env == k, ]
    p <- pava_np(s$N, s$r); tmax <- max(s$IT) * 1.05; chm <- if (isTRUE(input$trunc)) 100 else Inf
    gm <- gacula_fit(s$IT, s$r, CH_max = chm); go <- gacula_fit(s$IT, s$r, CH_max = chm, N = s$N)
    par(mar = c(4.2, 4.2, 1, 1))
    plot(NA, xlim = c(0, tmax), ylim = c(0, 1), xlab = "Storage time (days)",
         ylab = "Probability of acceptance", las = 1)
    if (!is.null(gm)) curve(gm$S(x), from = 0.5, to = tmax, add = TRUE, col = pal[2], lwd = 2, lty = 2)
    if (!is.null(go)) curve(go$S(x), from = 0.5, to = tmax, add = TRUE, col = pal[4], lwd = 2, lty = 3)
    step_line(s$IT, p, max(s$IT) + 0.02 * tmax, col = pal[1], lwd = 2)
    points(s$IT, 1 - s$r / s$N, pch = 21, bg = "white", col = ink, cex = 1.1)
    abline(h = 0.5, col = "grey60", lty = 3)
    legend("bottomleft", c("Non-parametric MLE", "Weibull hazard, modified", "Weibull hazard, original", "Observed"),
           col = c(pal[1], pal[2], pal[4], ink), lty = c(1, 2, 3, NA), lwd = c(2, 2, 2, NA),
           pch = c(NA, NA, NA, 21), pt.bg = "white", bty = "n")
  })
}

shinyApp(ui, server)
