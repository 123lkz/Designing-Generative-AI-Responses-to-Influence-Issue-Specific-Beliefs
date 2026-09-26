# Reproduce the manuscript analyses from the supplied participant-level totals
# and author-reviewed final response codes. No observations or codes are changed.
# Run: Rscript R_Analyse/run_analysis.R [project_root] [output_directory]
# Required packages: readxl, mvtnorm. See the project README.md.

options(stringsAsFactors = FALSE, warn = 1, digits = 12)
if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
}
set.seed(20260926)
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
default_root <- if (length(script_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)))
} else getwd()
root <- normalizePath(if (length(args)) args[1] else default_root, mustWork = TRUE)
out <- if (length(args) >= 2L) args[2] else file.path(root, "R_Analyse", "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
needed <- c("readxl", "mvtnorm")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
save_table <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")),
                                         row.names = FALSE, fileEncoding = "UTF-8", na = "")
checks <- list()
check <- function(label, ok) {
  checks[[length(checks) + 1L]] <<- data.frame(check = label, passed = isTRUE(ok))
  if (!isTRUE(ok)) stop("Input validation failed: ", label)
}
read_data <- function(experiment, name) {
  d <- read.csv(file.path(root, "Exp_Data", experiment, name), fileEncoding = "UTF-8-BOM")
  check(paste(name, "unique IDs and complete cells"), !anyDuplicated(d$id) && !anyNA(d))
  d
}
g1 <- c("control", "concept", "emotion", "behavior", "environment")
g2 <- c("control", "conceptual", "combined")
topics <- c("procrastination", "social_media", "multitasking")
covariates <- c("nfc_pre", "autonomy_sat_pre", "autonomy_frus_pre",
               "competence_sat_pre", "competence_frus_pre", "relatedness_sat_pre",
               "relatedness_frus_pre", "normative_influence_pre",
               "informational_influence_pre", "topic_involvement_pre",
               "info_processing_motivation_pre")
fit <- function(response, predictors, data) lm(reformulate(predictors, response), data = data,
                                               na.action = na.fail)
coef_table <- function(m) {
  a <- coef(summary(m)); ci <- confint(m)
  data.frame(term = rownames(a), estimate = a[, 1], se = a[, 2], t = a[, 3],
             p = a[, 4], ci_low = ci[, 1], ci_high = ci[, 2], row.names = NULL)
}
model_summary <- function(m) {
  s <- summary(m); f <- s$fstatistic
  data.frame(n = nobs(m), df_num = unname(f[2]), df_denom = unname(f[3]),
             F = unname(f[1]), p = pf(f[1], f[2], f[3], lower.tail = FALSE),
             R2 = s$r.squared, adjusted_R2 = s$adj.r.squared)
}
# HC3 sandwich covariance and finite-residual-df t/F inference match the paper.
hc3 <- function(m) {
  x <- model.matrix(m); bread <- solve(crossprod(x))
  u <- residuals(m) / (1 - hatvalues(m))
  bread %*% crossprod(x * u) %*% bread
}
wald <- function(m, L, V = vcov(m)) {
  L <- as.matrix(L); b <- drop(L %*% coef(m)); q <- nrow(L)
  stat <- drop(crossprod(b, solve(L %*% V %*% t(L), b))) / q
  data.frame(df_num = q, df_denom = df.residual(m), F = stat,
             p = pf(stat, q, df.residual(m), lower.tail = FALSE))
}
term_test <- function(m, term, V = vcov(m)) {
  labels <- attr(terms(m), "term.labels")
  k <- match(term, labels)
  if (is.na(k)) stop("Unknown model term: ", term)
  idx <- which(attr(model.matrix(m), "assign") == k)
  cbind(term = term, wald(m, diag(length(coef(m)))[idx, , drop = FALSE], V))
}
linear_estimate <- function(m, L, V = vcov(m)) {
  L <- as.matrix(L); b <- drop(L %*% coef(m))
  se <- sqrt(diag(L %*% V %*% t(L))); df <- df.residual(m)
  data.frame(estimate = b, se = se, df = df, ci_low = b - qt(.975, df) * se,
             ci_high = b + qt(.975, df) * se, t = b / se,
             p = 2 * pt(-abs(b / se), df))
}
# Two-sided single-step Dunnett adjustment over the four pathway-vs-control
# contrasts. Confidence intervals are simultaneous within this family.
dunnett <- function(m) {
  idx <- grep("^group", names(coef(m)))
  stopifnot(length(idx) == 4L)
  L <- diag(length(coef(m)))[idx, , drop = FALSE]
  z <- linear_estimate(m, L); v <- vcov(m)[idx, idx]; r <- cov2cor(v)
  algorithm <- mvtnorm::GenzBretz(maxpts = 100000, abseps = 1e-5)
  z$p_dunnett <- vapply(abs(z$t), function(tval) {
    val <- mvtnorm::pmvt(lower = rep(-tval, 4), upper = rep(tval, 4),
                         df = df.residual(m), corr = r, algorithm = algorithm,
                         seed = 20260926)
    max(0, min(1, 1 - as.numeric(val)))
  }, numeric(1))
  coverage <- function(cutoff) as.numeric(mvtnorm::pmvt(
    lower = rep(-cutoff, 4), upper = rep(cutoff, 4), df = df.residual(m),
    corr = r, algorithm = algorithm, seed = 20260926))
  critical <- uniroot(function(cutoff) coverage(cutoff) - .95,
                     c(qt(.975, df.residual(m)), qt(1 - .05 / 8, df.residual(m))),
                     tol = 1e-5)$root
  z$ci_low <- z$estimate - critical * z$se
  z$ci_high <- z$estimate + critical * z$se
  cbind(comparison = paste(sub("^group", "", names(coef(m))[idx]), "- control"), z)
}
describe <- function(d, variables, by = "group") {
  keys <- unique(d[by]); result <- list()
  for (i in seq_len(nrow(keys))) {
    keep <- rep(TRUE, nrow(d))
    for (k in by) keep <- keep & d[[k]] == keys[[k]][i]
    for (v in variables) {
      x <- d[[v]][keep]
      result[[length(result) + 1L]] <- cbind(keys[i, , drop = FALSE],
        data.frame(outcome = v, n = length(x), mean = mean(x), sd = sd(x),
                   se = sd(x) / sqrt(length(x)), min = min(x), max = max(x)))
    }
  }
  do.call(rbind, result)
}

message("Experiment 1: quantitative analyses")
exp1 <- list(
  PBS = read_data("Exp1", "PBS_exp1_for_R.csv"),
  CFI = read_data("Exp1", "CFI_exp1_for_R.csv"),
  BCIS = read_data("Exp1", "BCIS_exp1_for_R.csv"))
outcomes <- list(PBS = "pbs", CFI = c("cfi_alternatives", "cfi_control"),
                 BCIS = c("bcis_total", "bcis_reflection", "bcis_certainty"))
for (name in names(exp1)) {
  d <- exp1[[name]]
  check(paste(name, "160 participants and 32 per condition"),
        nrow(d) == 160L && setequal(d$group, g1) && all(table(d$group) == 32L))
  base <- exp1$PBS[match(d$id, exp1$PBS$id), c("group", covariates)]
  check(paste(name, "IDs, groups and covariates match PBS"),
        isTRUE(all.equal(unname(as.matrix(d[c("group", covariates)])), unname(as.matrix(base)))))
  d$group <- factor(d$group, levels = g1)
  contrasts(d$group) <- contr.treatment(length(g1))
  colnames(contrasts(d$group)) <- g1[-1]
  exp1[[name]] <- d
}
check("Experiment 1 PBS total range 10-70", all(as.matrix(exp1$PBS[c("pbs_pre", "pbs_post")]) >= 10 &
                                              as.matrix(exp1$PBS[c("pbs_pre", "pbs_post")]) <= 70))
baseline <- list(); anovas <- list(); ancovas <- list(); regressions <- list()
for (name in names(exp1)) {
  d <- exp1[[name]]
  for (y in outcomes[[name]]) {
    check(paste(y, "change = post - pre"),
          max(abs(d[[paste0(y, "_post")]] - d[[paste0(y, "_pre")]] - d[[paste0(y, "_change")]])) < 1e-8)
    save_table(describe(d, paste0(y, c("_pre", "_post", "_change"))), paste0("exp1_", y, "_descriptives"))
    for (phase in c("pre", "post", "change")) {
      m <- fit(paste0(y, "_", phase), "group", d)
      anovas[[length(anovas) + 1L]] <- cbind(outcome = y, phase = phase, model_summary(m))
      if (phase != "pre") save_table(dunnett(m), paste0("exp1_", y, "_", phase, "_dunnett"))
    }
    # Only PBS and CFI ANCOVAs are specified in the manuscript.
    if (name != "BCIS") {
      m <- fit(paste0(y, "_post"), c("group", paste0(y, "_pre")), d)
      ancovas[[length(ancovas) + 1L]] <- cbind(outcome = y, term_test(m, "group"))
      save_table(coef_table(m), paste0("exp1_", y, "_ancova_coefficients"))
      save_table(dunnett(m), paste0("exp1_", y, "_ancova_dunnett"))
    }
    for (phase in c("post", "change")) {
      m <- fit(paste0(y, "_", phase), c("group", if (phase == "post") paste0(y, "_pre"), covariates), d)
      regressions[[length(regressions) + 1L]] <- cbind(outcome = y, phase = phase, model_summary(m))
      save_table(coef_table(m), paste0("exp1_", y, "_", phase, "_regression_coefficients"))
      tests <- lapply(attr(terms(m), "term.labels"), function(k) term_test(m, k))
      save_table(do.call(rbind, tests), paste0("exp1_", y, "_", phase, "_regression_tests"))
    }
  }
}
for (v in covariates) baseline[[length(baseline) + 1L]] <- cbind(outcome = v, model_summary(fit(v, "group", exp1$PBS)))
save_table(do.call(rbind, baseline), "exp1_recipient_baseline_anova")
save_table(do.call(rbind, anovas), "exp1_outcome_anova")
save_table(do.call(rbind, ancovas), "exp1_ancova_condition_tests")
save_table(do.call(rbind, regressions), "exp1_regression_model_summaries")

message("Experiment 2: separate-slopes ANCOVA and gain-score models")
d <- read_data("Exp2", "PBS_exp2_for_R.csv")
check("Experiment 2 762 participants; 254 per condition", nrow(d) == 762L &&
        setequal(d$group, g2) && all(table(d$group) == 254L))
check("Experiment 2 topic allocation", setequal(d$theme, topics) &&
        all(table(factor(d$theme, levels = topics), factor(d$group, levels = g2)) == c(88, 83, 83)))
check("Experiment 2 change = post - pre", max(abs(d$pbs_post - d$pbs_pre - d$pbs_change)) < 1e-8)
check("Experiment 2 PBS total range 10-70", all(as.matrix(d[c("pbs_pre", "pbs_post")]) >= 10 &
                                              as.matrix(d[c("pbs_pre", "pbs_post")]) <= 70))
d$group <- factor(d$group, levels = g2); d$theme <- factor(d$theme, levels = topics)
contrasts(d$group) <- contr.sum(3); contrasts(d$theme) <- contr.sum(3)
d$pre_c <- d$pbs_pre - mean(d$pbs_pre)
common <- lm(pbs_post ~ pre_c + group * theme, d, na.action = na.fail)
full <- lm(pbs_post ~ pre_c * group * theme, d, na.action = na.fail)
gain <- lm(pbs_change ~ group * theme, d, na.action = na.fail)
save_table(describe(d, c("pbs_pre", "pbs_post", "pbs_change"), c("theme", "group")), "exp2_cell_descriptives")
save_table(data.frame(pretest_grand_mean = mean(d$pbs_pre)), "exp2_centering")
slopes <- anova(common, full)
# Koenker/studentized Breusch-Pagan LM = n * R^2 of squared-residual regression;
# this is the default robust=True calculation in the original Python analysis.
x <- model.matrix(common); e2 <- residuals(common)^2
aux <- lm(e2 ~ x[, -1, drop = FALSE]); bp <- nobs(common) * summary(aux)$r.squared
quad <- update(common, . ~ . + I(pre_c^2)); qtst <- anova(common, quad)
save_table(data.frame(test = c("slope_heterogeneity", "studentized_Breusch_Pagan", "quadratic_pretest"),
                     statistic = c(slopes$F[2], bp, qtst$F[2]),
                     df_num = c(slopes$Df[2], ncol(x) - 1L, qtst$Df[2]),
                     df_denom = c(df.residual(full), NA, df.residual(quad)),
                     p = c(slopes$`Pr(>F)`[2], pchisq(bp, ncol(x) - 1L, lower.tail = FALSE), qtst$`Pr(>F)`[2])),
           "exp2_diagnostics")
grid <- expand.grid(theme = topics, group = g2, KEEP.OUT.ATTRS = FALSE)
grid$group <- factor(grid$group, levels = g2); grid$theme <- factor(grid$theme, levels = topics)
grid$pre_c <- 0
pairs <- list(c("combined", "conceptual"), c("conceptual", "control"), c("combined", "control"))
for (name in c("primary", "gain", "common_slope_sensitivity")) {
  m <- switch(name, primary = full, gain = gain, common_slope_sensitivity = common)
  V <- hc3(m)
  X <- model.matrix(delete.response(terms(m)), grid, contrasts.arg = m$contrasts, xlev = m$xlevels)
  L <- do.call(rbind, lapply(g2, function(g) colMeans(X[grid$group == g, , drop = FALSE])))
  rownames(L) <- g2
  save_table(cbind(group = g2, linear_estimate(m, L, V)), paste0("exp2_", name, "_means"))
  C <- do.call(rbind, lapply(pairs, function(p) L[p[1], ] - L[p[2], ]))
  z <- cbind(comparison = vapply(pairs, paste, character(1), collapse = " - "), linear_estimate(m, C, V))
  z$p_holm <- p.adjust(z$p, "holm")
  save_table(z, paste0("exp2_", name, "_contrasts"))
  tests <- lapply(c("group", "theme", "group:theme"), function(k) term_test(m, k, V))
  save_table(do.call(rbind, tests), paste0("exp2_", name, "_omnibus"))
  if (name == "primary") {
    within <- list()
    for (topic in topics) for (p in pairs) {
      contrast <- X[grid$group == p[1] & grid$theme == topic, ] - X[grid$group == p[2] & grid$theme == topic, ]
      within[[length(within) + 1L]] <- cbind(theme = topic, comparison = paste(p, collapse = " - "),
                                            linear_estimate(m, matrix(contrast, nrow = 1), V))
    }
    z <- do.call(rbind, within); z$p_holm_9 <- p.adjust(z$p, "holm")
    save_table(z, "exp2_primary_within_topic_contrasts")
  }
}
png(file.path(out, "exp2_common_slope_diagnostics.png"), width = 1400, height = 650, res = 150)
par(mfrow = c(1, 2)); plot(fitted(common), residuals(common), xlab = "Fitted PBS posttest", ylab = "Residual")
abline(h = 0, lty = 2); qqnorm(residuals(common)); qqline(residuals(common)); dev.off()

message("Open-ended responses: read final codes without recoding")
read_codes <- function(e) {
  f <- file.path(root, "Exp_Data", paste0("Exp", e), paste0("exp", e, "_Open-ended_responses_coding.xlsx"))
  z <- as.data.frame(readxl::read_excel(f, sheet = "最终编码"))
  names(z)[match(c("编号", "条件", "题目", "参与者原回答", "最终主码"), names(z))] <- c("id", "group", "question", "response", "code")
  mapping <- c("基线" = "control", "概念路径" = "conceptual", "情感路径" = "emotional",
               "行为路径" = "behavioral", "环境路径" = "environmental", "组合路径" = "combined")
  z$group <- unname(mapping[z$group])
  if (e == 2L) {
    z$theme <- unname(c("拖延" = "procrastination", "社交媒体与心理健康" = "social_media",
                       "学习中的一心多用" = "multitasking")[z[["主题"]]])
    check("Experiment 2 open-response topic labels", !anyNA(z$theme))
  }
  z$response[is.na(z$response)] <- ""
  # Literal 'none' and 'no change' answers are retained; only empty cells excluded.
  z$nonblank <- nzchar(gsub("[[:space:]\u3000\u00a0]", "", z$response, perl = TRUE))
  check(paste("Experiment", e, "open-response groups and nonblank codes"),
        !anyNA(z$group) && all(!is.na(z$code[z$nonblank]) & nzchar(z$code[z$nonblank])))
  key <- if (e == 1L) paste(z$id, z$question) else paste(z$theme, z$group, z$id, z$question)
  check(paste("Experiment", e, "unique response-question keys"), !anyDuplicated(key))
  z
}
o1 <- read_codes(1); o2 <- read_codes(2)
q1 <- c("最有帮助的部分", "不适用的部分", "可能做出的改变", "AI再生成建议")
q2 <- c("最有帮助的内容", "不舒服、反感或不适用", "看待方式是否变化", "变化或不变的原因")
og1 <- c("control", "conceptual", "emotional", "behavioral", "environmental")
check("Experiment 1 final codes: 160 x 4 records", nrow(o1) == 640L && setequal(o1$question, q1) &&
        all(table(o1$question, o1$group) == 32L))
check("Experiment 2 final codes: 762 x 4 records", nrow(o2) == 3048L && setequal(o2$question, q2) &&
        all(table(o2$question, o2$group) == 254L))
code_summaries <- function(z, groups, questions, by_topic = FALSE) {
  rows <- list(); denoms <- list()
  for (q in questions) for (g in groups) for (topic in if (by_topic) topics else "all") {
    s <- z[z$question == q & z$group == g, ]
    if (by_topic) s <- s[s$theme == topic, ]
    valid <- s[s$nonblank, ]; n <- nrow(valid)
    denoms[[length(denoms) + 1L]] <- data.frame(question = q, group = g, theme = topic,
      total_n = nrow(s), nonblank_n = n, response_rate = n / nrow(s))
    counts <- table(valid$code)
    if (length(counts)) rows[[length(rows) + 1L]] <- data.frame(question = q, group = g, theme = topic,
      code = names(counts), count = as.integer(counts), nonblank_n = n, percent = 100 * as.integer(counts) / n)
  }
  list(summary = do.call(rbind, rows), denominators = do.call(rbind, denoms))
}
for (e in 1:2) {
  s <- code_summaries(if (e == 1) o1 else o2, if (e == 1) og1 else g2, if (e == 1) q1 else q2)
  save_table(s$summary, paste0("exp", e, "_open_code_frequencies"))
  save_table(s$denominators, paste0("exp", e, "_open_response_rates"))
}
s <- code_summaries(o2, g2, q2, TRUE)
save_table(s$summary, "exp2_open_topic_code_frequencies")
save_table(s$denominators, "exp2_open_topic_response_rates")
fisher_counts <- function(a_yes, a_n, b_yes, b_n) {
  tab <- rbind(c(a_yes, a_n - a_yes), c(b_yes, b_n - b_yes))
  test <- fisher.test(tab, alternative = "two.sided")
  data.frame(a_count = a_yes, a_n = a_n, b_count = b_yes, b_n = b_n,
             a_percent = 100 * a_yes / a_n, b_percent = 100 * b_yes / b_n,
             odds_ratio_conditional = unname(test$estimate), p = test$p.value)
}
targets <- c("理论概念或因果机制", "情绪接纳或自我宽恕", "行动步骤或自我提问", "环境压力或外部情境")
z <- o1[o1$question == q1[1] & o1$nonblank, ]; tests <- list()
for (i in seq_along(targets)) {
  a <- z[z$group == og1[i + 1L], ]; b <- z[z$group == "control", ]
  tests[[i]] <- cbind(condition = og1[i + 1L], code = targets[i],
    fisher_counts(sum(a$code == targets[i]), nrow(a), sum(b$code == targets[i]), nrow(b)))
}
z <- do.call(rbind, tests); z$p_holm_4 <- p.adjust(z$p, "holm")
save_table(z, "exp1_open_targeted_fisher")
positive <- function(q, code) switch(as.character(q),
  "1" = !code %in% c("无帮助", "笼统或无法判断"),
  "2" = code != "无明显问题",
  "3" = code %in% c("明确变化", "有限或保留性的变化", "表示有变化但未说明", "有限变化：原先已经认同"),
  "4" = !code %in% c("未提供原因", "笼统或无法判断"))
binary_counts <- list(); omnibus <- list(); pair_tests <- list()
for (family in c("response_rate", "content_nonblank")) {
  for (q in seq_along(q2)) {
    tab <- matrix(0, 3, 2, dimnames = list(g2, c("yes", "no")))
    for (g in g2) {
      s <- o2[o2$question == q2[q] & o2$group == g, ]
      if (family == "response_rate") vals <- s$nonblank else {
        s <- s[s$nonblank, ]; vals <- positive(q, s$code)
      }
      tab[g, ] <- c(sum(vals), length(vals) - sum(vals))
    }
    binary_counts[[length(binary_counts) + 1L]] <- data.frame(family = family, question = q2[q],
      group = g2, yes = tab[, 1], no = tab[, 2], denominator = rowSums(tab), percent = 100 * tab[, 1] / rowSums(tab))
    test <- chisq.test(tab, correct = FALSE)
    omnibus[[length(omnibus) + 1L]] <- data.frame(family = family, question = q2[q],
      chi2 = unname(test$statistic), df = unname(test$parameter), p = test$p.value,
      cramers_v = sqrt(unname(test$statistic) / sum(tab)), min_expected = min(test$expected))
    pt <- lapply(pairs, function(p) cbind(comparison = paste(p, collapse = " - "),
      fisher_counts(tab[p[1], 1], sum(tab[p[1], ]), tab[p[2], 1], sum(tab[p[2], ]))))
    pt <- do.call(rbind, pt); pt$p_holm_3 <- p.adjust(pt$p, "holm")
    pair_tests[[length(pair_tests) + 1L]] <- cbind(family = family, question = q2[q], pt)
  }
}
z <- do.call(rbind, omnibus)
z$p_holm_4 <- ave(z$p, z$family, FUN = function(p) p.adjust(p, "holm"))
save_table(z, "exp2_open_binary_omnibus")
save_table(do.call(rbind, binary_counts), "exp2_open_binary_counts")
save_table(do.call(rbind, pair_tests), "exp2_open_binary_pairwise")
# Exploratory, selected Q1 categories: combined vs control, one Holm family of 3.
categories <- c("行动建议与自我检验", "情绪理解与支持", "多种功能整合")
s <- o2[o2$question == q2[1] & o2$nonblank, ]
a <- s[s$group == "combined", ]; b <- s[s$group == "control", ]
z <- do.call(rbind, lapply(categories, function(code) cbind(code = code,
  fisher_counts(sum(a$code == code), nrow(a), sum(b$code == code), nrow(b)))))
z$p_holm_3 <- p.adjust(z$p, "holm")
save_table(z, "exp2_open_selected_categories")

save_table(do.call(rbind, checks), "input_validation")
inputs <- list.files(file.path(root, "Exp_Data"), pattern = "(for_R\\.csv|coding\\.xlsx)$",
                     recursive = TRUE, full.names = TRUE)
save_table(data.frame(file = substring(inputs, nchar(root) + 2L),
                       md5 = unname(tools::md5sum(inputs))), "input_manifest")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
message("Completed. Results saved to: ", normalizePath(out))
