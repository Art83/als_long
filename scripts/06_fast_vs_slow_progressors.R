###################################################################################
# Tracking of functional decline by the progression markers across rates of decline
# and across baseline serum NfL. 
# !!!Note: the scripts uses few .rds files that come from the repo-companion with 
# simulations on the same dataset (github.com/Art83/als_sim)
###################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lme4); library(mgcv)
  library(ggplot2); library(grid)
})

CLIN <- "../raw/als-metadata/answerals/clinical/"
NFL_FILE <- "../raw/NfL/NfL_Seimens_Atellica_Immunoassay_Proteomics_2022_07.csv"
NFL_ID <- "SubjID"
NFL_VAL <- "Result Value"
NFL_DAY <- "Days from Baseline"
PROT_LONG_FILE <- "../output/validation/prot_long.rds" # this one wasn't produced in the current repo , see above
TERMINAL_WINDOWS_DAYS <- c(182.625, 365.25)
N_CROSSFIT <- 100


loaded <- load("../output/validation/validation_results.RData")
prot <- readRDS("../output/validation/protein_blup_slopes.rds") # this one wasn't produced in the current repo , see above
mapper <- read.csv("../raw/als-metadata/answerals/aals_dataportal_datatable.csv",
                   check.names = FALSE)
alsfrs <- read.csv(paste0(CLIN, "ALSFRS_R.csv"))
mort <- read.csv(paste0(CLIN, "Mortality.csv"), check.names = FALSE)

als_guids <- mapper$GUID[mapper$SUBJECT_GROUP == "ALS"]
hits <- res_alsfrs$results$BLUP$SeqId[res_alsfrs$results$BLUP$p_fdr < 0.05]
cat("Progression markers:", length(hits), "\n")

af <- alsfrs %>%
  transmute(GUID = SubjectUID, day = suppressWarnings(as.numeric(Visit_Date)),
            score = suppressWarnings(as.numeric(alsfrst))) %>%
  filter(GUID %in% als_guids, is.finite(day), is.finite(score)) %>%
  group_by(GUID) %>% filter(n() >= 2) %>% ungroup() %>%
  mutate(t = day / 365.25)

als_blup <- function(dat) {
  m <- lmer(score ~ t + (t | GUID), data = dat, control = lmerControl(optimizer = "bobyqa"))
  ranef(m)$GUID %>% tibble::rownames_to_column("GUID") %>%
    transmute(GUID, als_slope = t + fixef(m)["t"])
}
als <- als_blup(af)

deaths <- mort %>%
  transmute(Participant_ID, death_day = suppressWarnings(as.numeric(dieddt))) %>%
  filter(is.finite(death_day)) %>% distinct(Participant_ID, .keep_all = TRUE)

pt <- tibble(GUID = longitudinal_guids) %>%
  inner_join(als, by = "GUID") %>%
  left_join(mapper %>% select(GUID, Participant_ID) %>% distinct(GUID, .keep_all = TRUE),
            by = "GUID") %>%
  left_join(deaths, by = "Participant_ID") %>%
  mutate(died = is.finite(death_day),
         tertile = cut(rank(als_slope, ties.method = "first"), 3,
                       labels = c("Fast", "Middle", "Slow")))

cat("Patients:", nrow(pt), "| recorded deceased during follow-up:", sum(pt$died), "\n")
print(table(third = pt$tertile, recorded_deceased = pt$died))
print(pt %>% group_by(tertile) %>%
        summarise(n = n(), median_slope = median(als_slope), sd_slope = sd(als_slope),
                  .groups = "drop") %>% as.data.frame(), digits = 3, row.names = FALSE)

# Helpers 
make_comp <- function(sl, sg) {
  sl %>% inner_join(sg, by = "SeqId") %>%
    group_by(SeqId) %>% mutate(z = sgn * as.numeric(scale(prot_slope))) %>% ungroup() %>%
    group_by(GUID) %>% summarise(composite = mean(z, na.rm = TRUE), .groups = "drop")
}
fit_beta <- function(x) {
  if (nrow(x) < 10) return(tibble(beta = NA_real_, lo = NA_real_, hi = NA_real_, n = nrow(x)))
  f <- lm(composite ~ als_slope, data = x); ci <- confint(f)["als_slope", ]
  tibble(beta = unname(coef(f)["als_slope"]), lo = ci[1], hi = ci[2], n = nrow(x))
}
by_grp <- function(dd, g) {
  dd %>% group_by(across(all_of(g))) %>% group_modify(~ fit_beta(.x)) %>% ungroup()
}
inter_p <- function(dd, g, ref, hi) {
  dd$grp <- relevel(factor(as.character(dd[[g]])), ref = ref)
  cf <- summary(lm(composite ~ als_slope * grp, data = dd))$coefficients
  nm <- paste0("als_slope:grp", hi)
  if (nm %in% rownames(cf)) unname(cf[nm, "Pr(>|t|)"]) else NA_real_
}
boot_r <- function(x, y, B = 2000) {
  r0 <- cor(x, y, method = "pearson")
  rb <- replicate(B, { i <- sample(length(x), replace = TRUE)
  cor(x[i], y[i], method = "pearson") })
  tibble(r = r0, lo = unname(quantile(rb, 0.025)), hi = unname(quantile(rb, 0.975)))
}

# 1. In-sample composite
d <- prot %>% filter(SeqId %in% hits, GUID %in% pt$GUID) %>%
  inner_join(pt %>% select(GUID, als_slope, tertile), by = "GUID")
signs <- d %>% group_by(SeqId) %>%
  summarise(sgn = sign(cor(prot_slope, als_slope)), .groups = "drop")
comp <- make_comp(d, signs) %>% inner_join(pt, by = "GUID")

# 2. Tracking by third of decline
tert_ins <- by_grp(comp, "tertile") %>% mutate(set = "In-sample markers")
p_ins <- inter_p(comp, "tertile", "Slow", "Fast")
cat("\nComposite change per ALSFRS-R point, by third of decline (in-sample markers):\n")
print(as.data.frame(tert_ins), digits = 3, row.names = FALSE)
cat("Interaction p, fast against slow:", signif(p_ins, 3), "\n")

# 3. Follow-up length
samp <- af %>% filter(GUID %in% pt$GUID) %>% group_by(GUID) %>%
  summarise(n_alsfrs = n(), span_years = diff(range(t)), .groups = "drop")
comp <- comp %>% left_join(samp, by = "GUID")
cat("\nFollow-up by third:\n")
print(comp %>% group_by(tertile) %>%
        summarise(mean_alsfrs_visits = mean(n_alsfrs), mean_span_years = mean(span_years),
                  .groups = "drop") %>% as.data.frame(), digits = 3, row.names = FALSE)

comp$grp <- relevel(factor(as.character(comp$tertile)), ref = "Slow")
cs <- summary(lm(composite ~ als_slope * (grp + span_years), data = comp))$coefficients
cat("\nTracking with follow-up length as a modifier:\n")
print(signif(cs[grepl("^als_slope:", rownames(cs)), , drop = FALSE], 3))
comp$grp <- NULL

comp_1y <- comp %>% filter(span_years >= 1)
tert_1y <- by_grp(comp_1y, "tertile")
cat("\nPatients with >= 1 year of follow-up (partly selects on survival):\n")
print(as.data.frame(tert_1y), digits = 3, row.names = FALSE)
cat("Interaction p, fast against slow:", signif(inter_p(comp_1y, "tertile", "Slow", "Fast"), 3), "\n")

# 4. Cross-fitted composite
W <- prot %>% filter(GUID %in% pt$GUID) %>% select(GUID, SeqId, prot_slope) %>%
  pivot_wider(names_from = SeqId, values_from = prot_slope)
ids <- W$GUID; M <- as.matrix(W[, -1])
y <- pt$als_slope[match(ids, pt$GUID)]
n <- length(ids); n_top <- length(hits)
acc <- matrix(NA_real_, n, N_CROSSFIT); overlap <- numeric(0)
set.seed(2)
for (b in seq_len(N_CROSSFIT)) {
  fold <- sample(rep(1:2, length.out = n))
  for (k in 1:2) {
    tr <- fold != k; te <- fold == k
    r_tr <- suppressWarnings(cor(M[tr, ], y[tr], use = "pairwise.complete.obs"))[, 1]
    top <- order(-abs(r_tr))[seq_len(n_top)]
    Z <- scale(M[te, top, drop = FALSE])
    acc[te, b] <- rowMeans(sweep(Z, 2, sign(r_tr[top]), "*"), na.rm = TRUE)
    overlap <- c(overlap, mean(colnames(M)[top] %in% hits))
  }
}
cv <- tibble(GUID = ids, composite = rowMeans(acc, na.rm = TRUE)) %>% inner_join(pt, by = "GUID")
tert_cv <- by_grp(cv, "tertile") %>% mutate(set = "Cross-fitted markers")
p_cv <- inter_p(cv, "tertile", "Slow", "Fast")
r_cv <- cor.test(cv$composite, cv$als_slope, method = "pearson", exact = FALSE)
r_in <- cor.test(comp$composite, comp$als_slope, method = "pearson", exact = FALSE)
cat(sprintf("\nCross-fitting: %d repeats, %d aptamers screened, median overlap with the 72 = %.2f\n",
            N_CROSSFIT, ncol(M), median(overlap)))
cat(sprintf("Pearson r with ALSFRS-R slope: cross-fitted %.3f (p = %.3g) | in-sample %.3f (p = %.3g)\n",
            r_cv$estimate, r_cv$p.value, r_in$estimate, r_in$p.value))
cat("Composite change per ALSFRS-R point, by third of decline (cross-fitted markers):\n")
print(as.data.frame(tert_cv), digits = 3, row.names = FALSE)
cat("Interaction p, fast against slow:", signif(p_cv, 3), "\n")

# 5. Per-marker tracking, fast against slow third
beta_of <- function(x) {
  if (nrow(x) < 10 || sd(x$als_slope) == 0) return(NA_real_)
  unname(coef(lm(prot_slope ~ als_slope, data = x))[2])
}
pm <- d %>% group_by(SeqId) %>% group_modify(~ tibble(
  b_all = beta_of(.x),
  b_fast = beta_of(filter(.x, tertile == "Fast")),
  b_slow = beta_of(filter(.x, tertile == "Slow")))) %>% ungroup() %>%
  mutate(ratio_fast = b_fast / b_all, ratio_slow = b_slow / b_all)
n_lower <- sum(pm$ratio_fast < pm$ratio_slow, na.rm = TRUE)
n_pm <- sum(is.finite(pm$ratio_fast) & is.finite(pm$ratio_slow))
cat(sprintf("\nPer-marker tracking relative to the whole cohort: median fast %.2f | median slow %.2f | %d of %d markers track less in fast\n",
            median(pm$ratio_fast, na.rm = TRUE), median(pm$ratio_slow, na.rm = TRUE), n_lower, n_pm))
cat("Sign test p:", signif(binom.test(n_lower, n_pm)$p.value, 3), "(markers are correlated, so this overstates certainty)\n")

# 6. Curvature (supplementary)
g_lin <- gam(composite ~ als_slope, data = comp, method = "REML")
g_smo <- gam(composite ~ s(als_slope, k = 6), data = comp, method = "REML")
nl_p <- anova(g_lin, g_smo, test = "F")$`Pr(>F)`[2]
cat(sprintf("\nCurvature (supplementary): smooth edf = %.2f | p = %.3g\n", summary(g_smo)$edf, nl_p))
grid_x <- tibble(als_slope = seq(min(comp$als_slope), max(comp$als_slope), length.out = 200))
pr <- predict(g_smo, newdata = grid_x, se.fit = TRUE)
grid_x <- grid_x %>% mutate(fit = pr$fit, lo = fit - 1.96 * pr$se.fit, hi = fit + 1.96 * pr$se.fit)

# 7. Baseline NfL and NfL change
nfl_tab <- NULL; nfl_r <- NULL; nfl_diff <- NULL; both <- NULL; p_nfl <- NA_real_
strip <- function(x) sub("^[A-Za-z]+-", "", as.character(x))
if (file.exists(NFL_FILE)) {
  nf <- read.csv(NFL_FILE, check.names = FALSE)
  if ("Analyte" %in% names(nf)) nf <- nf[grepl("nfl|neurofilament", nf$Analyte, ignore.case = TRUE), ]
  nf <- nf %>%
    transmute(id = strip(.data[[NFL_ID]]),
              day = suppressWarnings(as.numeric(.data[[NFL_DAY]])),
              nfl = suppressWarnings(as.numeric(gsub("[^0-9.]", "", .data[[NFL_VAL]])))) %>%
    filter(is.finite(nfl), nfl > 0, is.finite(day))
  nf_slope <- nf %>% group_by(id) %>% filter(n() >= 2, sd(day) > 0) %>%
    summarise(nfl_slope = unname(coef(lm(log2(nfl) ~ I(day / 365.25)))[2]), .groups = "drop")
  nf_bl <- nf %>% group_by(id) %>% arrange(day, .by_group = TRUE) %>% slice(1) %>% ungroup() %>%
    transmute(id, nfl_bl = nfl)
  
  both <- cv %>% mutate(id = strip(Participant_ID)) %>%
    inner_join(nf_bl, by = "id") %>% inner_join(nf_slope, by = "id") %>%
    mutate(nfl_z = as.numeric(scale(log2(nfl_bl))),
           nfl_slope_z = as.numeric(scale(nfl_slope)),
           comp_z = as.numeric(scale(composite)),
           nfl_third = cut(rank(nfl_z, ties.method = "first"), 3,
                           labels = c("Low", "Middle", "High")))
  cat("\nPatients with baseline NfL and >= 2 NfL draws:", nrow(both), "\n")
  
  r_ns <- cor.test(-both$nfl_slope, both$als_slope, method = "pearson", exact = FALSE)
  r_cs <- cor.test(both$composite, both$als_slope, method = "pearson", exact = FALSE)
  cat(sprintf("Pearson r with ALSFRS-R slope: NfL change %.3f (p = %.3g) | cross-fitted composite %.3f (p = %.3g)\n",
              r_ns$estimate, r_ns$p.value, r_cs$estimate, r_cs$p.value))
  c0 <- lm(als_slope ~ nfl_z + nfl_slope_z, data = both)
  c1 <- lm(als_slope ~ nfl_z + nfl_slope_z + comp_z, data = both)
  cat(sprintf("Variance in ALSFRS-R slope explained: baseline NfL + NfL change %.3f | plus cross-fitted composite %.3f | p = %.3g\n",
              summary(c0)$r.squared, summary(c1)$r.squared, anova(c0, c1)$`Pr(>F)`[2]))
  
  nfl_tab <- both %>% group_by(nfl_third) %>%
    summarise(median_nfl_pg_ml = median(nfl_bl), median_als_slope = median(als_slope), .groups = "drop") %>%
    inner_join(by_grp(both, "nfl_third"), by = "nfl_third")
  p_nfl <- inter_p(both, "nfl_third", "Low", "High")
  cat("\nCross-fitted composite change per ALSFRS-R point, by third of baseline NfL:\n")
  print(as.data.frame(nfl_tab), digits = 3, row.names = FALSE)
  cat("Interaction p, high against low NfL:", signif(p_nfl, 3), "\n")
  
  set.seed(3)
  nfl_r <- both %>% group_by(nfl_third) %>% group_modify(~ bind_rows(
    boot_r(.x$composite, .x$als_slope) %>% mutate(measure = "Plasma marker composite change"),
    boot_r(-.x$nfl_slope, .x$als_slope) %>% mutate(measure = "Serum NfL change"))) %>%
    ungroup()
  cat("\nWithin thirds of baseline NfL, Pearson r of each change measure with ALSFRS-R slope (NfL sign-reversed):\n")
  print(as.data.frame(nfl_r), digits = 3, row.names = FALSE)
  
  # Paired bootstrap of composite r minus NfL change r within each NfL third
  r_diff <- function(x) cor(x$composite, x$als_slope, method = "pearson") -
    cor(-x$nfl_slope, x$als_slope, method = "pearson")
  set.seed(4)
  nfl_diff <- both %>% group_by(nfl_third) %>% group_modify(~ {
    d0 <- r_diff(.x)
    db <- replicate(2000, r_diff(.x[sample(nrow(.x), replace = TRUE), ]))
    tibble(diff = d0, lo = unname(quantile(db, 0.025)), hi = unname(quantile(db, 0.975)),
           p = min(1, 2 * min(mean(db <= 0), mean(db >= 0))), n = nrow(.x))
  }) %>% ungroup()
  d_all <- r_diff(both)
  db_all <- replicate(2000, r_diff(both[sample(nrow(both), replace = TRUE), ]))
  cat("\nComposite r minus NfL change r (paired bootstrap):\n")
  print(as.data.frame(nfl_diff), digits = 3, row.names = FALSE)
  cat(sprintf("All %d patients: difference %.3f | 95%% CI %.3f to %.3f\n", nrow(both), d_all,
              quantile(db_all, 0.025), quantile(db_all, 0.975)))
} else cat("\nNfL file not found:", NFL_FILE, "\n")

# 8. Excluding samples close to death
term_tab <- NULL
if (file.exists(PROT_LONG_FILE)) {
  pl <- readRDS(PROT_LONG_FILE)
  tcol <- intersect(c("time", "time_years", "years", "t", "day", "days"), names(pl))
  vcol <- intersect(c("y", "value", "log2_rfu", "rfu"), names(pl))
  if (length(tcol) == 0 || length(vcol) == 0 || !all(c("GUID", "SeqId") %in% names(pl))) {
    cat("\nprot_long.rds columns not recognised:", paste(names(pl), collapse = ", "), "\n")
  } else {
    pl <- pl %>% filter(SeqId %in% hits, GUID %in% pt$GUID) %>%
      transmute(GUID, SeqId, yv = .data[[vcol[1]]], tv = .data[[tcol[1]]])
    to_days <- if (max(abs(pl$tv), na.rm = TRUE) <= 20) 365.25 else 1
    pl <- pl %>% mutate(day = tv * to_days) %>%
      left_join(pt %>% select(GUID, death_day), by = "GUID")
    gap <- pl %>% filter(is.finite(death_day)) %>% group_by(GUID) %>%
      summarise(gap = first(death_day) - max(day), .groups = "drop")
    cat(sprintf("\nDeath minus last draw (days), recorded deceased: median %.0f | min %.0f | n = %d\n",
                median(gap$gap), min(gap$gap), nrow(gap)))
    if (min(gap$gap) < -30) {
      cat("Draws after the recorded death date, so the clocks differ. Terminal check skipped.\n")
    } else {
      prot_blup <- function(dat) {
        dat %>% group_by(SeqId) %>% group_modify(~ {
          x <- .x %>% group_by(GUID) %>% filter(n() >= 2) %>% ungroup() %>% mutate(ty = day / 365.25)
          m <- tryCatch(lmer(yv ~ ty + (ty | GUID), data = x,
                             control = lmerControl(optimizer = "bobyqa")), error = function(e) NULL)
          if (is.null(m)) return(tibble(GUID = character(0), prot_slope = numeric(0)))
          ranef(m)$GUID %>% tibble::rownames_to_column("GUID") %>%
            transmute(GUID, prot_slope = ty + fixef(m)["ty"])
        }) %>% ungroup()
      }
      term_tab <- bind_rows(lapply(TERMINAL_WINDOWS_DAYS, function(w) {
        pl_w <- pl %>% filter(!is.finite(death_day) | day <= death_day - w)
        af_w <- af %>% left_join(pt %>% select(GUID, death_day), by = "GUID") %>%
          filter(!is.finite(death_day) | day <= death_day - w) %>%
          group_by(GUID) %>% filter(n() >= 2) %>% ungroup()
        als_w <- als_blup(af_w)
        sl <- prot_blup(pl_w) %>% inner_join(als_w, by = "GUID")
        cw <- make_comp(sl, signs) %>% inner_join(als_w, by = "GUID") %>%
          inner_join(pt %>% select(GUID, tertile), by = "GUID")
        by_grp(cw, "tertile") %>%
          mutate(window_months = round(w / 30.44), p_fast_vs_slow = inter_p(cw, "tertile", "Slow", "Fast"))
      }))
      cat("\nExcluding protein and ALSFRS-R samples within the window before death (in-sample markers):\n")
      print(as.data.frame(term_tab), digits = 3, row.names = FALSE)
    }
  }
} else {
  cat("\nTerminal check skipped: output/validation/prot_long.rds not found.\n")
}

# Figures
base <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

a_dat <- bind_rows(tert_ins, tert_cv) %>%
  mutate(set = factor(set, levels = c("In-sample markers", "Cross-fitted markers")))
pa <- ggplot(a_dat %>% filter(is.finite(beta)), aes(tertile, beta, colour = set)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c(`In-sample markers` = "black", `Cross-fitted markers` = "#2166AC"),
                      name = NULL) +
  labs(tag = "a", x = "Third of ALSFRS-R decline",
       y = "Composite change per ALSFRS-R point",
       subtitle = sprintf("Fast against slow interaction\np = %s in-sample, p = %s cross-fitted",
                          format.pval(p_ins, digits = 2), format.pval(p_cv, digits = 2))) +
  base + theme(legend.position = "bottom")

pb <- ggplot(pm %>% filter(is.finite(ratio_fast), is.finite(ratio_slow)),
             aes(ratio_fast, ratio_slow)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40") +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey70") +
  geom_point(alpha = 0.6, size = 1.6) +
  coord_cartesian(xlim = c(-1, 3), ylim = c(-1, 3)) +
  labs(tag = "b", x = "Tracking in fast third (relative to whole cohort)",
       y = "Tracking in slow third (relative to whole cohort)",
       subtitle = sprintf("%d of %d markers above the line of equality", n_lower, n_pm),
       caption = "One point per marker. Axes clipped at -1 and 3") +
  base

if (!is.null(nfl_r)) {
  pc <- ggplot(nfl_r, aes(nfl_third, r, colour = measure)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = c(`Plasma marker composite change` = "black",
                                   `Serum NfL change` = "#B2182B"), name = NULL) +
    labs(tag = "c", x = "Third of baseline serum NfL",
         y = "Pearson r with ALSFRS-R slope",
         subtitle = sprintf("Same %d patients, composite cross-fitted\nDifference p = %s (low), %s (middle), %s (high)",
                            nrow(both), format.pval(nfl_diff$p[1], digits = 2),
                            format.pval(nfl_diff$p[2], digits = 2), format.pval(nfl_diff$p[3], digits = 2))) +
    base + theme(legend.position = "bottom")
} else {
  pc <- ggplot() + theme_void() + labs(tag = "c")
}

fig_file <- "../output/figures/fig_fast_slow_tracking.tiff"
tiff(fig_file, width = 14, height = 5, units = "in", res = 300, compression = "lzw")
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
print(pa, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pb, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pc, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
dev.off()
cat("\nSaved", fig_file, "\n")

ps1 <- ggplot() +
  geom_point(data = comp, aes(als_slope, composite, shape = died), alpha = 0.6, size = 1.6) +
  geom_ribbon(data = grid_x, aes(als_slope, ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(data = grid_x, aes(als_slope, fit), linewidth = 0.9) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                     labels = c("Not recorded deceased", "Recorded deceased"), name = NULL) +
  labs(tag = "a", x = "ALSFRS-R slope (points per year)",
       y = "Composite marker change (z, in-sample markers)",
       subtitle = sprintf("Smooth fit, curvature p = %s", format.pval(nl_p, digits = 2))) +
  base + theme(legend.position = "bottom")

if (!is.null(nfl_tab)) {
  ps2 <- ggplot(nfl_tab %>% filter(is.finite(beta)), aes(nfl_third, beta)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_pointrange(aes(ymin = lo, ymax = hi), colour = "#2166AC") +
    labs(tag = "b", x = "Third of baseline serum NfL",
         y = "Composite change per ALSFRS-R point (cross-fitted)",
         subtitle = sprintf("High against low NfL interaction p = %s",
                            format.pval(p_nfl, digits = 2))) +
    base
} else {
  ps2 <- ggplot() + theme_void() + labs(tag = "b")
}

supp_file <- "../output/figures/fig_fast_slow_supp.tiff"
tiff(supp_file, width = 10, height = 4.8, units = "in", res = 300, compression = "lzw")
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(ps1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ps2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
cat("Saved", supp_file, "\n")

out <- "../output/validation/"
write.csv(bind_rows(tert_ins, tert_cv), paste0(out, "fast_slow_by_third.csv"), row.names = FALSE)
write.csv(pm, paste0(out, "fast_slow_marker_ratios.csv"), row.names = FALSE)
if (!is.null(nfl_tab)) write.csv(nfl_tab, paste0(out, "fast_slow_by_nfl_third.csv"), row.names = FALSE)
if (!is.null(nfl_r))   write.csv(nfl_r, paste0(out, "fast_slow_nfl_change_r.csv"), row.names = FALSE)
if (!is.null(term_tab)) write.csv(term_tab, paste0(out, "fast_slow_terminal_window.csv"), row.names = FALSE)
save(pt, comp, cv, tert_ins, tert_cv, p_ins, p_cv, tert_1y, pm, nfl_tab, nfl_r, nfl_diff, p_nfl, term_tab,
     file = paste0(out, "fast_slow_tracking.RData"))