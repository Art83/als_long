###################################################################################
# Clinical profile of fast, middle and slow progressors
###################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lme4); library(ggplot2); library(grid)
})

CLIN <- "../raw/als-metadata/answerals/clinical/"
NFL_FILE <- "../raw/NfL/NfL_Seimens_Atellica_Immunoassay_Proteomics_2022_07.csv"
C9_REPEAT_CUTOFF <- 30
N_BOOT <- 2000

num <- function(x) suppressWarnings(as.numeric(x))
txt <- function(x) ifelse(is.na(x), "", trimws(as.character(x)))
rd  <- function(f) read.csv(paste0(CLIN, f), check.names = FALSE)

load("../output/validation/fast_slow_tracking.RData")   # pt, comp, cv
portal <- read.csv("../raw/als-metadata/answerals/aals_dataportal_datatable.csv",
                   check.names = FALSE)
als_guids <- portal$GUID[portal$SUBJECT_GROUP == "ALS"]
ids <- pt %>% select(GUID, Participant_ID, tertile, als_slope, died) %>%
  mutate(tertile = factor(tertile, levels = c("Fast", "Middle", "Slow")))
cat("Patients:", nrow(ids), "\n")

# Baseline characteristics
base <- portal %>% filter(GUID %in% ids$GUID) %>% distinct(GUID, .keep_all = TRUE) %>%
  transmute(GUID,
            male = SEX == "Male",
            bulbar_onset = ifelse(txt(Site_of_Onset) == "", NA, txt(Site_of_Onset) == "Bulbar"),
            el_escorial_definite = ifelse(txt(REVISED_EL_ESCORIAL_CRITERIA) == "", NA,
                                          txt(REVISED_EL_ESCORIAL_CRITERIA) == "Definite"),
            age_onset = num(AGE_AT_SYMPTOM_ONSET),
            cbs_baseline = num(CBS_Baseline_Value),
            c9_repeat = num(C9orf72_repeat_length),
            gene_info = txt(HAS_GENOMICS) == "Yes" | txt(ClinReport_Mutations_Details) != "" |
              txt(HasVariant) != "",
            c9 = (is.finite(c9_repeat) & c9_repeat >= C9_REPEAT_CUTOFF) |
              grepl("C9", ClinReport_Mutations_Details, ignore.case = TRUE) |
              grepl("C9", HasVariant, ignore.case = TRUE),
            any_als_gene = c9 | txt(ClinReport_Mutations_Details) != "" | txt(HasVariant) != "") %>%
  mutate(c9 = ifelse(gene_info, c9, NA), any_als_gene = ifelse(gene_info, any_als_gene, NA)) %>%
  select(-gene_info, -c9_repeat)

dem <- rd("Demographics.csv") %>% transmute(GUID = SubjectUID, age_screening = num(age)) %>%
  filter(is.finite(age_screening)) %>% distinct(GUID, .keep_all = TRUE)

hx <- rd("AALSHXFX.csv") %>% transmute(GUID = SubjectUID, onset_day = num(onsetdt)) %>%
  filter(is.finite(onset_day), onset_day <= 0) %>% distinct(GUID, .keep_all = TRUE) %>%
  mutate(symptom_years = -onset_day / 365.25) %>% select(GUID, symptom_years)

# ALSFRS-R with subdomains
fr <- rd("ALSFRS_R.csv")
if ("alsfdn" %in% names(fr)) fr <- fr %>% filter(is.na(alsfdn) | alsfdn != 2)
items <- c(paste0("alsfrs", 1:4), "alsfrs5a", "alsfrs5b", paste0("alsfrs", 6:9), paste0("alsfrsr", 1:3))
fr <- fr %>% mutate(across(all_of(items), num)) %>%
  mutate(item5 = ifelse(!is.na(alsfrs5a), alsfrs5a, alsfrs5b),
         bulbar = alsfrs1 + alsfrs2 + alsfrs3,
         fine_motor = alsfrs4 + item5 + alsfrs6,
         gross_motor = alsfrs7 + alsfrs8 + alsfrs9,
         respiratory = alsfrsr1 + alsfrsr2 + alsfrsr3,
         total = bulbar + fine_motor + gross_motor + respiratory,
         total = ifelse(is.na(total), num(alsfrst), total),
         GUID = SubjectUID, day = num(Visit_Date)) %>%
  filter(GUID %in% als_guids, is.finite(day), is.finite(total)) %>%
  select(GUID, day, total, bulbar, fine_motor, gross_motor, respiratory) %>%
  arrange(GUID, day) %>% distinct(GUID, day, .keep_all = TRUE)

fr_bl <- fr %>% group_by(GUID) %>% slice(1) %>% ungroup() %>%
  transmute(GUID, alsfrs_baseline = total, bulbar_baseline = bulbar, fine_baseline = fine_motor,
            gross_baseline = gross_motor, resp_baseline = respiratory)

# Vital capacity, best % predicted per visit (FVC, else SVC)
vc <- rd("Vital_Capacity.csv") %>%
  mutate(across(c(fvcp1, fvcp2, fvcp3, svcp1, svcp2, svcp3), num)) %>%
  mutate(fvc = suppressWarnings(pmax(fvcp1, fvcp2, fvcp3, na.rm = TRUE)),
         svc = suppressWarnings(pmax(svcp1, svcp2, svcp3, na.rm = TRUE)),
         fvc = ifelse(is.finite(fvc), fvc, NA), svc = ifelse(is.finite(svc), svc, NA),
         vc_pct = ifelse(!is.na(fvc), fvc, svc),
         GUID = SubjectUID, day = num(Visit_Date)) %>%
  filter(GUID %in% als_guids, is.finite(day), is.finite(vc_pct)) %>%
  select(GUID, day, vc_pct) %>% arrange(GUID, day) %>% distinct(GUID, day, .keep_all = TRUE)
vc_bl <- vc %>% group_by(GUID) %>% slice(1) %>% ungroup() %>% transmute(GUID, vc_baseline = vc_pct)

# BMI
vs <- rd("Vital_Signs.csv") %>%
  transmute(GUID = SubjectUID, day = num(Visit_Date), bmi = num(bmi)) %>%
  filter(GUID %in% als_guids, is.finite(day), is.finite(bmi), bmi > 10, bmi < 60) %>%
  arrange(GUID, day) %>% distinct(GUID, day, .keep_all = TRUE)
bmi_bl <- vs %>% group_by(GUID) %>% slice(1) %>% ungroup() %>% transmute(GUID, bmi_baseline = bmi)

# Serum NfL at first draw
nfl_bl <- NULL
if (file.exists(NFL_FILE)) {
  nf <- read.csv(NFL_FILE, check.names = FALSE)
  if ("Analyte" %in% names(nf)) nf <- nf[grepl("nfl|neurofilament", nf$Analyte, ignore.case = TRUE), ]
  nfl_bl <- nf %>%
    transmute(GUID = as.character(SubjID), day = num(`Days from Baseline`),
              nfl = num(gsub("[^0-9.]", "", `Result Value`))) %>%
    filter(is.finite(nfl), nfl > 0, is.finite(day)) %>%
    arrange(GUID, day) %>% group_by(GUID) %>% slice(1) %>% ungroup() %>%
    transmute(GUID, nfl_baseline = nfl)
}

# Disease-modifying drugs ever recorded
meds <- rd("ANSWER_ALS_Medications_Log.csv") %>%
  transmute(GUID = SubjectUID, m = tolower(txt(med))) %>%
  group_by(GUID) %>%
  summarise(riluzole = any(grepl("rilu|rilutek|tiglutik|exservan", m)),
            edaravone = any(grepl("edaravone|radicava", m)), .groups = "drop")

# Clinical events, day relative to screening
ev_ft <- rd("Feeding_Tube_Placement.csv") %>%
  transmute(GUID = SubjectUID, d = coalesce(num(ftpadmdt), num(ftpdchdt))) %>%
  group_by(GUID) %>% summarise(ft_day = suppressWarnings(min(d, na.rm = TRUE)), .groups = "drop")
ev_niv <- rd("NIV_Log.csv") %>%
  transmute(GUID = SubjectUID, d = num(nivstrdt), use = num(nivusg) == 1) %>%
  group_by(GUID) %>% summarise(niv_day = suppressWarnings(min(d, na.rm = TRUE)),
                               niv_any = any(is.finite(d) | (!is.na(use) & use)), .groups = "drop")
ev_tr <- rd("Tracheostomy.csv") %>%
  transmute(GUID = SubjectUID, d = num(trchdt)) %>%
  group_by(GUID) %>% summarise(trach_day = suppressWarnings(min(d, na.rm = TRUE)), .groups = "drop")

df <- ids %>%
  left_join(base, by = "GUID") %>% left_join(dem, by = "GUID") %>% left_join(hx, by = "GUID") %>%
  left_join(fr_bl, by = "GUID") %>% left_join(vc_bl, by = "GUID") %>% left_join(bmi_bl, by = "GUID") %>%
  left_join(meds, by = "GUID") %>% left_join(ev_ft, by = "GUID") %>% left_join(ev_niv, by = "GUID") %>%
  left_join(ev_tr, by = "GUID") %>%
  mutate(dfs_pre = ifelse(is.finite(symptom_years) & symptom_years > 0,
                          (48 - alsfrs_baseline) / (symptom_years * 12), NA),
         riluzole = coalesce(riluzole, FALSE), edaravone = coalesce(edaravone, FALSE),
         gastrostomy_at_baseline = is.finite(ft_day) & ft_day < 0,
         gastrostomy_ever = is.finite(ft_day),
         niv_at_baseline = is.finite(niv_day) & niv_day < 0,
         niv_ever = coalesce(niv_any, FALSE),
         trach_ever = is.finite(trach_day))
if (!is.null(nfl_bl)) df <- df %>% left_join(nfl_bl, by = "GUID")

# Tests
fmt <- function(z) formatC(z, digits = 3, format = "fg")
sum_cont <- function(v) {
  x <- df[[v]]; g <- df$tertile; ok <- is.finite(x)
  s <- sapply(levels(g), function(l) {
    y <- x[g == l & ok]
    sprintf("%s [%s, %s] (n = %d)", fmt(median(y)), fmt(quantile(y, 0.25)), fmt(quantile(y, 0.75)), length(y))
  })
  p3 <- tryCatch(kruskal.test(x[ok], g[ok])$p.value, error = function(e) NA_real_)
  p2 <- tryCatch(wilcox.test(x[ok & g == "Fast"], x[ok & g == "Slow"], exact = FALSE)$p.value,
                 error = function(e) NA_real_)
  tibble(variable = v, Fast = s[1], Middle = s[2], Slow = s[3], p_three_groups = p3, p_fast_vs_slow = p2)
}
sum_bin <- function(v) {
  x <- as.logical(df[[v]]); g <- df$tertile; ok <- !is.na(x)
  s <- sapply(levels(g), function(l) {
    y <- x[g == l & ok]
    sprintf("%d/%d (%.0f%%)", sum(y), length(y), 100 * mean(y))
  })
  xf <- factor(x, levels = c(FALSE, TRUE))
  p3 <- tryCatch(fisher.test(table(g[ok], xf[ok]))$p.value, error = function(e) NA_real_)
  fs <- ok & g %in% c("Fast", "Slow")
  p2 <- tryCatch(fisher.test(table(droplevels(g[fs]), xf[fs]))$p.value, error = function(e) NA_real_)
  tibble(variable = v, Fast = s[1], Middle = s[2], Slow = s[3], p_three_groups = p3, p_fast_vs_slow = p2)
}

cont_base <- c("age_onset", "age_screening", "symptom_years", "dfs_pre", "alsfrs_baseline",
               "bulbar_baseline", "fine_baseline", "gross_baseline", "resp_baseline",
               "vc_baseline", "bmi_baseline", "cbs_baseline", "nfl_baseline")
bin_base  <- c("male", "bulbar_onset", "el_escorial_definite", "c9", "any_als_gene",
               "riluzole", "edaravone", "gastrostomy_at_baseline", "niv_at_baseline")
cont_base <- intersect(cont_base, names(df))
tab_base <- bind_rows(lapply(cont_base, sum_cont), lapply(bin_base, sum_bin)) %>%
  mutate(fdr_three_groups = p.adjust(p_three_groups, "BH"))
cat("\nBaseline characteristics by third of ALSFRS-R decline (median [IQR] or n/N):\n")
print(as.data.frame(tab_base %>% mutate(across(starts_with("p_") | starts_with("fdr"), ~ signif(.x, 2)))),
      row.names = FALSE, right = FALSE)

# Course during follow-up
blup_slope <- function(dat, var) {
  x <- dat %>% transmute(GUID, t = day / 365.25, v = .data[[var]]) %>%
    filter(is.finite(v), is.finite(t)) %>% group_by(GUID) %>% filter(n() >= 2, sd(t) > 0) %>% ungroup()
  m <- tryCatch(suppressMessages(lmer(v ~ t + (t | GUID), data = x,
                                      control = lmerControl(optimizer = "bobyqa"))),
                error = function(e) NULL)
  if (is.null(m)) return(tibble(GUID = character(0), slope = numeric(0)))
  ranef(m)$GUID %>% tibble::rownames_to_column("GUID") %>%
    transmute(GUID, slope = t + fixef(m)["t"])
}
slopes <- list(
  bulbar_slope = blup_slope(fr, "bulbar"),
  fine_motor_slope = blup_slope(fr, "fine_motor"),
  gross_motor_slope = blup_slope(fr, "gross_motor"),
  respiratory_slope = blup_slope(fr, "respiratory"),
  vc_slope = blup_slope(vc, "vc_pct"),
  bmi_slope = blup_slope(vs, "bmi"))
for (nm in names(slopes)) df <- df %>% left_join(slopes[[nm]] %>% rename(!!nm := slope), by = "GUID")

tab_course <- bind_rows(lapply(names(slopes), sum_cont),
                        lapply(c("gastrostomy_ever", "niv_ever", "trach_ever", "died"), sum_bin)) %>%
  mutate(fdr_three_groups = p.adjust(p_three_groups, "BH"))
cat("\nCourse during follow-up by third (slopes per year, BLUP; events ever recorded):\n")
print(as.data.frame(tab_course %>% mutate(across(starts_with("p_") | starts_with("fdr"), ~ signif(.x, 2)))),
      row.names = FALSE, right = FALSE)

# Does the panel track other measures of decline within each third?
boot_r <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 15) return(tibble(r = NA_real_, lo = NA_real_, hi = NA_real_, n = length(x)))
  r0 <- cor(x, y, method = "spearman")
  rb <- replicate(N_BOOT, { i <- sample(length(x), replace = TRUE); cor(x[i], y[i], method = "spearman") })
  tibble(r = r0, lo = unname(quantile(rb, 0.025, na.rm = TRUE)),
         hi = unname(quantile(rb, 0.975, na.rm = TRUE)), n = length(x))
}
outcomes <- c(`ALSFRS-R total` = "als_slope", Bulbar = "bulbar_slope", `Fine motor` = "fine_motor_slope",
              `Gross motor` = "gross_motor_slope", `Respiratory (ALSFRS-R)` = "respiratory_slope",
              `Vital capacity` = "vc_slope", BMI = "bmi_slope")
dt <- df %>% inner_join(cv %>% select(GUID, composite), by = "GUID")
set.seed(5)
track <- bind_rows(lapply(names(outcomes), function(lab) {
  v <- outcomes[[lab]]
  bind_rows(
    dt %>% group_by(tertile) %>% group_modify(~ boot_r(.x$composite, .x[[v]])) %>% ungroup() %>%
      mutate(tertile = as.character(tertile)),
    boot_r(dt$composite, dt[[v]]) %>% mutate(tertile = "All")) %>%
    mutate(outcome = lab)
}))
cat("\nSpearman r of cross-fitted composite change with each slope, by third:\n")
print(as.data.frame(track %>% select(outcome, tertile, r, lo, hi, n)), digits = 3, row.names = FALSE)

# sensitivity for the fast progressor blind spot
cut_tab <- bind_rows(lapply(2:5, function(k) {
  x <- cv %>% select(GUID, composite, als_slope) %>%
    mutate(grp = cut(rank(als_slope, ties.method = "first"), k, labels = FALSE))
  ext <- x %>% filter(grp %in% c(1, k)) %>%
    mutate(grp = factor(ifelse(grp == 1, "Fastest", "Slowest"), levels = c("Slowest", "Fastest")))
  one <- function(d) {
    f <- lm(composite ~ als_slope, data = d); ci <- confint(f)["als_slope", ]
    c(beta = unname(coef(f)[2]), lo = unname(ci[1]), hi = unname(ci[2]),
      r = cor(d$composite, d$als_slope, method = "spearman"), n = nrow(d))
  }
  fa <- one(filter(ext, grp == "Fastest")); sl <- one(filter(ext, grp == "Slowest"))
  p <- summary(lm(composite ~ als_slope * grp, data = ext))$coefficients["als_slope:grpFastest", "Pr(>|t|)"]
  tibble(groups = k, n_per_group = unname(fa["n"]),
         beta_fastest = unname(fa["beta"]), beta_fastest_lo = unname(fa["lo"]), beta_fastest_hi = unname(fa["hi"]),
         beta_slowest = unname(sl["beta"]), r_fastest = unname(fa["r"]), r_slowest = unname(sl["r"]),
         p_fastest_vs_slowest = p)
}))
cat("\nBlind spot across cut-points (cross-fitted composite, fastest against slowest group):\n")
print(as.data.frame(cut_tab), digits = 3, row.names = FALSE)

# Genetic make-up by third (descriptive given that N is not impressive)
genes <- portal %>% filter(GUID %in% ids$GUID) %>% distinct(GUID, .keep_all = TRUE) %>%
  transmute(GUID,
            has_info = txt(HAS_GENOMICS) == "Yes" | txt(ClinReport_Mutations_Details) != "" |
              txt(HasVariant) != "",
            reported = paste(txt(ClinReport_Mutations_Details), txt(HasVariant), sep = ";"),
            c9_rep = num(C9orf72_repeat_length), atxn2_rep = num(ATXN2_repeat_length)) %>%
  inner_join(ids %>% select(GUID, tertile), by = "GUID") %>% filter(has_info)
gene_long <- genes %>%
  mutate(reported = gsub("\\(WGS\\)", "", reported)) %>%
  separate_rows(reported, sep = "[,;/]+") %>%
  mutate(gene = trimws(reported),
         gene = ifelse(grepl("^c9", gene, ignore.case = TRUE), "C9orf72", toupper(gene))) %>%
  filter(gene != "") %>% select(GUID, tertile, gene) %>%
  bind_rows(genes %>% filter(is.finite(c9_rep), c9_rep >= C9_REPEAT_CUTOFF) %>%
              transmute(GUID, tertile, gene = "C9orf72"),
            genes %>% filter(is.finite(atxn2_rep), atxn2_rep >= 27, atxn2_rep <= 33) %>%
              transmute(GUID, tertile, gene = "ATXN2 intermediate (27-33)")) %>%
  distinct(GUID, gene, .keep_all = TRUE)
denom <- genes %>% count(tertile, name = "N")
tab_genes <- gene_long %>% count(gene, tertile) %>%
  complete(gene, tertile = factor(levels(ids$tertile), levels = levels(ids$tertile)), fill = list(n = 0)) %>%
  left_join(denom, by = "tertile") %>%
  mutate(cell = sprintf("%d/%d", n, N)) %>% select(gene, tertile, cell) %>%
  pivot_wider(names_from = tertile, values_from = cell) %>%
  left_join(gene_long %>% count(gene, name = "total"), by = "gene") %>% arrange(desc(total))
cat("\nReported ALS gene variants and repeat expansions by third (patients with genetic information):\n")
print(as.data.frame(tab_genes), row.names = FALSE)
cat("C9orf72 expansion is a WGS repeat length >=", C9_REPEAT_CUTOFF,
    "or a clinical report. Other variants are as reported, without pathogenicity review.\n")

# Identifying fast progressors at enrolment: pre-baseline Delta FS against NfL
third_of <- function(x, labs) cut(rank(x, ties.method = "first", na.last = "keep"), 3, labels = labs)
ds <- dt %>% filter(is.finite(dfs_pre), is.finite(nfl_baseline)) %>%
  mutate(dfs_third = third_of(dfs_pre, c("Low", "Middle", "High")),
         nfl_third = third_of(nfl_baseline, c("Low", "Middle", "High")),
         is_fast = tertile == "Fast")
auc <- function(score, y) {
  r <- rank(score); n1 <- sum(y); n0 <- sum(!y)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
g_both <- glm(is_fast ~ log(dfs_pre) + log(nfl_baseline), data = ds, family = binomial)
cat(sprintf("\nPatients with both Delta FS and NfL: %d\n", nrow(ds)))
cat(sprintf("AUC for the fast third: Delta FS %.2f | baseline NfL %.2f | both (in-sample, optimistic) %.2f\n",
            auc(ds$dfs_pre, ds$is_fast), auc(ds$nfl_baseline, ds$is_fast),
            auc(fitted(g_both), ds$is_fast)))
cat(sprintf("Spearman r of Delta FS with baseline NfL: %.2f\n",
            cor(ds$dfs_pre, ds$nfl_baseline, method = "spearman")))
cat("\nFast-third patients in each third of the stratifier:\n")
print(rbind(`Delta FS` = table(ds$dfs_third[ds$is_fast]), NfL = table(ds$nfl_third[ds$is_fast])))

fit_b <- function(x) {
  f <- lm(composite ~ als_slope, data = x); ci <- confint(f)["als_slope", ]
  tibble(beta = unname(coef(f)[2]), beta_lo = ci[1], beta_hi = ci[2])
}
set.seed(6)
strat <- bind_rows(
  ds %>% group_by(stratum = dfs_third) %>%
    group_modify(~ bind_cols(boot_r(.x$composite, .x$als_slope), fit_b(.x),
                             tibble(median_als_slope = median(.x$als_slope)))) %>%
    ungroup() %>% mutate(stratifier = "Pre-baseline Delta FS"),
  ds %>% group_by(stratum = nfl_third) %>%
    group_modify(~ bind_cols(boot_r(.x$composite, .x$als_slope), fit_b(.x),
                             tibble(median_als_slope = median(.x$als_slope)))) %>%
    ungroup() %>% mutate(stratifier = "Baseline serum NfL"))
strat_p <- sapply(c(dfs_third = "dfs_third", nfl_third = "nfl_third"), function(g) {
  x <- ds; x$grp <- relevel(factor(as.character(x[[g]])), ref = "Low")
  summary(lm(composite ~ als_slope * grp, data = x))$coefficients["als_slope:grpHigh", "Pr(>|t|)"]
})
cat("\nCross-fitted composite tracking ALSFRS-R slope within thirds of each enrolment stratifier:\n")
print(as.data.frame(strat), digits = 3, row.names = FALSE)
cat(sprintf("Interaction p, high against low: Delta FS %.3g | NfL %.3g\n", strat_p[1], strat_p[2]))

fo <- dt %>% filter(tertile == "Fast", !is.na(bulbar_onset)) %>%
  mutate(onset = factor(ifelse(bulbar_onset, "Bulbar", "Limb or other"),
                        levels = c("Limb or other", "Bulbar")))
set.seed(9)
onset_tab <- fo %>% group_by(onset) %>%
  group_modify(~ bind_cols(
    boot_r(.x$composite, .x$als_slope),
    tibble(p = cor.test(.x$composite, .x$als_slope, method = "spearman", exact = FALSE)$p.value))) %>%
  ungroup()
p_onset <- summary(lm(composite ~ als_slope * onset, data = fo))$coefficients["als_slope:onsetBulbar", "Pr(>|t|)"]
cat("\nWithin the fast third, composite r with ALSFRS-R slope by onset site:\n")
print(as.data.frame(onset_tab), digits = 3, row.names = FALSE)
cat("Bulbar against limb or other onset, interaction p:", signif(p_onset, 3), "\n")

# Figures
base_theme <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
tert_cols <- c(Fast = "#B2182B", Middle = "#999999", Slow = "#2166AC", All = "black")

dom <- df %>% select(tertile, bulbar_slope, fine_motor_slope, gross_motor_slope, respiratory_slope) %>%
  pivot_longer(-tertile, names_to = "domain", values_to = "slope") %>% filter(is.finite(slope)) %>%
  mutate(domain = factor(domain, levels = c("bulbar_slope", "fine_motor_slope", "gross_motor_slope",
                                            "respiratory_slope"),
                         labels = c("Bulbar", "Fine motor", "Gross motor", "Respiratory"))) %>%
  group_by(tertile, domain) %>%
  summarise(med = median(slope), q1 = quantile(slope, 0.25), q3 = quantile(slope, 0.75), .groups = "drop")

pa <- ggplot(dom, aes(domain, med, colour = tertile)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = q1, ymax = q3), position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = tert_cols[1:3], name = "Third of decline") +
  labs(tag = "a", x = NULL, y = "Subdomain slope (points per year, median and IQR)") +
  base_theme + theme(legend.position = "bottom")

pb <- ggplot(track %>% filter(is.finite(r)) %>%
               mutate(outcome = factor(outcome, levels = rev(names(outcomes))),
                      tertile = factor(tertile, levels = c("All", "Slow", "Middle", "Fast"))),
             aes(r, outcome, colour = tertile)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = tert_cols, name = NULL, breaks = c("All", "Fast", "Middle", "Slow")) +
  labs(tag = "b", x = "Spearman r with cross-fitted marker composite change", y = NULL) +
  base_theme + theme(legend.position = "bottom")

pc <- ggplot(strat %>% mutate(stratum = factor(stratum, levels = c("Low", "Middle", "High"))),
             aes(stratum, r, colour = stratifier)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c(`Pre-baseline Delta FS` = "#B2182B", `Baseline serum NfL` = "black"),
                      name = NULL) +
  labs(tag = "c", x = "Third of the enrolment measure",
       y = "Spearman r, composite change with ALSFRS-R slope") +
  base_theme + theme(legend.position = "bottom")

fig_file <- "../output/figures/fig_fast_progressor_clinical.tiff"
tiff(fig_file, width = 16, height = 5.5, units = "in", res = 300, compression = "lzw")
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3, widths = unit(c(1, 1.1, 0.8), "null"))))
print(pa, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pb, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pc, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
dev.off()
cat("\nSaved", fig_file, "\n")

out <- "../output/validation/"
write.csv(tab_base, paste0(out, "fast_progressor_baseline.csv"), row.names = FALSE)
write.csv(tab_course, paste0(out, "fast_progressor_course.csv"), row.names = FALSE)
write.csv(track, paste0(out, "fast_progressor_tracking.csv"), row.names = FALSE)
write.csv(tab_genes, paste0(out, "fast_progressor_genes.csv"), row.names = FALSE)
write.csv(strat, paste0(out, "fast_progressor_enrolment_strata.csv"), row.names = FALSE)
write.csv(cut_tab, paste0(out, "fast_progressor_cutpoints.csv"), row.names = FALSE)
write.csv(onset_tab, paste0(out, "fast_progressor_onset_site.csv"), row.names = FALSE)
save(df, tab_base, tab_course, track, tab_genes, strat, cut_tab, file = paste0(out, "fast_progressor_clinical.RData"))

