##############################################################
# Script for validation of progressors results in Target ALS #
##############################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lme4); library(readxl)
  library(ggplot2); library(grid);
})

TMT_FILE <- "../raw/target_als/ms/lc_ms_proteomics_tmt_plasma.xlsx"
N_PERM <- 2000
N_NULL <- 1000
MIN_SUBJ <- 15

num   <- function(x) { x <- suppressWarnings(as.numeric(x)); x[!is.na(x) & x <= -6666] <- NA; x }
score <- function(x) { x <- suppressWarnings(as.numeric(x)); x[!is.na(x) & x < 0] <- NA; x }

# Answer ALS markers and directions
load("../output/validation/validation_results.RData")
aals <- res_alsfrs$results$BLUP %>%
  left_join(seq_map %>% select(SeqId, uniprot), by = "SeqId") %>%
  transmute(SeqId, uniprot = as.character(uniprot), est, pval, p_fdr) %>%
  filter(!is.na(uniprot), !grepl("\\|", uniprot)) %>%
  mutate(uniprot = sub("-\\d+$", "", uniprot)) %>%
  group_by(uniprot) %>% slice_min(pval, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(sign_aals = sign(est), marker = p_fdr < 0.05)
cat("Answer ALS: single-target proteins", nrow(aals), "| progression markers", sum(aals$marker), "\n")

# Target ALS clinical
clin <- read.csv("../raw/target_als/nhs_clinical_data.csv", check.names = FALSE) %>%
  transmute(subject_id = `Subject ID`, visit = as.integer(`Visit Number`),
            group = `Clinical Group`, alsfrs = score(`Alsfrs R Score Total`)) %>%
  filter(group == "ALS")
fx <- clin %>% filter(is.finite(alsfrs), alsfrs <= 48) %>% mutate(t = visit - 1) %>%
  group_by(subject_id) %>% filter(n() >= 2) %>% ungroup()
m_fx <- lmer(alsfrs ~ t + (t | subject_id), data = fx, control = lmerControl(optimizer = "bobyqa"))
als <- ranef(m_fx)$subject_id %>% tibble::rownames_to_column("subject_id") %>%
  transmute(subject_id, als_slope = t + fixef(m_fx)["t"])
cat(sprintf("ALSFRS-R scores used: range %g to %g | per-patient slopes: median %.2f, IQR %.2f to %.2f\n",
            min(fx$alsfrs), max(fx$alsfrs), median(als$als_slope),
            quantile(als$als_slope, 0.25), quantile(als$als_slope, 0.75)))
cat("Target ALS: ALS patients with >= 2 ALSFRS-R:", nrow(als),
    "| mean slope", round(fixef(m_fx)["t"], 2), "points per visit (no visit dates in the release)\n")

# TMT plasma 
info <- read_excel(TMT_FILE, sheet = "1-plasma_sample_info", col_types = "text") %>% as.data.frame()
ratios <- read_excel(TMT_FILE, sheet = "4-plasma_protein_TMT_ratios") %>% as.data.frame()
names(ratios)[1] <- "uniprot"
info <- info %>% filter(technical.replicate != "-8888") %>%
  mutate(key1 = paste0(set, ".", label), key2 = set.label)
key <- if (sum(info$key1 %in% names(ratios)) >= sum(info$key2 %in% names(ratios))) "key1" else "key2"
info$key <- info[[key]]
info <- info %>% filter(key %in% names(ratios))
cat("TMT channels matched to sample info:", nrow(info), "\n")

mat <- as.matrix(ratios[, info$key]); mode(mat) <- "numeric"
mat <- log2(mat); mat[!is.finite(mat)] <- NA
rownames(mat) <- sub("-\\d+$", "", ratios$uniprot)
mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
info$grp <- paste(info$subject.id, info$visit.number, sep = "|")
grps <- unique(info$grp)
mat_v <- sapply(grps, function(g) rowMeans(mat[, info$grp == g, drop = FALSE], na.rm = TRUE))
mat_v[!is.finite(mat_v)] <- NA
samp <- info %>% distinct(grp, .keep_all = TRUE) %>%
  transmute(grp, subject_id = subject.id, visit = as.integer(visit.number)) %>%
  filter(subject_id %in% als$subject_id) %>%
  semi_join(clin %>% distinct(subject_id, visit), by = c("subject_id", "visit")) %>%
  group_by(subject_id) %>% filter(n() >= 2) %>% ungroup() %>% mutate(t = visit - 1)
cat("TMT plasma:", nrow(mat_v), "proteins | ALS patients with >= 2 visits and an ALSFRS-R slope:",
    n_distinct(samp$subject_id), "\n")

prot_blup <- function(u) {
  d <- samp %>% mutate(v = mat_v[u, grp]) %>% filter(is.finite(v)) %>%
    group_by(subject_id) %>% filter(n() >= 2) %>% ungroup()
  if (n_distinct(d$subject_id) < MIN_SUBJ) return(NULL)
  m <- tryCatch(suppressMessages(suppressWarnings(
    lmer(v ~ t + (t | subject_id), data = d, control = lmerControl(calc.derivs = FALSE)))),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  ranef(m)$subject_id %>% tibble::rownames_to_column("subject_id") %>%
    transmute(subject_id, uniprot = u, prot_slope = t + fixef(m)["t"])
}
shared <- intersect(rownames(mat_v), aals$uniprot)
cat("Fitting per-protein slopes for", length(shared), "proteins also measured in Answer ALS...\n")
sl <- bind_rows(lapply(shared, prot_blup)) %>% inner_join(als, by = "subject_id")

# Composite with directions fixed from Answer ALS
make_comp <- function(sl_sub, signs) {
  sl_sub %>% inner_join(signs, by = "uniprot") %>%
    group_by(uniprot) %>% mutate(z = sign_aals * as.numeric(scale(prot_slope))) %>% ungroup() %>%
    group_by(subject_id) %>% summarise(composite = mean(z, na.rm = TRUE), n_prot = n(), .groups = "drop")
}
mk <- aals %>% filter(marker, uniprot %in% unique(sl$uniprot)) %>% select(uniprot, sign_aals, SeqId)
cat("Answer ALS markers fitted in Target ALS:", nrow(mk), "\n")
comp <- make_comp(sl %>% filter(uniprot %in% mk$uniprot), mk) %>%
  filter(n_prot >= ceiling(nrow(mk) / 2)) %>% inner_join(als, by = "subject_id")

r_obs <- cor(comp$composite, comp$als_slope, method = "pearson")
set.seed(7)
r_boot <- replicate(N_PERM, { i <- sample(nrow(comp), replace = TRUE)
                              cor(comp$composite[i], comp$als_slope[i], method = "pearson") })
r_perm <- replicate(N_PERM, cor(comp$composite, sample(comp$als_slope), method = "pearson"))
p_perm <- (1 + sum(r_perm >= r_obs)) / (N_PERM + 1)
cat(sprintf("\nComposite (%d markers, %d patients): Pearson r with ALSFRS-R slope %.3f | 95%% CI %.3f to %.3f | one-sided permutation p = %.3g\n",
            nrow(mk), nrow(comp), r_obs, quantile(r_boot, 0.025), quantile(r_boot, 0.975), p_perm))

# Per-marker correlations and direction agreement
gene_tab <- read_excel(TMT_FILE, sheet = "2-plasma_protein_info") %>% as.data.frame()
genes <- tibble(uniprot = sub("-\\d+$", "", gene_tab[[1]]), gene = sub(";.*$", "", gene_tab$Gene.Symbol)) %>%
  distinct(uniprot, .keep_all = TRUE)
per_marker <- sl %>% filter(uniprot %in% mk$uniprot) %>% group_by(uniprot) %>%
  summarise(r_ta = cor(prot_slope, als_slope, method = "pearson"),
            p_ta = cor.test(prot_slope, als_slope, method = "pearson", exact = FALSE)$p.value,
            n = n(), .groups = "drop") %>%
  inner_join(mk, by = "uniprot") %>% left_join(genes, by = "uniprot") %>%
  mutate(concordant = sign(r_ta) == sign_aals, gene = coalesce(gene, uniprot)) %>%
  arrange(desc(sign_aals * r_ta))
bt <- binom.test(sum(per_marker$concordant), nrow(per_marker), 0.5)
cat(sprintf("Markers in the Answer ALS direction: %d of %d (binomial p = %.3g)\n",
            sum(per_marker$concordant), nrow(per_marker), bt$p.value))
print(as.data.frame(per_marker %>% select(gene, uniprot, sign_aals, r_ta, p_ta, n)), digits = 3, row.names = FALSE)

# Leave one marker out
loo <- bind_rows(lapply(mk$uniprot, function(u) {
  k <- mk %>% filter(uniprot != u)
  cc <- make_comp(sl %>% filter(uniprot %in% k$uniprot), k) %>% inner_join(als, by = "subject_id")
  tibble(dropped = u, r = cor(cc$composite, cc$als_slope, method = "pearson"))
})) %>% left_join(genes, by = c("dropped" = "uniprot"))
cat(sprintf("Leave one marker out: r from %.3f (without %s) to %.3f\n",
            min(loo$r), loo$gene[which.min(loo$r)], max(loo$r)))

# random non-marker proteins with directions from Answer ALS
cat("Fitting slopes for the null set...\n")
nm <- aals %>% filter(!marker, uniprot %in% unique(sl$uniprot)) %>% select(uniprot, sign_aals)
set.seed(8)
r_null <- replicate(N_NULL, {
  k <- nm %>% slice_sample(n = nrow(mk))
  cc <- make_comp(sl %>% filter(uniprot %in% k$uniprot), k) %>% inner_join(als, by = "subject_id")
  cor(cc$composite, cc$als_slope, method = "pearson")
})
p_null <- (1 + sum(r_null >= r_obs)) / (N_NULL + 1)
cat(sprintf("Null of %d random sets of %d non-marker proteins (directions from Answer ALS): median r %.3f | 95th percentile %.3f | empirical p = %.3g\n",
            N_NULL, nrow(mk), median(r_null), quantile(r_null, 0.95), p_null))

# Composite beyond baseline NEFL (Olink)
nefl <- NULL; nefl_res <- NULL
olk_file <- "../raw/target_als/olink_matrix.csv"
if (file.exists(olk_file)) {
  olk <- read.csv(olk_file, check.names = FALSE, row.names = 1)
  if ("NEFL" %in% names(olk)) {
    nefl <- tibble(subject_id = rownames(olk), nefl = num(olk$NEFL)) %>% filter(is.finite(nefl))
    cn <- comp %>% inner_join(nefl, by = "subject_id") %>%
      mutate(nefl_z = as.numeric(scale(nefl)), comp_z = as.numeric(scale(composite)))
    if (nrow(cn) >= 20) {
      r_n <- cor.test(-cn$nefl, cn$als_slope, method = "pearson", exact = FALSE)
      m0 <- lm(als_slope ~ nefl_z, data = cn); m1 <- lm(als_slope ~ nefl_z + comp_z, data = cn)
      nefl_res <- tibble(n = nrow(cn), r_nefl = unname(r_n$estimate), p_nefl = r_n$p.value,
                         r_nefl_pearson = unname(r_n$estimate), p_nefl_pearson = r_n$p.value,
                         r2_nefl = summary(m0)$r.squared, r2_both = summary(m1)$r.squared,
                         p_added = anova(m0, m1)$`Pr(>F)`[2])
      cat(sprintf("\nBaseline NEFL (Olink), %d patients: r with ALSFRS-R slope (sign-reversed) %.3f (p = %.3g)\n",
                  nrow(cn), r_n$estimate, r_n$p.value))
      cat(sprintf("Variance in ALSFRS-R slope explained: NEFL %.3f | NEFL + composite %.3f | composite adds p = %.3g\n",
                  nefl_res$r2_nefl, nefl_res$r2_both, nefl_res$p_added))
    } else cat("\nFewer than 20 patients with both NEFL and the composite\n")
  }
}

# Fast against slow, direction only
half <- comp %>% mutate(half = factor(ifelse(rank(als_slope, ties.method = "first") <= nrow(comp) / 2,
                                             "Faster half", "Slower half"),
                                      levels = c("Slower half", "Faster half")))
hb <- half %>% group_by(half) %>% group_modify(~ {
  f <- lm(composite ~ als_slope, data = .x); ci <- confint(f)["als_slope", ]
  tibble(beta = unname(coef(f)[2]), lo = unname(ci[1]), hi = unname(ci[2]),
         r = cor(.x$composite, .x$als_slope, method = "pearson"), n = nrow(.x),
         median_slope = median(.x$als_slope))
}) %>% ungroup()
p_half <- summary(lm(composite ~ als_slope * half, data = half))$coefficients["als_slope:halfFaster half", "Pr(>|t|)"]
cat("\nFaster against slower half (underpowered, direction only):\n")
print(as.data.frame(hb), digits = 3, row.names = FALSE)
cat("Interaction p:", signif(p_half, 3), "\n")

# Figures
base <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
pa <- ggplot(comp, aes(als_slope, composite)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#B2182B", linewidth = 0.8) +
  labs(tag = "a", x = "ALSFRS-R slope (points per visit)", y = "Marker composite change (z)",
       subtitle = sprintf("Target ALS, %d patients, %d markers\nr = %.2f, permutation p = %s",
                          nrow(comp), nrow(mk), r_obs, format.pval(p_perm, digits = 2))) + base
pb <- ggplot(per_marker %>% mutate(gene = factor(gene, levels = rev(unique(gene))), r_dir = sign_aals * r_ta),
             aes(r_dir, gene, colour = concordant)) +
  geom_vline(xintercept = 0, colour = "grey50") + geom_point(size = 2) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "#B2182B"),
                      labels = c(`TRUE` = "Same direction", `FALSE` = "Opposite"), name = NULL) +
  labs(tag = "b", x = "Pearson r in Target ALS, signed by Answer ALS direction", y = NULL,
       subtitle = sprintf("%d of %d in the Answer ALS direction", sum(per_marker$concordant), nrow(per_marker))) +
  base + theme(legend.position = "bottom")
pc <- ggplot(tibble(r = r_null), aes(r)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "white") +
  geom_vline(xintercept = r_obs, colour = "#B2182B", linewidth = 1) +
  labs(tag = "c", x = "Pearson r of random non-marker composites", y = "Count",
       subtitle = sprintf("Markers against %d random protein sets, p = %s", N_NULL,
                          format.pval(p_null, digits = 2))) + base

fig_file <- "../output/figures/fig_target_als_replication.tiff"
tiff(fig_file, width = 14, height = 5, units = "in", res = 300, compression = "lzw")
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
print(pa, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pb, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pc, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
dev.off()
cat("\nSaved", fig_file, "\n")

out <- "../output/validation/"
write.csv(per_marker, paste0(out, "target_als_per_marker.csv"), row.names = FALSE)
write.csv(loo, paste0(out, "target_als_leave_one_out.csv"), row.names = FALSE)
write.csv(hb, paste0(out, "target_als_halves.csv"), row.names = FALSE)
save(comp, mk, per_marker, loo, r_obs, r_boot, p_perm, r_null, p_null, nefl_res, hb, p_half,
     file = paste0(out, "target_als_replication.RData"))
cat("Saved output/validation/target_als_*.csv and target_als_replication.RData\n")
