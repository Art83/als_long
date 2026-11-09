################################################################################
#  Supplementary figures
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

OUT <- "../output/manuscript_figures/sups/"

cat("================================================================\n")
cat("  10: Supplementary figures\n")
cat("================================================================\n\n")

base <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 7, colour = "grey25"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_rect(fill = "grey92", colour = NA))

tier_levels <- c("confirmed_full", "confirmed_partial", "rna_only_ms_gap", "ms_only_rna_gap", "no_mn_expression")
tier_labels <- c(confirmed_full = "Expressed, promoter-supported",
                 confirmed_partial = "Expressed, promoter-unconfirmed",
                 rna_only_ms_gap = "Transcript only",
                 ms_only_rna_gap = "Protein only",
                 no_mn_expression = "Not expressed")
tier_colors <- c(`Expressed, promoter-supported` = "#2166AC",
                 `Expressed, promoter-unconfirmed` = "#4393C3",
                 `Transcript only` = "grey50",
                 `Protein only` = "#B2ABD2",
                 `Not expressed` = "grey80")
col_up <- "#B2182B"; col_down <- "#2166AC"
route_cols <- c(Secreted = "#1B7837", Membrane = "#A6DBA0", Intracellular = "#C2A5CF", Unknown = "grey85")

fmt_p_plain <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "fg", digits = 2))
fisher_ci <- function(r, n) {
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  data.frame(lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se))
}
save_fig <- function(file, plots, layout, width, height) {
  tiff(paste0(OUT, file), width = width, height = height, units = "in", res = 300, compression = "lzw")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
  for (i in seq_along(plots)) {
    cells <- which(layout == i, arr.ind = TRUE)
    print(plots[[i]], vp = viewport(layout.pos.row = min(cells[, 1]):max(cells[, 1]),
                                    layout.pos.col = min(cells[, 2]):max(cells[, 2])))
  }
  invisible(dev.off())
  cat("  Saved", paste0(OUT, file), "\n")
}

ta_env <- new.env(); load("../output/validation/target_als_replication.RData", envir = ta_env)
fs_env <- new.env(); load("../output/validation/fast_slow_tracking.RData", envir = fs_env)
fc_env <- new.env(); load("../output/validation/fast_progressor_clinical.RData", envir = fc_env)


# =============================================================================
# fig. S1. Target ALS composite: leave one marker out, and baseline NEFL
# =============================================================================
cat("-- fig. S1\n")

loo <- ta_env$loo %>% mutate(gene = reorder(gene, r))
ps1a <- ggplot(loo, aes(r, gene)) +
  geom_vline(xintercept = ta_env$r_obs, colour = col_up, linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_point(size = 1.6) +
  labs(tag = "", x = "Composite correlation with ALSFRS-R slope", y = "Marker removed",
       subtitle = sprintf("All %d markers: r = %.2f (dashed)\nRange after removal: %.2f to %.2f",
                          nrow(ta_env$mk), ta_env$r_obs, min(loo$r), max(loo$r))) +
  base + theme(axis.text.y = element_text(face = "italic", size=9),
               axis.title.x = element_text(size=12),
               axis.title.y = element_text(size=12))

nr <- ta_env$nefl_res
nd <- data.frame(model = factor(c("Baseline NEFL", "Baseline NEFL\n+ composite change"),
                                levels = c("Baseline NEFL", "Baseline NEFL\n+ composite change")),
                 r2 = c(nr$r2_nefl, nr$r2_both))
ps1b <- ggplot(nd, aes(model, 100 * r2)) +
  geom_col(width = 0.55, fill = c("grey55", col_down)) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * r2)), vjust = -0.4, size = 3) +
  scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, 0.02))) +
  labs(tag = "", x = NULL, y = "Variance in ALSFRS-R slope explained (%)",
       subtitle = sprintf("Olink baseline NEFL, %d patients: r = %.2f (sign-reversed), p = %s\nComposite added p = %s",
                          nr$n, nr$r_nefl, fmt_p_plain(nr$p_nefl), fmt_p_plain(nr$p_added))) +
  base + theme(
               axis.text.x = element_text(size=9),
               axis.title.y = element_text(size=12))


save_fig("FigS1.tiff", list(ps1a, ps1b), matrix(1:2, nrow = 1), width = 7.2, height = 3.6)


tiff(paste0(OUT, "FigS1A.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(ps1a)
dev.off()

tiff(paste0(OUT, "FigS1B.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(ps1b)
dev.off()


# =============================================================================
# fig. S2. Stability of the motor neuron expression map
# =============================================================================
cat("-- fig. S2\n")

sens <- read.csv("../output/summaries/fb_tier_threshold_sensitivity.csv", check.names = FALSE) %>%
  pivot_longer(-threshold, names_to = "tier", values_to = "n") %>%
  mutate(tier = factor(unname(tier_labels[tier]), levels = unname(tier_labels)))
ne <- sens %>% filter(tier == "Not expressed")
ps2a <- ggplot(sens, aes(100 * threshold, n, colour = tier)) +
  geom_line(size=0.8) + geom_point(size = 2) +
  scale_colour_manual(values = c(tier_colors[1:4], `Not expressed` = "black"), name = NULL) +
  scale_x_continuous(breaks = c(30, 50, 70, 90)) +
  labs(tag = "", x = "Protein detection threshold\n(% of iMN lines)", y = "Plasma proteins (aptamers)") +
  base + theme(legend.position = "bottom",
               legend.text = element_text(size=9),
               axis.text.x = element_text(size=9),
               axis.title.y = element_text(size=12),
               axis.title.x = element_text(size=12))


# Probably not going to include
rs <- read.csv("../output/tables/rna_threshold_sensitivity.csv", check.names = FALSE) %>%
  mutate(rule = gsub("pct>=", "", rule), rule = gsub(" \\+ vst floor", "\nwith floor", rule),
         rule = gsub(" only", "\nno floor", rule), rule = gsub("edgeR::filterByExpr", "edgeR\nfilterByExpr", rule)) %>%
  mutate(rule = factor(rule, levels = rule)) %>%
  pivot_longer(c(confirmed_full, confirmed_partial, rna_only_ms_gap, ms_only_rna_gap, no_mn_expression),
               names_to = "tier", values_to = "n") %>%
  mutate(tier = factor(unname(tier_labels[tier]), levels = unname(tier_labels)))
ps2b <- ggplot(rs, aes(rule, n, fill = tier)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = tier_colors, name = NULL) +
  labs(tag = "", x = "Transcript detection rule (proportion of iMN lines)", y = "Plasma proteins (aptamers)",
       subtitle = "Primary rule: 0.20 with expression floor") +
  guides(fill = guide_legend(nrow = 2)) +
  base + theme(legend.position = "bottom", 
               legend.text = element_text(size=9),
               axis.text.x = element_text(size = 9),
               axis.title.y = element_text(size=12),
               axis.title.x = element_text(size=12))

gm <- read.csv("../output/grounding/grounding_map.csv", stringsAsFactors = FALSE)
ms_sum <- readRDS("../output/grounding/ms_summary.rds")
imputed_only <- ms_sum$uniprot[ms_sum$ms_in_imputed %in% TRUE & !(ms_sum$ms_in_nonimputed %in% TRUE)]
to <- gm %>% filter(grounding_tier == "rna_only_ms_gap") %>%
  mutate(ms = case_when(ms_measured %in% TRUE ~ "Below threshold",
                        uniprot %in% imputed_only ~ "Imputed only",
                        TRUE ~ "Never quantified"),
         ms = factor(ms, levels = c("Never quantified", "Imputed only", "Below threshold"))) %>%
  count(ms) %>% mutate(pct = 100 * n / sum(n))
print(as.data.frame(to))
to_labels <- setNames(sprintf("%s\n%s (%.1f%%)", as.character(to$ms), format(to$n, big.mark = ",", trim = TRUE), to$pct), as.character(to$ms))
ps2c <- ggplot(to, aes(x = "", y = pct, fill = ms)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("grey85", "grey65", "grey40"), labels = to_labels, name = NULL) +
  labs(tag = "", x = "Transcript-only\nplasma proteins", y = "Proteins (%)",
       subtitle = sprintf("%s aptamers", format(sum(to$n), big.mark = ","))) +
  base + theme(legend.position = "bottom", legend.direction = "vertical")

save_fig("FigS2.tiff", list(ps2a, ps2b, ps2c), matrix(c(1, 2, 2, 3), nrow = 1), width = 7.2, height = 4.2)


tiff(paste0(OUT, "FigS2A.tiff"), width = 8, height = 5, units = "in", res = 300, compression = "lzw")
print(ps2a)
dev.off()

tiff(paste0(OUT, "FigS2B.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(ps2b)
dev.off()


# =============================================================================
# fig. S3. Grounding and tracking with external blood concentration
# =============================================================================
cat("-- fig. S3\n")

ext <- read.csv("../output/summaries/fb_external_abundance.csv", check.names = FALSE) %>%
  mutate(model = recode(model,
                        `naive (covered subset)` = "Unadjusted",
                        `adjusted: SomaScan signal + CV` = "SomaScan signal + CV",
                        `adjusted: HPA MS conc + CV` = "HPA MS concentration + CV",
                        `adjusted: HPA MS conc + SomaScan signal + CV` = "HPA MS concentration +\nSomaScan signal + CV",
                        `adjusted: HPA IM conc + CV` = "HPA immunoassay concentration + CV"),
         model = factor(model, levels = rev(model)),
         txt = sprintf("p = %s, n = %s", fmt_p_plain(p_grounding), format(n, big.mark = ",")))
ps3 <- ggplot(ext, aes(beta_grounding, model)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_point(size = 2) +
  geom_text(aes(label = txt), vjust = -1.1, size = 3.5) +
  scale_x_continuous(limits = c(-0.0035, 0.0045)) +
  labs(x = "Grounding score coefficient", y = NULL,
       subtitle = sprintf("Proteins with an HPA blood concentration estimate\nPermutation p, HPA MS model = %s",
                          fmt_p_plain(ext$permutation_p_external[!is.na(ext$permutation_p_external)][1]))) +
  base + theme(
               axis.text.y = element_text(size = 9),
               axis.title.x = element_text(size=12))

save_fig("FigS3.tiff", list(ps3), matrix(1), width = 4.8, height = 3.4)

tiff(paste0(OUT, "FigS3.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(ps3)
dev.off()



# =============================================================================
# fig. S4. Route of the intracellular progression markers, and top-ranked tissue
# =============================================================================
cat("-- fig. S4\n")

tri_fdr <- read.csv("../output/landscape/triage_progression_fdr.csv", stringsAsFactors = FALSE)
muscle <- c("CKM", "ACTN2", "MYOM2")
ic <- tri_fdr %>% filter(route == "intracellular") %>%
  mutate(class = case_when(gene_symbol %in% muscle ~ "Muscle",
                           specificity == "specific" ~ "Tissue-specific",
                           TRUE ~ "Broadly expressed"),
         class = factor(class, levels = c("Muscle", "Broadly expressed", "Tissue-specific")))
ic_lab <- ic %>% group_by(class) %>% summarise(n = n(), genes = paste(sort(gene_symbol), collapse = ", "), .groups = "drop")
ps4a <- ggplot(ic_lab, aes(class, n, fill = class)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = genes), vjust = -0.5, size = 3, lineheight = 0.9) +
  scale_fill_manual(values = c(Muscle = col_up, `Broadly expressed` = "grey55", `Tissue-specific` = "#C2A5CF"), guide = "none") +
  scale_y_continuous(limits = c(0, 9), breaks = 0:8, expand = expansion(mult = c(0, 0.02))) +
  labs(tag = "", x = NULL, y = "Intracellular progression markers",
       subtitle = sprintf("%d of %d progression markers are intracellular", nrow(ic), nrow(tri_fdr))) +
  base + theme(
              axis.text.x = element_text(size = 9),
              axis.title.y = element_text(size=12))

top <- read.csv("../output/summaries/fb_top_tissue_by_set.csv", stringsAsFactors = FALSE) %>%
  filter(set %in% c("prog_fdr", "background"), !is.na(top_tissue))
cx <- top %>% group_by(set) %>% summarise(k = sum(n[top_tissue == "Cerebral cortex"]), N = sum(n), .groups = "drop")
ft_cx <- fisher.test(matrix(c(cx$k[cx$set == "prog_fdr"], cx$N[cx$set == "prog_fdr"] - cx$k[cx$set == "prog_fdr"],
                              cx$k[cx$set == "background"], cx$N[cx$set == "background"] - cx$k[cx$set == "background"]), 2))
keep_t <- top %>% filter(set == "prog_fdr") %>% arrange(desc(n)) %>% head(8) %>% pull(top_tissue)
td <- top %>% filter(top_tissue %in% keep_t) %>%
  mutate(set = factor(ifelse(set == "background", "Panel", "Progression markers"), levels = c("Panel", "Progression markers")),
         top_tissue = factor(top_tissue, levels = rev(keep_t)))
ps4b <- ggplot(td, aes(pct, top_tissue, fill = set)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c(Panel = "grey60", `Progression markers` = col_up), name = NULL) +
  labs(tag = "", x = "Proteins for which the tissue ranks first (%)", y = NULL,
       subtitle = sprintf("Cerebral cortex: %.0f%% of progression markers, %.0f%% of panel\nOR = %.2f, p = %s",
                          100 * cx$k[cx$set == "prog_fdr"] / cx$N[cx$set == "prog_fdr"],
                          100 * cx$k[cx$set == "background"] / cx$N[cx$set == "background"],
                          ft_cx$estimate, fmt_p_plain(ft_cx$p.value))) +
  base + theme(legend.position = "bottom", 
               legend.text = element_text(size=9),
               axis.text.y = element_text(size = 9),
               axis.title.x = element_text(size=12))

save_fig("FigS4.tiff", list(ps4a, ps4b), matrix(1:2, nrow = 1), width = 7.2, height = 3.6)

tiff(paste0(OUT, "FigS4A.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(ps4a)
dev.off()

tiff(paste0(OUT, "FigS4B.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(ps4b)
dev.off()


# =============================================================================
# fig. S5. NfL slopes within thirds of baseline NfL, all draws and follow-up only
# =============================================================================
cat("-- fig. S5\n")

t_all <- read.csv("../output/summaries/fb_nfl_plateau_by_tertile.csv")
t_fu  <- read.csv("../output/summaries/fb_nfl_plateau_followup_only.csv")
z_het <- function(r1, n1, r2, n2) 2 * pnorm(-abs((atanh(r1) - atanh(r2)) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))))
p_all <- z_het(t_all$r_slope_slope[1], t_all$n[1], t_all$r_slope_slope[3], t_all$n[3])
p_fu  <- z_het(t_fu$r_slope_slope[1], t_fu$n[1], t_fu$r_slope_slope[3], t_fu$n[3])
s5 <- bind_rows(t_all %>% transmute(base_tertile, n, r = r_slope_slope, p, set = sprintf("All draws (%d patients)", sum(t_all$n))),
                t_fu %>% transmute(base_tertile, n, r = r_slope_slope, p, set = sprintf("Follow-up draws only (%d patients)", sum(t_fu$n)))) %>%
  bind_cols(fisher_ci(.$r, .$n)) %>%
  mutate(third = factor(c("Low", "Intermediate", "High")[base_tertile], levels = c("Low", "Intermediate", "High")),
         set = factor(set, levels = unique(set)))
ps5 <- ggplot(s5, aes(third, r, colour = set)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.45), size = 0.6) +
  scale_colour_manual(values = c("black", col_up), name = NULL) +
  labs(x = "Baseline serum NfL", y = "Correlation of NfL slope with ALSFRS-R slope (95% CI)",
       subtitle = sprintf("Low against high: p = %s (all draws), p = %s (follow-up only)",
                          fmt_p_plain(p_all), fmt_p_plain(p_fu))) +
  base + theme(legend.position = "bottom",
               legend.text = element_text(size = 9),
               axis.text.x = element_text(size = 9),
               axis.title.x = element_text(size=12),
               axis.title.y = element_text(size=12))

save_fig("FigS5.tiff", list(ps5), matrix(1), width = 4.2, height = 3.8)


tiff(paste0(OUT, "FigS5.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(ps5)
dev.off()

# =============================================================================
# fig. S6. Robustness of the loss of tracking in the fastest progressors
# =============================================================================
cat("-- fig. S6\n")

ct <- fc_env$cut_tab %>%
  mutate(split = factor(c(`2` = "Halves", `3` = "Thirds", `4` = "Quartiles", `5` = "Quintiles")[as.character(groups)],
                        levels = c("Halves", "Thirds", "Quartiles", "Quintiles")))
ct_long <- bind_rows(ct %>% transmute(split, group = "Fastest group", beta = beta_fastest, lo = beta_fastest_lo, hi = beta_fastest_hi),
                     ct %>% transmute(split, group = "Slowest group", beta = beta_slowest, lo = NA_real_, hi = NA_real_))
ps6a <- ggplot(ct_long, aes(split, beta, colour = group, group = group)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.4) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3, na.rm = TRUE) +
  geom_point(size = 1.6) +
  geom_text(data = ct, aes(split, 0.1, label = paste0("p = ", fmt_p_plain(p_fastest_vs_slowest))),
            inherit.aes = FALSE, size = 3) +
  scale_colour_manual(values = c(`Fastest group` = col_up, `Slowest group` = col_down), name = NULL) +
  labs(tag = "", x = "Patients divided into", y = "Composite change per ALSFRS-R point",
       subtitle = "Cross-fitted composite. p: fastest against slowest") +
  base + theme(legend.position = "bottom",
               legend.text = element_text(size = 9),
               axis.text.x = element_text(size = 9),
               axis.title.x = element_text(size=12),
               axis.title.y = element_text(size=12))

inter_p <- function(dd, ref = "Slow", hi = "Fast") {
  dd$grp <- relevel(factor(as.character(dd$tertile)), ref = ref)
  summary(lm(composite ~ als_slope * grp, data = dd))$coefficients[paste0("als_slope:grp", hi), "Pr(>|t|)"]
}
p_1y <- inter_p(fs_env$comp %>% filter(span_years >= 1))
cat(sprintf("  Interaction p, at least 1 year of follow-up: %.3f\n", p_1y))
rb <- bind_rows(
  fs_env$tert_ins %>% mutate(analysis = "All patients", p = fs_env$p_ins),
  fs_env$tert_1y %>% mutate(analysis = "Followed for at least 1 year", p = p_1y),
  fs_env$term_tab %>% filter(window_months == 6) %>% mutate(analysis = "Excluding samples within\n6 months of death", p = p_fast_vs_slow),
  fs_env$term_tab %>% filter(window_months == 12) %>% mutate(analysis = "Excluding samples within\n12 months of death", p = p_fast_vs_slow)) %>%
  filter(tertile %in% c("Fast", "Slow")) %>%
  mutate(tertile = factor(ifelse(tertile == "Fast", "Fastest third", "Slowest third"), levels = c("Fastest third", "Slowest third")),
         analysis = factor(analysis, levels = rev(unique(analysis))))
rb_lab <- rb %>% distinct(analysis, p) %>% mutate(txt = paste0("p = ", fmt_p_plain(p)))
ps6b <- ggplot(rb, aes(beta, analysis, colour = tertile)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.55), size = 0.6) +
  geom_text(data = rb_lab, aes(0.27, analysis, label = txt), inherit.aes = FALSE, size = 3, hjust = 1) +
  scale_x_continuous(limits = c(-0.03, 0.27)) +
  scale_colour_manual(values = c(`Fastest third` = col_up, `Slowest third` = col_down), name = NULL) +
  labs(tag = "", x = "Change per ALSFRS-R point (95% CI)", y = NULL,
       subtitle = "In-sample composite. p: fastest against slowest") +
  base + theme(legend.position = "bottom",
               legend.text = element_text(size = 9),
               axis.text.y = element_text(size = 9),
               axis.title.x = element_text(size=12))

save_fig("FigS6.tiff", list(ps6a, ps6b), matrix(1:2, nrow = 1), width = 7.2, height = 3.6)

tiff(paste0(OUT, "FigS6A.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(ps6a)
dev.off()

tiff(paste0(OUT, "FigS6B.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(ps6b)
dev.off()
cat("\n== SCRIPT 10 COMPLETE ==\n")
