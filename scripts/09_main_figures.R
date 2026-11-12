################################################################################
#  Main manuscript figures
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

OUT <- "../output/manuscript_figures/"

cat("================================================================\n")
cat("  09: Main manuscript figures\n")
cat("================================================================\n\n")

# ---- Shared settings ---------------------------------------------------------

base <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 7, colour = "grey25"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_rect(fill = "grey92", colour = NA))

has_repel <- requireNamespace("ggrepel", quietly = TRUE)
label_layer <- function(data, mapping, size = 2.3) {
  if (has_repel) ggrepel::geom_text_repel(data = data, mapping = mapping, size = size,
                                          min.segment.length = 0, max.overlaps = Inf,
                                          segment.colour = "grey50", segment.size = 0.25)
  else geom_text(data = data, mapping = mapping, size = size, vjust = -0.6, check_overlap = TRUE)
}

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
tier_lab <- function(x) factor(unname(tier_labels[x]), levels = unname(tier_labels))
tier_axis <- c("Expressed,\npromoter-\nsupported", "Expressed,\npromoter-\nunconfirmed", "Transcript\nonly", "Protein\nonly", "Not\nexpressed")

col_up <- "#B2182B"; col_down <- "#2166AC"; col_ns <- "grey75"
set_colors <- c(Panel = "grey60", `Differential plasma proteins` = "#2166AC", `Progression markers` = "#B2182B")

fmt_p <- function(p) {
  ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p)) %>% sub("e-0?", " x 10^-", .)
}
fmt_p_plain <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "fg", digits = 2))

fisher_ci <- function(r, n) {
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  data.frame(lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se))
}

blank_panel <- function(tag, text) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = text, size = 3, colour = "grey50") +
    xlim(0, 1) + ylim(0, 1) + labs(tag = tag) + theme_void() +
    theme(plot.tag = element_text(face = "bold", size = 11))
}

# layout: a matrix of panel indices (NA = empty)
save_fig <- function(file, plots, layout, width, height) {
  tiff(paste0(OUT, file), width = width, height = height, units = "in", res = 300, compression = "lzw")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
  for (i in seq_along(plots)) {
    cells <- which(layout == i, arr.ind = TRUE)
    print(plots[[i]], vp = viewport(layout.pos.row = range(cells[, 1])[1]:range(cells[, 1])[2],
                                    layout.pos.col = range(cells[, 2])[1]:range(cells[, 2])[2]))
  }
  invisible(dev.off())
  cat("  Saved", paste0(OUT, file), "\n")
}

# ---- Load ---------------------------------------------------------------------

gdedup <- readRDS("../output/grounding/grounding_dedup.rds")
prog <- read.csv("../output/grounding/progression_with_grounding.csv", stringsAsFactors = FALSE)
daps <- read.csv("../output/landscape/daps_full.csv", stringsAsFactors = FALSE) %>%
  left_join(gdedup %>% select(uniprot, gene_symbol, grounding_tier), by = "uniprot")
stage <- read.csv("../output/summaries/fb_dap_stage_association.csv", stringsAsFactors = FALSE)
ctx <- read.csv("../output/summaries/fb_concordr_context_presence.csv", stringsAsFactors = FALSE)
route_or <- read.csv("../output/summaries/fb_concordr_route_enrichment.csv", stringsAsFactors = FALSE)
tri_bg <- read.csv("../output/landscape/triage_background.csv", stringsAsFactors = FALSE)
tri_dap <- read.csv("../output/landscape/triage_daps.csv", stringsAsFactors = FALSE)
tri_fdr <- read.csv("../output/landscape/triage_progression_fdr.csv", stringsAsFactors = FALSE)
ent_fdr <- read.csv("../output/landscape/tissue_profile_gene_summary_fdr.csv", stringsAsFactors = FALSE)
ent_dap <- read.csv("../output/landscape/tissue_profile_gene_summary_dap.csv", stringsAsFactors = FALSE)
ent_tbl <- read.csv("../output/summaries/fb_tissue_entropy.csv", stringsAsFactors = FALSE)
tiers <- read.csv("../output/summaries/fb_grounding_tiers.csv", stringsAsFactors = FALSE)
known <- read.csv("../output/summaries/fb_known_marker_grounding.csv", stringsAsFactors = FALSE)
primary <- read.csv("../output/summaries/fb_primary_result.csv", stringsAsFactors = FALSE)
tier_enr <- read.csv("../output/summaries/fb_tier_enrichment.csv", stringsAsFactors = FALSE)
nfl_des <- read.csv("../output/summaries/fb_nfl_designs.csv", stringsAsFactors = FALSE)
nfl_ia <- read.csv("../output/summaries/fb_nfl_immunoassay.csv", stringsAsFactors = FALSE)
nfl_rel <- read.csv("../output/summaries/fb_nfl_reliability.csv", stringsAsFactors = FALSE)
nfl_conc <- read.csv("../output/summaries/fb_nfl_concordance.csv", stringsAsFactors = FALSE)
rep_dap <- read.csv("../output/summaries/fb_replication_dap_concordance.csv", stringsAsFactors = FALSE)
tmt_dap <- read.csv("../output/replication/target_als_tmt_plasma_dap.csv", stringsAsFactors = FALSE)
olk_dap <- read.csv("../output/replication/target_als_olink_dap.csv", stringsAsFactors = FALSE)
tmt_prog <- read.csv("../output/replication/target_als_tmt_plasma_progression.csv", stringsAsFactors = FALSE)

nfl_fig_file <- "../output/nfl_check/nfl_figure_data.RData"
nfl_env <- new.env(); load(nfl_fig_file, envir = nfl_env)

fs_env <- new.env(); load("../output/validation/fast_slow_tracking.RData", envir = fs_env)
fc_env <- new.env(); load("../output/validation/fast_progressor_clinical.RData", envir = fc_env)
ta_env <- new.env(); load("../output/validation/target_als_replication.RData", envir = ta_env)

# Gene-level sets as used throughout (one row per gene, aptamer with the smallest p)
dap_sig  <- daps %>% filter(sig %in% TRUE)
dap_gene <- dap_sig %>% filter(!is.na(gene_symbol)) %>% group_by(gene_symbol) %>%
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>% ungroup()
prog_sig  <- prog %>% filter(p_fdr < 0.05)
prog_gene <- prog_sig %>% filter(!is.na(gene_symbol)) %>% group_by(gene_symbol) %>%
  slice_min(p_value, n = 1, with_ties = FALSE) %>% ungroup()
cat("DAPs:", nrow(dap_sig), "aptamers,", nrow(dap_gene), "genes | progression markers:",
    nrow(prog_sig), "aptamers,", nrow(prog_gene), "genes\n\n")


# =============================================================================
# FIGURE 1. Plasma proteins that distinguish ALS from controls
# =============================================================================
cat("-- Figure 1\n")

f1_labels <- c("PDLIM3", "ACTN2", "MYOM2", "CA3", "TNNI2", "MYL11", "CKM", "MSTN",
               "ART3", "ANTXR2", "CRTAC1", "HPGDS")
vd <- daps %>% mutate(class = case_when(sig %in% TRUE & logFC > 0 ~ "Higher in ALS",
                                        sig %in% TRUE & logFC < 0 ~ "Lower in ALS",
                                        TRUE ~ "Not significant"))
lab1 <- vd %>% filter(sig %in% TRUE, gene_symbol %in% f1_labels) %>%
  group_by(gene_symbol) %>% slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>% ungroup()

p1b <- ggplot(vd, aes(logFC, -log10(P.Value), colour = class)) +
  geom_point(data = vd %>% filter(class == "Not significant"), size = 0.4, alpha = 0.5) +
  geom_point(data = vd %>% filter(class != "Not significant"), size = 1) +
  label_layer(lab1, aes(logFC, -log10(P.Value), label = gene_symbol)) +
  scale_colour_manual(values = c(`Higher in ALS` = col_up, `Lower in ALS` = col_down, `Not significant` = col_ns),
                      name = NULL) +
  labs(tag = "", x = expression(log[2]~"fold change, ALS against controls"), y = expression(-log[10]~p),
       subtitle = sprintf("%d aptamers at FDR 5%%: %d higher, %d lower in ALS",
                          sum(vd$sig %in% TRUE), sum(vd$class == "Higher in ALS"), sum(vd$class == "Lower in ALS"))) +
  base + theme(legend.position = "bottom",
               axis.title.x = element_text(size=12),
               axis.title.y = element_text(size=12),
               legend.text = element_text(size =11))

tiff(paste0(OUT, "Fig1B.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(p1b)
dev.off()


ctx_test <- function(set_a, set_b = "background") {
  d <- ctx %>% filter(set %in% c(set_a, set_b)) %>%
    mutate(n_det = round(pct_detected / 100 * n))
  d %>% group_by(context_tissue) %>% group_modify(~ {
    a <- .x[.x$set == set_a, ]; b <- .x[.x$set == set_b, ]
    ft <- fisher.test(matrix(c(a$n_det, a$n - a$n_det, b$n_det, b$n - b$n_det), 2))
    tibble(or = unname(ft$estimate), p = ft$p.value)
  }) %>% ungroup()
}
ctx_plot <- function(set_a, label_a, tag) {
  tt <- ctx_test(set_a)
  d <- ctx %>% filter(set %in% c(set_a, "background")) %>%
    mutate(set = factor(ifelse(set == "background", "Panel", label_a), levels = c("Panel", label_a)),
           context_tissue = factor(context_tissue,
                                   levels = c("Skeletal muscle", "Heart muscle", "Cerebral cortex", "Liver", "Adipose tissue")))
  lab <- tt %>% mutate(context_tissue = factor(context_tissue, levels = levels(d$context_tissue)),
                       txt = sprintf("OR %.2f\np = %s", or, fmt_p_plain(p)))
  set_colors1 <- c(Panel = "grey60", `Differential plasma proteins` = "#660066", `Progression markers` = "#B2182B")
  ggplot(d, aes(context_tissue, pct_detected, fill = set)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(data = lab, aes(context_tissue, 92, label = txt), inherit.aes = FALSE, size = 3, lineheight = 0.9) +
    scale_fill_manual(values = set_colors1, name = NULL) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
    labs(tag = tag, x = NULL, y = "Detected in tissue (% of proteins)") +
    base + theme(legend.position = "bottom", 
                 axis.text.x = element_text(angle = 25, hjust = 1, size=11),
                 axis.title.x = element_text(size=12),
                 axis.title.y = element_text(size=12),
                 legend.text = element_text(size =11))
}
p1c <- ctx_plot("dap", "Differential plasma proteins", "")

tiff(paste0(OUT, "Fig1C.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(p1c)
dev.off()


st <- stage %>% filter(covariate == "alsfrs_base") %>%
  mutate(sig = ifelse(p_fdr < 0.05, "FDR < 0.05", "Not significant"),
         dirn = ifelse(dap_logFC > 0, "Higher in ALS", "Lower in ALS"))
lab1d <- st %>% filter(gene_symbol %in% c("ANTXR2", "CLEC3B", "MXRA8", "CRTAC1", "EHMT2", "COL28A1",
                                          "TNNT2", "DCLK1", "CALCB", "PDLIM3", "CA3", "PPP1R10", "THBS4"))

p1d <- ggplot(st, aes(dap_logFC, r)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(aes(colour = dirn, alpha = sig), size = 1.4) +
  label_layer(lab1d, aes(dap_logFC, r, label = gene_symbol)) +
  scale_colour_manual(values = c(`Higher in ALS` = col_up, `Lower in ALS` = col_down), name = NULL) +
  scale_alpha_manual(values = c(`FDR < 0.05` = 1, `Not significant` = 0.3), name = NULL) +
  labs(tag = "", x = expression(log[2]~"fold change, ALS against controls"),
       y = "Correlation with baseline ALSFRS-R",
       subtitle = sprintf("%d of %d differential plasma proteins at FDR 5%% (n = %d patients)",
                          sum(st$p_fdr < 0.05), nrow(st), 244)) +
  base + theme(legend.position = "bottom",
               axis.title.x = element_text(size=12),
               axis.title.y = element_text(size=12),
               legend.text = element_text(size =11))

tiff(paste0(OUT, "Fig1D.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(p1d)
dev.off()


save_fig("Fig1.tiff", list(blank_panel("A", "Study design (drawn separately)"), p1b, p1c, p1d),
         matrix(c(1, 1, 2, 2, 2, 2, 3, 4, 3, 4), nrow = 5, byrow = TRUE), width = 7.2, height = 9)


# =============================================================================
# FIGURE 2. Progression markers
# =============================================================================
cat("-- Figure 2\n")

f2_labels <- c("CLEC3B", "ANTXR2", "ART3", "HPGDS", "MYOC", "CRTAC1", "IGFBP6", "EHMT2",
               "CXCL12", "ANOS1", "TNNT2", "FTL", "C9", "HAMP", "S100A12")
pv <- prog %>% mutate(class = ifelse(p_fdr < 0.05, ifelse(cor_r > 0, "Falls with faster decline", "Rises with faster decline"),
                                     "Not significant"))
lab2 <- pv %>% filter(p_fdr < 0.05, gene_symbol %in% f2_labels) %>% group_by(gene_symbol) %>%
  slice_min(p_value, n = 1, with_ties = FALSE) %>% ungroup()
fdr_line <- max(pv$p_value[pv$p_fdr < 0.05])
p2a <- ggplot(pv, aes(cor_r, -log10(p_value), colour = class)) +
  geom_hline(yintercept = -log10(fdr_line), linetype = 2, colour = "grey50") +
  geom_point(data = pv %>% filter(class == "Not significant"), size = 0.4, alpha = 0.5) +
  geom_point(data = pv %>% filter(class != "Not significant"), size = 1) +
  label_layer(lab2, aes(cor_r, -log10(p_value), label = gene_symbol)) +
  scale_colour_manual(values = c(`Falls with faster decline` = col_down, `Rises with faster decline` = col_up,
                                 `Not significant` = col_ns), name = NULL) +
  labs(tag = "", x = "Correlation of protein slope with ALSFRS-R slope", y = expression(-log[10]~p),
       subtitle = sprintf("229 patients | %d aptamers nominal, %d at FDR 5%% (dashed line)",
                          sum(pv$p_value < 0.05), sum(pv$p_fdr < 0.05))) +
  base + theme(legend.position = "bottom",
               axis.title.x = element_text(size=12),
               axis.title.y = element_text(size=12),
               legend.text = element_text(size =11))

tiff(paste0(OUT, "Fig2A.tiff"), width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(p2a)
dev.off()


ov <- dap_gene %>% mutate(dirn = ifelse(logFC > 0, "Higher in ALS", "Lower in ALS"),
                          is_prog = gene_symbol %in% prog_gene$gene_symbol)
ov_tab <- table(ov$dirn, ov$is_prog)
ov_ft <- fisher.test(ov_tab[c("Lower in ALS", "Higher in ALS"), c("TRUE", "FALSE")])
ov_or <- (ov_tab["Lower in ALS", "TRUE"] * ov_tab["Higher in ALS", "FALSE"]) /
  (ov_tab["Lower in ALS", "FALSE"] * ov_tab["Higher in ALS", "TRUE"])
cat(sprintf("  Overlap: lower %d of %d, higher %d of %d | OR = %.2f, p = %.3f\n",
            ov_tab["Lower in ALS", "TRUE"], sum(ov_tab["Lower in ALS", ]),
            ov_tab["Higher in ALS", "TRUE"], sum(ov_tab["Higher in ALS", ]), ov_or, ov_ft$p.value))
ovd <- ov %>% group_by(dirn) %>% summarise(n = n(), k = sum(is_prog), .groups = "drop") %>%
  mutate(pct = 100 * k / n, txt = sprintf("%d of %d", k, n))
p2b <- ggplot(ovd, aes(dirn, pct, fill = dirn)) +
  geom_col(width = 0.6) + geom_text(aes(label = txt), vjust = -0.4, size = 2.6) +
  scale_fill_manual(values = c(`Higher in ALS` = col_up, `Lower in ALS` = col_down), guide = "none") +
  scale_y_continuous(limits = c(0, 70), expand = expansion(mult = c(0, 0.02))) +
  labs(tag = "", x = "Differential plasma proteins", y = "Progression markers (%)",
       subtitle = sprintf("OR %.2f, Fisher p = %s", ov_or, fmt_p_plain(ov_ft$p.value))) +
  base + theme(
               axis.title.x = element_text(size=12),
               axis.text = element_text(size =11),
               axis.title.y = element_text(size=12))


tiff(paste0(OUT, "Fig2B.tiff"), width = 5, height = 5, units = "in", res = 300, compression = "lzw")
print(p2b)
dev.off()


save_fig("Fig2.tiff", list(p2a, p2b), matrix(c(1, 1, 2), nrow = 1), width = 7.2, height = 3.6)


# =============================================================================
# FIGURE 3. Replication in Target ALS
# =============================================================================
cat("-- Figure 3\n")

aals_all <- daps %>% filter(!is.na(gene_symbol)) %>% group_by(gene_symbol) %>%
  summarise(aals_logFC = logFC[which.min(adj.P.Val)], .groups = "drop")
bg_rho <- function(res) {
  d <- aals_all %>% inner_join(res %>% select(gene_symbol, est = logFC), by = "gene_symbol")
  c(rho = unname(cor(d$aals_logFC, d$est, method = "pearson", use = "complete.obs")), n = nrow(d))
}
rep_plot <- function(platform, res, tag, ylab) {
  d <- rep_dap %>% filter(platform == !!platform)
  same <- sign(d$aals_logFC) == sign(d$est)
  bt <- binom.test(sum(same), nrow(d), 0.5)
  rho <- cor(d$aals_logFC, d$est, method = "pearson")
  bg <- bg_rho(res)
  lab <- if (nrow(d) <= 15) d else d %>% filter(adj.P.Val < 0.05 | abs(aals_logFC) > 0.8)
  ggplot(d, aes(aals_logFC, est)) +
    geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
    geom_point(aes(colour = adj.P.Val < 0.05), size = 1.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey40", linewidth = 0.4) +
    label_layer(lab, aes(aals_logFC, est, label = gene_symbol), size = 2.1) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"),
                        labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "NS"), name = NULL) +
    labs(tag = tag, x = expression("Answer ALS "*log[2]*"FC (SomaScan)"), y = ylab,
         subtitle = sprintf("%d of %d same direction, p = %s\nr = %.2f (all shared: r = %.2f, n = %d)",
                            sum(same), nrow(d), fmt_p_plain(bt$p.value), rho, bg["rho"], bg["n"])) +
    base + theme(legend.position = "bottom",
                 axis.title.x = element_text(size=12),
                 axis.title.y = element_text(size=12),
                 legend.text = element_text(size =11))
}
p3a <- rep_plot("TMT plasma", tmt_dap, "", "Target ALS fold change (TMT)")
p3b <- rep_plot("Olink HT", olk_dap, "", "Target ALS fold change (Olink)")


tiff(paste0(OUT, "Fig3A.tiff"), width = 5, height = 5, units = "in", res = 300, compression = "lzw")
print(p3a)
dev.off()

tiff(paste0(OUT, "Fig3B.tiff"), width = 5, height = 5, units = "in", res = 300, compression = "lzw")
print(p3b)
dev.off()


aals_prog <- prog_gene %>% transmute(gene_symbol, aals_cor_r = cor_r)
pc <- aals_prog %>% inner_join(tmt_prog %>% select(gene_symbol, ta_r = cor_r, ta_p = p_value), by = "gene_symbol")
same_p <- sign(pc$aals_cor_r) == sign(pc$ta_r)
bt_p <- binom.test(sum(same_p), nrow(pc), 0.5)
rho_p <- cor.test(pc$aals_cor_r, pc$ta_r, method = "pearson")
bgp <- prog %>% filter(!is.na(gene_symbol)) %>% group_by(gene_symbol) %>%
  summarise(a = cor_r[which.min(p_value)], .groups = "drop") %>%
  inner_join(tmt_prog %>% select(gene_symbol, b = cor_r), by = "gene_symbol")
rho_bg <- cor(bgp$a, bgp$b, method = "pearson", use = "complete.obs")
cat(sprintf("  Progression replication: %d present, %d same direction, r = %.2f, background r = %.2f (n = %d)\n",
            nrow(pc), sum(same_p), rho_p$estimate, rho_bg, nrow(bgp)))
p3c <- ggplot(pc, aes(aals_cor_r, ta_r)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(aes(colour = ta_p < 0.05), size = 1.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey40", linewidth = 0.4) +
  label_layer(pc, aes(aals_cor_r, ta_r, label = gene_symbol), size = 2.1) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"),
                      labels = c(`TRUE` = "p < 0.05", `FALSE` = "NS"), name = NULL) +
  labs(tag = "C", x = "Answer ALS r with ALSFRS-R slope",
       y = "Target ALS r with ALSFRS-R slope",
       subtitle = sprintf("%d of %d same direction, p = %s\nr = %.2f (all shared: r = %.2f, n = %d)",
                          sum(same_p), nrow(pc), fmt_p_plain(bt_p$p.value), rho_p$estimate, rho_bg, nrow(bgp))) +
  base + theme(legend.position = "bottom",
               axis.title.x = element_text(size=12),
               axis.title.y = element_text(size=12),
               legend.text = element_text(size =11))


tiff(paste0(OUT, "Fig3C.tiff"), width = 5, height = 5, units = "in", res = 300, compression = "lzw")
print(p3c)
dev.off()



tc <- ta_env$comp
ci_r <- quantile(ta_env$r_boot, c(0.025, 0.975), na.rm = TRUE)
p3d <- ggplot(tc, aes(als_slope, composite)) +
  geom_point(size = 1.4, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, colour = col_down, fill = col_down, alpha = 0.15, linewidth = 0.6) +
  labs(tag = "", x = "ALSFRS-R slope (points per visit)", y = "Composite change (z)",
       subtitle = sprintf("%d patients, %d markers\nr = %.2f (95%% CI %.2f to %.2f)\np = %s against shuffled ALSFRS-R slopes",
                                     nrow(tc), nrow(ta_env$mk), ta_env$r_obs, ci_r[1], ci_r[2], fmt_p_plain(ta_env$p_perm))) +
  base+theme(
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12),
            legend.text = element_text(size =11))

p3e <- ggplot(data.frame(r = ta_env$r_null), aes(r)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "grey45", linewidth = 0.2) +
  geom_vline(xintercept = ta_env$r_obs, colour = col_up, linewidth = 0.8) +
  labs(tag = "", x = "Composite correlation with ALSFRS-R slope", y = "Random sets",
       subtitle = sprintf("1,000 random sets of %d non-marker proteins\nobserved r = %.2f (red)\np = %s against random protein sets",
                          nrow(ta_env$mk), ta_env$r_obs, fmt_p_plain(ta_env$p_null))) +
  base + theme(
               axis.title.x = element_text(size=12),
               axis.title.y = element_text(size=12),
               legend.text = element_text(size =11))

save_fig("Fig3.tiff", list(p3a, p3b, p3c, p3d, p3e),
         matrix(c(1, 2, 3, 4, 5, NA), nrow = 2, byrow = TRUE), width = 7.2, height = 5.8)

tiff(paste0(OUT, "Fig3D.tiff"), width = 5, height = 5, units = "in", res = 300, compression = "lzw")
print(p3d)
dev.off()

tiff(paste0(OUT, "Fig3E.tiff"), width = 5, height = 5, units = "in", res = 300, compression = "lzw")
print(p3e)
dev.off()



# =============================================================================
# FIGURE 4. Motor neuron expression map
# =============================================================================
cat("-- Figure 4\n")

td <- tiers %>% mutate(tier = tier_lab(grounding_tier))
p4b <- ggplot(td, aes(tier, n, fill = tier)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(n, big.mark = ","), pct)), vjust = -0.2, size = 3, lineheight = 0.9) +
  scale_fill_manual(values = tier_colors, guide = "none") +
  scale_y_continuous(limits = c(0, max(td$n) * 1.2), expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(labels = tier_axis) +
  labs(tag = "", x = NULL, y = "Plasma proteins (aptamers)",
       subtitle = sprintf("%s aptamers mapped to an annotated gene", format(sum(td$n), big.mark = ","))) +
  base+theme(axis.text.x = element_text(size=9),
             axis.title.y = element_text(size=12))


tiff(paste0(OUT, "Fig4B.tiff"), width = 5, height = 5, units = "in", res = 300, compression = "lzw")
print(p4b)
dev.off()


# ---- Fig 4C. Known motor neuron and glial markers, per iMN line ----
# Same samples, normalisation and detection rules as 01_grounding_map.R
mn_genes <- c("NEFL", "NEFH", "STMN2", "TARDBP", "SOD1", "UCHL1")
glia_genes <- c("GFAP", "CHIT1", "CHI3L1", "CHI3L2", "MMP9")
grp_lab <- c("Motor neuron", "Glial and inflammatory")

ids4 <- read.csv("../output/grounding/grounding_map.csv", stringsAsFactors = FALSE) %>%
  filter(gene_symbol %in% c(mn_genes, glia_genes)) %>% distinct(gene_symbol, ensembl_id, uniprot)
mapper4 <- read.csv("../raw/transcriptomics/4_matrix/Sample Mapping Information/Sample Mapping File Feb 2024.csv",
                    check.names = FALSE)
meta4 <- read.csv("../proc/meta.csv") %>% filter(SUBJECT_GROUP %in% c("ALS", "Healthy Control"))

# transcript: VST across iMN lines, plus % of lines with > 5 reads
cnt4 <- read.csv("../raw/transcriptomics/4_matrix/AnswerALS-290-T-v1-release5_rna-counts-matrix.csv", check.names = FALSE)
rownames(cnt4) <- cnt4[, 1]; cnt4 <- cnt4[, -1]
meta_rna4 <- meta4 %>% inner_join(mapper4[, c("Participant_ID", "Sample_ID", "Batch")], by = "Participant_ID") %>%
  filter(Sample_ID %in% colnames(cnt4))
cnt4 <- cnt4[, meta_rna4$Sample_ID]
cnt4 <- cnt4[rowSums(cnt4) > 0, ]
dds4 <- DESeq2::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(cnt4, meta_rna4, design = ~ 1))
vst4 <- SummarizedExperiment::assay(DESeq2::vst(dds4, blind = TRUE))
rna_pct <- data.frame(ensembl_id = intersect(ids4$ensembl_id, rownames(cnt4))) %>%
  mutate(pct = 100 * rowMeans(cnt4[ensembl_id, , drop = FALSE] > 5)) %>%
  left_join(ids4 %>% distinct(ensembl_id, gene_symbol), by = "ensembl_id")
rna4 <- as.data.frame(t(vst4[intersect(ids4$ensembl_id, rownames(vst4)), , drop = FALSE])) %>%
  mutate(line = rownames(.)) %>% pivot_longer(-line, names_to = "ensembl_id", values_to = "value") %>%
  left_join(ids4 %>% distinct(ensembl_id, gene_symbol), by = "ensembl_id") %>%
  mutate(layer = "Transcript")
rm(cnt4, dds4, vst4)

# protein: non-imputed iMN matrix, lines restricted to those in the imputed matrix, zero = not quantified
ms4 <- read.csv("../raw/proteomics/4_matrix/AnswerALS-447-P_proteomics-protein-matrix_prebatchcorrected.txt",
                sep = "\t", check.names = FALSE)
ms4 <- ms4[, !colnames(ms4) %in% c("nFragment", "nPeptide")]
ms4$uniprot <- vapply(strsplit(as.character(ms4$Protein), "[|]"),
                      function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
ms4 <- ms4[!is.na(ms4$uniprot) & !duplicated(ms4$uniprot), ]
scols4 <- setdiff(colnames(ms4), c("Protein", "uniprot"))
sclean4 <- gsub("_", "-", gsub("(-P).+$", "\\1", scols4))
imp_cols4 <- names(read.delim("../raw/proteomics/4_matrix/AnswerALS-436-P_proteomics-protein-matrix_correctedImputed.txt",
                              nrows = 1, check.names = FALSE))
imp_clean4 <- gsub("_", "-", gsub("(-P).+$", "\\1", setdiff(imp_cols4, c("Protein", "nFragment", "nPeptide"))))
ms_ids4 <- meta4 %>% inner_join(mapper4[, c("Participant_ID", "Sample_ID")], by = "Participant_ID") %>%
  filter(Sample_ID %in% imp_clean4) %>% pull(Sample_ID)
keep4 <- scols4[sclean4 %in% ms_ids4]
prot4 <- ms4 %>% filter(uniprot %in% ids4$uniprot) %>% select(uniprot, all_of(keep4)) %>%
  pivot_longer(-uniprot, names_to = "line", values_to = "value") %>%
  mutate(value = suppressWarnings(as.numeric(value))) %>%
  left_join(ids4 %>% distinct(uniprot, gene_symbol), by = "uniprot") %>%
  mutate(layer = "Protein")
n_ms_lines <- length(keep4)
rm(ms4)
prot_pct <- prot4 %>% group_by(gene_symbol) %>%
  summarise(pct = 100 * mean(is.finite(value) & value > 0), .groups = "drop")
prot4 <- prot4 %>% filter(is.finite(value), value > 0) %>% mutate(value = log2(value))

lev4 <- c(mn_genes, glia_genes)
layer_lab <- c(Transcript = "Transcript (VST)", Protein = "Protein (log2 intensity)")
fd <- bind_rows(rna4 %>% select(gene_symbol, layer, value), prot4 %>% select(gene_symbol, layer, value)) %>%
  mutate(group = factor(ifelse(gene_symbol %in% mn_genes, grp_lab[1], grp_lab[2]), levels = grp_lab),
         gene_symbol = factor(gene_symbol, levels = lev4),
         layer = factor(unname(layer_lab[layer]), levels = unname(layer_lab)))

# annotations: % of lines above each box, one label in the empty glial protein panel
rna_y <- range(fd$value[fd$layer == layer_lab["Transcript"]])
prot_y <- range(fd$value[fd$layer == layer_lab["Protein"]])
rpct_d <- rna_pct %>%
  mutate(group = factor(ifelse(gene_symbol %in% mn_genes, grp_lab[1], grp_lab[2]), levels = grp_lab),
         gene_symbol = factor(gene_symbol, levels = lev4),
         layer = factor(unname(layer_lab["Transcript"]), levels = unname(layer_lab)),
         y = rna_y[2] + 0.06 * diff(rna_y), lab = sprintf("%.0f%%", pct))
pct_d <- prot_pct %>%
  mutate(group = factor(grp_lab[1], levels = grp_lab), gene_symbol = factor(gene_symbol, levels = lev4),
         layer = factor(unname(layer_lab["Protein"]), levels = unname(layer_lab)),
         y = prot_y[2] + 0.06 * diff(prot_y), lab = sprintf("%.0f%%", pct))
nm_d <- data.frame(gene_symbol = factor("CHI3L1", levels = lev4), group = factor(grp_lab[2], levels = grp_lab),
                   layer = factor(unname(layer_lab["Protein"]), levels = unname(layer_lab)), y = mean(prot_y))

p4c <- ggplot(fd, aes(gene_symbol, value)) +
  geom_jitter(aes(colour = group), width = 0.15, height = 0, size = 0.25, alpha = 0.35) +
  geom_boxplot(aes(fill = group), width = 0.6, outlier.shape = NA, alpha = 0.55, linewidth = 0.3) +
  geom_text(data = rpct_d, aes(gene_symbol, y, label = lab), size = 2, colour = "grey25") +
  geom_text(data = pct_d, aes(gene_symbol, y, label = lab), size = 2, colour = "grey25") +
  geom_text(data = nm_d, aes(gene_symbol, y), label = "",
            size = 2.4, colour = "grey45", lineheight = 0.85) +
  facet_grid(layer ~ group, scales = "free", space = "free_x", switch = "y") +
  scale_fill_manual(values = c(`Motor neuron` = col_down, `Glial and inflammatory` = "grey55"), guide = "none") +
  scale_colour_manual(values = c(`Motor neuron` = col_down, `Glial and inflammatory` = "grey40"), guide = "none") +
  labs(tag = "", x = NULL, y = NULL) +
  base + theme(strip.placement = "outside",
               axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 9),
               strip.text = element_text(size = 9))

cat(sprintf("  Fig 4C: %d iMN lines with RNA-seq, %d with mass spectrometry\n", length(unique(rna4$line)), n_ms_lines))

tiff(paste0(OUT, "Fig4C.tiff"), width = 5, height = 5, units = "in", res = 300, compression = "lzw")
print(p4c)
dev.off()

save_fig("Fig4.tiff", list(blank_panel("A", "Map construction (drawn separately)"), p4b, p4c),
         matrix(c(1, 1, 2, 3, 2, 3), nrow = 3, byrow = TRUE), width = 7.2, height = 6.5)


# =============================================================================
# FIGURE 5. Motor neuron expression and progression tracking
# =============================================================================
cat("-- Figure 5\n")

mg <- prog %>% filter(grounding_tier %in% tier_levels) %>% mutate(tier = tier_lab(grounding_tier), abs_r = abs(cor_r))
mean_tab <- mg %>% group_by(tier) %>% summarise(m = mean(abs_r), n = n(), .groups = "drop")
p5a <- ggplot(mg, aes(tier, abs_r, fill = tier)) +
  geom_violin(colour = NA, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_point(data = mean_tab, aes(tier, m), inherit.aes = FALSE, shape = 23, size = 1.8, fill = "black") +
  geom_text(data = mean_tab, aes(tier, 0.52, label = sprintf("mean %.3f\nn = %s", m, format(n, big.mark = ","))),
            inherit.aes = FALSE, size = 3, lineheight = 0.9) +
  scale_fill_manual(values = tier_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 0.56)) +
  scale_x_discrete(labels = tier_axis) +
  labs(tag = "", x = NULL, y = "|Correlation with ALSFRS-R slope|") +
  base + theme(axis.text.x = element_text(size=9),
               axis.title.y = element_text(size=12))

# permutation null of the adjusted grounding coefficient, same model and seed as script 02
td5 <- prog %>% filter(!is.na(grounding_score))
m_naive <- lm(abs(cor_r) ~ grounding_score, data = td5)
m_adj <- lm(abs(cor_r) ~ grounding_score + signal_mean + signal_cv, data = td5)
obs_beta <- coef(m_adj)["grounding_score"]
set.seed(42)
perm_betas <- replicate(2000, {
  d <- td5; d$grounding_score <- sample(d$grounding_score)
  coef(lm(abs(cor_r) ~ grounding_score + signal_mean + signal_cv, data = d))["grounding_score"]
})
p_perm <- mean(abs(perm_betas) >= abs(obs_beta))
cat(sprintf("  Grounding: naive beta %.4f (p = %.2g) | adjusted beta %.4f (p = %.2g) | permutation p = %.3f\n",
            coef(m_naive)[2], summary(m_naive)$coefficients[2, 4], obs_beta,
            summary(m_adj)$coefficients[2, 4], p_perm))
p5b <- ggplot(data.frame(beta = perm_betas), aes(beta)) +
  geom_histogram(bins = 50, fill = "grey75", colour = "grey45", linewidth = 0.2) +
  geom_vline(xintercept = obs_beta, colour = col_up, linewidth = 0.8) +
  labs(tag = "B", x = "Grounding score coefficient", y = "Permutations",
       subtitle = sprintf("Unadjusted: beta = %.4f, p = %s\nAdjusted for signal and CV (red line): beta = %.4f, p = %s\nPermutation null of the adjusted coefficient, p = %s",
                          coef(m_naive)[2], fmt_p_plain(summary(m_naive)$coefficients[2, 4]),
                          obs_beta, fmt_p_plain(summary(m_adj)$coefficients[2, 4]), fmt_p_plain(p_perm))) +
  base + theme(axis.text.x = element_text(size=9),
               axis.title.y = element_text(size=12),
               axis.title.x = element_text(size=12))

panel_u <- gdedup %>% filter(grounding_tier %in% tier_levels)
class_bar <- function(hit_uniprots, label, tag, enr_row) {
  d <- bind_rows(panel_u %>% mutate(set = "Panel"),
                 panel_u %>% filter(uniprot %in% hit_uniprots) %>% mutate(set = label)) %>%
    mutate(set = factor(set, levels = c("Panel", label)), tier = tier_lab(grounding_tier)) %>%
    count(set, tier) %>% group_by(set) %>% mutate(pct = 100 * n / sum(n), N = sum(n)) %>% ungroup()
  nlab <- d %>% distinct(set, N)
  ggplot(d, aes(set, pct, fill = tier)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
    geom_text(data = nlab, aes(set, 103, label = paste0("n = ", format(N, big.mark = ","))), inherit.aes = FALSE, size = 3) +
    scale_fill_manual(values = tier_colors, name = NULL) +
    scale_y_continuous(limits = c(0, 106), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(labels = function(x) gsub("Differential plasma proteins", "Differential\nplasma proteins", gsub("Progression markers", "Progression\nmarkers", x))) +
    labs(tag = tag, x = NULL, y = "Plasma proteins (%)",
         subtitle = sprintf("Expressed against all other classes:\nOR %.2f (95%% CI %.2f to %.2f), p = %s\nAll five classes: Fisher p = %s",
                            enr_row$confirmed_or, enr_row$confirmed_ci_lo, enr_row$confirmed_ci_hi,
                            fmt_p_plain(enr_row$confirmed_p), fmt_p_plain(enr_row$fisher_5tier_p))) +
    base + theme(legend.position = "right",
                 legend.text = element_text(size=9),
                 axis.text.x = element_text(size=9),
                 axis.title.y = element_text(size=12),
                 axis.title.x = element_text(size=12))
}
p5c <- class_bar(unique(dap_sig$uniprot), "Differential plasma proteins", "", tier_enr %>% filter(set == "DAPs"))
p5d <- class_bar(unique(prog_sig$uniprot), "Progression markers", "", tier_enr %>% filter(set == "Progression FDR"))

rk <- prog_sig %>% arrange(desc(abs(cor_r))) %>% mutate(rank = row_number(),
                                                         tier = factor(ifelse(grounding_tier %in% tier_levels, unname(tier_labels[grounding_tier]), "Not annotated gene"),
                                                                       levels = c(unname(tier_labels), "Not annotated gene")),
                                                         lab = ifelse(is.na(gene_symbol), SeqId, gene_symbol))
lab5e <- rk
p5e <- ggplot(rk, aes(rank, cor_r, fill = tier)) +
  geom_col(width = 0.8, colour = "grey40", linewidth = 0.1) +
  geom_text(data = lab5e, aes(rank, cor_r + sign(cor_r) * 0.02, label = lab, hjust = ifelse(cor_r > 0, 0, 1)),
            angle = 90, size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = c(tier_colors, `Not annotated gene` = "white"), name = NULL, drop = TRUE) +
  scale_y_continuous(limits = c(-0.45, 0.85)) +
  labs(tag = "", x = "Progression markers ranked by strength of association",
       y = "Correlation with ALSFRS-R slope", subtitle = sprintf("%d progression markers (aptamers)", nrow(rk))) +
  base + theme(legend.position = "bottom",
               legend.text = element_text(size=9),
               axis.title.y = element_text(size=12),
               axis.title.x = element_text(size=12))

save_fig("Fig5.tiff", list(p5a, p5b, p5c, p5d, p5e),
         matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), width = 7.2, height = 8.6)


tiff(paste0(OUT, "Fig5A.tiff"), width = 5, height = 5, units = "in", res = 300, compression = "lzw")
print(p5a)
dev.off()

tiff(paste0(OUT, "Fig5B.tiff"), width = 5, height = 5, units = "in", res = 300, compression = "lzw")
print(p5b)
dev.off()

tiff(paste0(OUT, "Fig5C.tiff"), width = 5, height = 6, units = "in", res = 300, compression = "lzw")
print(p5c)
dev.off()

tiff(paste0(OUT, "Fig5D.tiff"), width = 5, height = 6, units = "in", res = 300, compression = "lzw")
print(p5d)
dev.off()


tiff(paste0(OUT, "Fig5E.tiff"), width = 8, height = 6, units = "in", res = 300, compression = "lzw")
print(p5e)
dev.off()


# =============================================================================
# FIGURE 6. Route to plasma and tissue distribution
# =============================================================================
cat("-- Figure 6\n")

route_levels <- c("secreted", "membrane", "intracellular", "unknown")
route_cols <- c(Secreted = "#1B7837", Membrane = "#A6DBA0", Intracellular = "#C2A5CF", Unknown = "grey85")
rt <- bind_rows(tri_bg %>% mutate(set = "Panel"),
                tri_dap %>% mutate(set = "Differential plasma proteins"),
                tri_fdr %>% mutate(set = "Progression markers")) %>%
  mutate(set = factor(set, levels = names(set_colors)),
         route = factor(tools::toTitleCase(route), levels = tools::toTitleCase(route_levels))) %>%
  count(set, route) %>% group_by(set) %>% mutate(pct = 100 * n / sum(n), N = sum(n)) %>% ungroup()
or_prog <- route_or %>% filter(set == "Progression FDR")
or_dap  <- route_or %>% filter(set == "DAPs")
p6a <- ggplot(rt, aes(set, pct, fill = route)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
  geom_text(data = distinct(rt, set, N), aes(set, 103, label = paste0("n = ", format(N, big.mark = ","))),
            inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = route_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 106), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(labels = c("Panel", "Differential\nplasma proteins", "Progression\nmarkers")) +
  labs(tag = "", x = NULL, y = "Proteins (%)",
       subtitle = sprintf("Intracellular against secreted or membrane, versus panel:\nprogression markers OR %.2f (95%% CI %.2f to %.2f), p = %s\ndifferential proteins OR %.2f, p = %s",
                          or_prog$or, or_prog$ci_lo, or_prog$ci_hi, fmt_p_plain(or_prog$p), or_dap$or, fmt_p_plain(or_dap$p))) +
  base + theme(legend.position = "right",
               legend.text = element_text(size=8),
               axis.title.y = element_text(size=12),
               axis.text = element_text(size=9),
               axis.title.x = element_text(size=12))

expr_prog <- prog_gene %>% filter(grounding_tier %in% c("confirmed_full", "confirmed_partial")) %>%
  left_join(tri_fdr %>% select(gene_symbol, route), by = "gene_symbol") %>%
  mutate(route = factor(tools::toTitleCase(coalesce(route, "unknown")), levels = tools::toTitleCase(route_levels)),
         gene_symbol = reorder(gene_symbol, cor_r))
cat("  Expressed progression markers:", nrow(expr_prog), "\n")
p6b <- ggplot(expr_prog, aes(cor_r, gene_symbol, fill = route)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey50") +
  scale_fill_manual(values = route_cols, name = NULL, drop = TRUE) +
  labs(tag = "", x = "Correlation with ALSFRS-R slope", y = NULL,
       subtitle = "Progression markers expressed in iMNs") +
  base + theme(axis.text.y = element_text(face = "italic"), 
               legend.position = "right",
               legend.text = element_text(size=9),
               axis.title.x = element_text(size=12))

p6c <- ctx_plot("prog_fdr", "Progression markers", "")

ent <- bind_rows(ent_fdr %>% mutate(set = "Progression markers"), ent_dap %>% mutate(set = "Differential plasma proteins")) %>%
  mutate(set = factor(set, levels = c("Differential plasma proteins", "Progression markers")))
bg_med <- ent_tbl$entropy_norm_median[ent_tbl$set == "background"]
ent_lab <- ent_tbl %>% filter(set %in% c("dap", "prog_fdr")) %>%
  mutate(set = factor(ifelse(set == "dap", "Differential plasma proteins", "Progression markers"), levels = levels(ent$set)),
         txt = sprintf("median %.2f\nn = %d\nWilcoxon p = %s", entropy_norm_median, n_genes, fmt_p_plain(wilcox_vs_bg_norm)))
p6d <- ggplot(ent, aes(set, entropy_norm, fill = set)) +
  geom_hline(yintercept = bg_med, linetype = 2, colour = "grey40") +
  geom_violin(colour = NA, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_text(data = ent_lab, aes(set, 0.08, label = txt), inherit.aes = FALSE, size = 3, lineheight = 0.9) +
  scale_fill_manual(values = set_colors, guide = "none") +
  scale_x_discrete(labels = c("Differential\nplasma proteins", "Progression\nmarkers")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(tag = "", x = NULL, y = "Normalized tissue entropy",
       subtitle = sprintf("Dashed line: panel median %.2f (n = %s genes)",
                          bg_med, format(ent_tbl$n_genes[ent_tbl$set == "background"], big.mark = ","))) +
  base + theme(
    axis.text.x = element_text(size = 9),
    axis.title.y = element_text(size = 12) 
  )

save_fig("Fig6.tiff", list(p6a, p6b, p6c, p6d), matrix(1:4, nrow = 2, byrow = TRUE), width = 7.2, height = 6.6)

tiff(paste0(OUT, "Fig6A.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p6a)
dev.off()

tiff(paste0(OUT, "Fig6B.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p6b)
dev.off()

tiff(paste0(OUT, "Fig6C.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p6c)
dev.off()

tiff(paste0(OUT, "Fig6D.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p6d)
dev.off()

# =============================================================================
# FIGURE 7. NEFL aptamer and serum NfL
# =============================================================================
cat("-- Figure 7\n")

dz <- nfl_des %>% filter(gene %in% c("NEFL", "NEFH", "CLEC3B", "ANTXR2")) %>%
  transmute(gene, r_slope = slope_r, p_slope = slope_p, n_slope, r_base = base_r, p_base = base_p, n_base) %>%
  pivot_longer(-gene, names_to = c(".value", "design"), names_sep = "_") %>%
  bind_cols(fisher_ci(.$r, .$n)) %>%
  mutate(design = factor(ifelse(design == "slope", "Change over time (slope design)", "Baseline level (baseline design)"),
                         levels = c("Change over time (slope design)", "Baseline level (baseline design)")),
         gene = factor(gene, levels = rev(c("NEFL", "NEFH", "CLEC3B", "ANTXR2"))))
p7a <- ggplot(dz, aes(r, gene, colour = design)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5), size = 0.4) +
  scale_colour_manual(values = c("black", "#2166AC"), name = NULL) +
  labs(tag = "", x = "Correlation with ALSFRS-R slope (95% CI)", y = "SomaScan aptamer",
       subtitle = sprintf("n = %d (slope), %d (baseline)", dz$n[1], dz$n[2])) +
  base + theme(legend.position = "bottom", 
               axis.text.y = element_text(face = "italic", size=9),
               axis.title.x = element_text(size = 12),
               axis.title.y = element_text(size = 12),
               legend.text = element_text(size = 9))

mt <- nfl_env$matched
p7b <- ggplot(mt, aes(value, apt)) +
  geom_point(size = 0.7, alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, colour = col_up, fill = col_up, alpha = 0.15, linewidth = 0.6) +
  labs(tag = "", x = expression("Serum NfL immunoassay ("*log[2]~"pg/mL)"), y = expression("Plasma NEFL aptamer ("*log[2]~"RFU)"),
       subtitle = sprintf("%d draws from %d patients: r = %.2f\nWithin-patient change (n = %d): r = %.2f",
                          nfl_conc$n_matched, nfl_conc$n_subjects, nfl_conc$pearson_r, nfl_conc$within_n, nfl_conc$within_r)) +
  base + theme(
               axis.title.x = element_text(size = 12),
               axis.title.y = element_text(size = 12))
  

bb <- nfl_env$bb_ia %>% mutate(third = factor(c("Low", "Intermediate", "High")[tertile], levels = c("Low", "Intermediate", "High")))
bb_med <- bb %>% group_by(third) %>% summarise(nfl = 2^median(base_value), m = mean(ALSFRS_slope), n = n(), .groups = "drop")
kw <- kruskal.test(ALSFRS_slope ~ third, data = bb)
cat(sprintf("  Baseline NfL thirds: median %s pg/mL | mean slope %s | KW p = %.2g\n",
            paste(round(bb_med$nfl), collapse = ", "), paste(round(bb_med$m, 1), collapse = ", "), kw$p.value))
p7c <- ggplot(bb, aes(third, ALSFRS_slope, fill = third)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_violin(colour = NA, alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_text(data = bb_med, aes(third, max(bb$ALSFRS_slope) + 2, label = sprintf("%.0f pg/mL\nn = %d", nfl, n)),
            inherit.aes = FALSE, size = 3, lineheight = 0.9) +
  scale_fill_manual(values = c(Low = "#FDDBC7", Intermediate = "#EF8A62", High = "#B2182B"), guide = "none") +
  labs(tag = "", x = "Baseline serum NfL", y = "ALSFRS-R slope (points per year)",
       subtitle = sprintf("Kruskal-Wallis p = %s", fmt_p_plain(kw$p.value))) +
  base + theme(
               axis.text.x = element_text(size=9),
               axis.title.x = element_text(size = 12),
               axis.title.y = element_text(size = 12),
               legend.text = element_text(size = 9))
  

ia1 <- nfl_ia %>% filter(marker == "NEFL_immunoassay_allALS")
rd <- data.frame(measure = c("Baseline\nNfL", "Change\nin NfL", "Change in NfL,\nreliability-\ncorrected"),
                 r = c(ia1$base_r, ia1$slope_r, nfl_rel$disattenuated_r),
                 n = c(ia1$n_base, ia1$n_slope, ia1$n_slope)) %>%
  mutate(measure = factor(measure, levels = measure))

rd <- bind_cols(rd, fisher_ci(rd$r, rd$n))
k_rel <- sqrt(nfl_rel$reliability_nfl_slope * nfl_rel$reliability_alsfrs_slope)
rd[3, c("lo", "hi")] <- rd[2, c("lo", "hi")] / k_rel

p7d <- ggplot(rd, aes(measure, r)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_col(width = 0.6, fill = c("#B2182B", "grey55", "grey80")) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, na.rm = TRUE) +
  geom_text(aes(label = sprintf("r = %.2f\nn = %d", r, n), y = pmin(r, lo, na.rm = TRUE) - 0.07), size = 3, lineheight = 0.9) +
  scale_y_continuous(limits = c(-0.75, 0.1)) +
  labs(tag = "", x = NULL, y = "Correlation with ALSFRS-R slope",
       subtitle = sprintf("Slope reliability: NfL %.2f, ALSFRS-R %.2f", nfl_rel$reliability_nfl_slope, nfl_rel$reliability_alsfrs_slope)) +
  base + theme(
    axis.text.x = element_text(size=9),
    axis.title.y = element_text(size = 12))

pt7 <- nfl_env$pt %>% filter(!is.na(base_tertile)) %>%
  mutate(third = factor(c("Low", "Intermediate", "High")[base_tertile], levels = c("Low", "Intermediate", "High")))
tt7 <- nfl_env$tert_tbl %>% mutate(third = factor(c("Low", "Intermediate", "High")[base_tertile], levels = levels(pt7$third)),
                                   txt = sprintf("median %.0f pg/mL\nr = %.2f, p = %s", median_nfl_pg, r_slope_slope, fmt_p_plain(p)))
p7e <- ggplot(pt7, aes(ALSFRS_slope, nfl_slope)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(size = 0.7, alpha = 0.45) +
  geom_smooth(method = "lm", formula = y ~ x, colour = col_up, fill = col_up, alpha = 0.15, linewidth = 0.6) +
  facet_wrap(~ third, nrow = 1) +
  geom_text(data = tt7, aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.2, size = 3, lineheight = 0.9, inherit.aes = FALSE) +
  labs(tag = "", x = "ALSFRS-R slope (points per year)", y = expression("NfL slope ("*log[2]~"pg/mL per year)"),
       subtitle = "Within thirds of baseline serum NfL") +
  base + theme(axis.title.y = element_text(size = 12),
               axis.title.x = element_text(size = 12))

save_fig("Fig7.tiff", list(p7a, p7b, p7c, p7d, p7e),
         matrix(c(1, 1, 2, 2, 3, 3,
                  4, 4, 5, 5, 5, 5), nrow = 2, byrow = TRUE), width = 7.2, height = 6.2)


tiff(paste0(OUT, "Fig7A.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p7a)
dev.off()

tiff(paste0(OUT, "Fig7B.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p7b)
dev.off()

tiff(paste0(OUT, "Fig7C.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p7c)
dev.off()


tiff(paste0(OUT, "Fig7D.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p7d)
dev.off()


tiff(paste0(OUT, "Fig7E.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p7e)
dev.off()

# =============================================================================
# FIGURE 8. Composite against NfL, and loss of tracking in fast progressors
# =============================================================================
cat("-- Figure 8\n")

nr <- fs_env$nfl_r %>%
  mutate(nfl_third = factor(recode(nfl_third, Middle = "Intermediate"), levels = c("Low", "Intermediate", "High")),
         measure = recode(measure, `Plasma marker composite change` = "Change in progression markers (cross-fitted)",
                          `Serum NfL change` = "Change in serum NfL"))
nd <- fs_env$nfl_diff %>% mutate(nfl_third = factor(recode(nfl_third, Middle = "Intermediate"), levels = levels(nr$nfl_third)),
                                 txt = sprintf("diff. %.2f\np = %s", diff, fmt_p_plain(p)))
p8a <- ggplot(nr, aes(nfl_third, r, colour = measure)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.45), size = 0.5) +
  geom_text(data = nd, aes(nfl_third, 0.72, label = txt), inherit.aes = FALSE, size = 3, lineheight = 0.9) +
  scale_colour_manual(values = c("black", col_up), name = NULL) +
  scale_y_continuous(limits = c(-0.25, 0.8)) +
  labs(tag = "", x = "Baseline serum NfL", y = "Correlation with ALSFRS-R slope (95% CI)",
       subtitle = "Same 229 patients") +
  base + theme(legend.position = "bottom",
               axis.text.x = element_text(size=9),
               axis.title.x = element_text(size = 12),
               axis.title.y = element_text(size = 12),
               legend.text = element_text(size = 9))

cb <- fs_env$comp %>% mutate(tertile = factor(tertile, levels = c("Slow", "Middle", "Fast")))
ti <- fs_env$tert_ins %>% mutate(tertile = factor(tertile, levels = levels(cb$tertile)),
                                 txt = sprintf("beta = %.3f\n(%.3f to %.3f)", beta, lo, hi))
p8b <- ggplot(cb, aes(als_slope, composite)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, colour = col_down, fill = col_down, alpha = 0.15, linewidth = 0.6) +
  facet_wrap(~ tertile, nrow = 1, scales = "free_x",
             labeller = as_labeller(c(Slow = "Slowest third",Middle = "Middle third", Fast = "Fastest third" ))) +
  geom_text(data = ti, aes(x = -Inf, y = Inf, label = txt), hjust = -0.08, vjust = 1.3, size = 3, inherit.aes = FALSE) +
  labs(tag = "", x = "ALSFRS-R slope (points per year)", y = "Composite change (z)",
       subtitle = sprintf("Fastest against slowest interaction p = %s", fmt_p_plain(fs_env$p_ins))) +
  base + theme(
               axis.text.x = element_text(size=9),
               axis.title.x = element_text(size = 12),
               axis.title.y = element_text(size = 12))

pm <- fs_env$pm %>% filter(is.finite(ratio_fast), is.finite(ratio_slow))
n_lower <- sum(pm$ratio_fast < pm$ratio_slow)
p8c <- ggplot(pm, aes(ratio_fast, ratio_slow)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40") +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey70") +
  geom_point(alpha = 0.6, size = 1.3) +
  coord_cartesian(xlim = c(-1, 3), ylim = c(-1, 3)) +
  labs(tag = "", x = "Tracking in fastest third\n(relative to all patients)", y = "Tracking in slowest third\n(relative to all patients)",
       subtitle = sprintf("%d of %d markers above the line\nof equality (axes clipped)", n_lower, nrow(pm))) +
  base + theme(
    axis.text.x = element_text(size=9),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12))

tr <- fc_env$track %>% filter(tertile %in% c("Fast", "Slow")) %>%
  mutate(tertile = factor(ifelse(tertile == "Fast", "Fastest third", "Slowest third"),
                          levels = c("Slowest third", "Fastest third")),
         outcome = recode(outcome, `Respiratory (ALSFRS-R)` = "Respiratory"),
         outcome = factor(outcome, levels = rev(c("ALSFRS-R total", "Bulbar", "Fine motor", "Gross motor",
                                                  "Respiratory", "Vital capacity", "BMI"))))
p8d <- ggplot(tr, aes(r, outcome, colour = tertile)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.55), size = 0.4) +
  scale_colour_manual(values = c(`Fastest third` = col_up, `Slowest third` = col_down), name = NULL) +
  labs(tag = "", x = "Pearson r of composite change with decline (95% CI)", y = NULL,
       subtitle = "Cross-fitted composite") +
  base + theme(legend.position = "bottom",
               axis.text.x = element_text(size=9),
               axis.title.x = element_text(size = 12),
               axis.text.y = element_text(size = 12),
               legend.text = element_text(size = 9))

save_fig("Fig8.tiff", list(p8a, p8b, p8c, p8d),
         matrix(c(1, 2, 2, 3, 4, 4), nrow = 2, byrow = TRUE), width = 7.2, height = 6.4)

tiff(paste0(OUT, "Fig8A.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p8a)
dev.off()

tiff(paste0(OUT, "Fig8B.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p8b)
dev.off()

tiff(paste0(OUT, "Fig8C.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p8c)
dev.off()

tiff(paste0(OUT, "Fig8D.tiff"), width = 6, height = 7, units = "in", res = 300, compression = "lzw")
print(p8d)
dev.off()


cat("\n== SCRIPT 09 COMPLETE ==\n")
