# Shared code of the notebooks in analysis/ that make the figures and tables of the paper, with the penalised
# regression spline (REML) as the trajectory model.
#
# Every notebook sets the repository root as its working directory and sources this file. It sources the gene colours
# (Set_colors.R) and the spline code (spline_functions.R: curves, derivatives, summary statistics, ages at the VAF
# thresholds, consensus clusters). Figures are written to output/ (FigureN/, Supplementary_Figures/) and to
# figures/paper/, tables to results/paper/ and output/Tables/.
# Participant ids are always the ids of the paper (resources/SHCS_change_ids.txt: SHCS_1, ...).

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
  library(gridExtra)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(mgcv)
})

source("R/Set_colors.R")
source("R/spline_functions.R")


# ---- Folders ---------------------------------------------------------------------------------------------------------

paper_dir <- "output"   # the figures and tables as in the paper
rev <- list(results = "results/paper", figures = "figures/paper",
            lod_resources = "resources/LoD", lod_results = "results/LoD",
            cross_resources = "resources/LoD/cross_sectional")
paper <- list(supp = file.path(paper_dir, "Supplementary_Figures"), tables = file.path(paper_dir, "Tables"),
              figure = function(n) file.path(paper_dir, paste0("Figure", n)))
for (d in c(rev$results, rev$figures, paper$supp, paper$tables)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# a main figure (or a panel of it): output/FigureN/ and figures/paper/
save_main <- function(p, figure, name, width, height, dpi = 300, svg = TRUE) {
  dir.create(paper$figure(figure), showWarnings = FALSE, recursive = TRUE)
  for (d in c(paper$figure(figure), rev$figures)) {
    ggsave(file.path(d, paste0(name, ".png")), p, width = width, height = height, dpi = dpi, bg = "white")
    if (svg) ggsave(file.path(d, paste0(name, ".svg")), p, width = width, height = height, bg = "white")
  }
  invisible(p)
}

# a supplementary figure: output/Supplementary_Figures/ and figures/paper/
save_supp <- function(p, name, width, height, dpi = 300, svg = TRUE) {
  for (d in c(paper$supp, rev$figures)) {
    ggsave(file.path(d, paste0(name, ".png")), p, width = width, height = height, dpi = dpi, bg = "white")
    if (svg) ggsave(file.path(d, paste0(name, ".svg")), p, width = width, height = height, bg = "white")
  }
  invisible(p)
}

# a figure of an additional analysis that is not a display item of the paper: figures/paper/additional/ only
save_extra <- function(p, name, width, height, dpi = 300) {
  d <- file.path(rev$figures, "additional")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(d, paste0(name, ".png")), p, width = width, height = height, dpi = dpi, bg = "white")
  invisible(p)
}

# the separate panels of a multi-panel figure: FigureN/panels/<name>.png and .svg (and the same under figures/paper)
save_panel <- function(p, figure, name, width, height, dpi = 300) {
  d1 <- file.path(if (is.numeric(figure)) paper$figure(figure) else paper$supp, "panels")
  d2 <- file.path(rev$figures, "panels")
  for (d in c(d1, d2)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(d, paste0(name, ".png")), p, width = width, height = height, dpi = dpi, bg = "white")
    ggsave(file.path(d, paste0(name, ".svg")), p, width = width, height = height, bg = "white")
  }
  invisible(p)
}

write_res <- function(x, name, paper_copy = FALSE) {
  write.table(x, file.path(rev$results, name), sep = "\t", quote = FALSE, row.names = FALSE)
  if (paper_copy) write.table(x, file.path(paper$tables, name), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(x)
}

read_tsv_ <- function(file) read.table(file, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "") %>%
  as_tibble() %>% mutate(across(any_of("ID"), as.character))


# ---- Formatting ------------------------------------------------------------------------------------------------------

fmt_p <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), ifelse(p < 0.1, sprintf("%.2g", p), sprintf("%.2f", p))))
# p value as printed in the paper: p < 0.001 below 0.001, otherwise two significant digits
p_label <- function(p) ifelse(is.na(p), "p = -", ifelse(p < 0.001, "p < 0.001", paste("p =", fmt_p(p))))
stars_of <- function(p) case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
iqr_txt <- function(x, digits = 1) sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), median(x, na.rm = TRUE),
                                           quantile(x, 0.25, na.rm = TRUE), quantile(x, 0.75, na.rm = TRUE))
n_pct <- function(x) sprintf("%d (%.0f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))

# "DNMT3A p.R882C"; splice variants without a protein change get the cDNA change; `gene` overrides the gene field of
# the trajectory id (a splice variant carries its transcript there)
short_variant <- function(traj_id, gene = NULL) {
  g <- sub("^[^:]+:([^:]+):.*", "\\1", traj_id)
  if (!is.null(gene)) g <- gene
  p <- sub(".*:p\\.", "p.", traj_id)
  p <- ifelse(grepl("^p\\.", p), p, sub(".*:(c\\.[^:]+)$", "\\1", traj_id))
  paste(g, p)
}

source_shapes <- c(original = 16, bam_readcount = 17, imputed = 1)
source_labels <- c(original = "Mutect2 call", bam_readcount = "variant reads in the pile-up", imputed = "not detected (LoD/2)")
cohort_colors <- c(SardiNIA = "grey44", SHCS = "red3", Sardinia = "grey44")
cluster_labels <- c("1" = "\"not-growing\"\ncluster", "2" = "\"growing\"\ncluster")
cluster_fill <- c("#9ad7f5", "#ffb266")   # lighter versions of cluster_colors (Set_colors.R)
# colours of up to four clusters
cluster_colors_k <- c("1" = "#65CCFF", "2" = "#FF7F00", "3" = "#6A3D9A", "4" = "#33A02C")
cluster_fill_k <- c("1" = "#9ad7f5", "2" = "#ffb266", "3" = "#b39ddb", "4" = "#a5d6a7")

# the final consensus clusters (03_consensus_clusters.Rmd): traj_id, cluster (1..k), label
load_final_clusters <- function() {
  f <- file.path(rev$results, "clusters_final.txt")
  if (!file.exists(f)) stop("run 03_consensus_clusters.Rmd first")
  read_tsv_(f) %>% mutate(cluster = as.integer(cluster))
}
load_final_sard_clusters <- function() read_tsv_(file.path(rev$results, "clusters_final_sardinia.txt")) %>% mutate(cluster = as.integer(cluster))
# the cluster labels as a named vector ("1" = "shrinking", ...) and as factor levels in cluster order
cluster_label_map <- function(cl) { d <- cl %>% distinct(cluster, label) %>% arrange(cluster); setNames(d$label, as.character(d$cluster)) }
cluster_factor <- function(cluster, cl) { m <- cluster_label_map(cl); factor(m[as.character(cluster)], levels = m) }
# the same with the number of trajectories of each cluster in the label ("growing (n = ...)" or, with sep = "\n", on a second line)
cluster_factor_n <- function(cluster, cl, sep = " ") {
  f <- cluster_factor(cluster, cl); n_all <- table(cluster_factor(cl$cluster, cl))
  lab <- paste0(levels(f), sep, if (sep == " ") "(n = " else "n = ", n_all[levels(f)], if (sep == " ") ")" else "")
  factor(lab[as.integer(f)], levels = lab)
}
tag_theme <- theme(plot.tag = element_text(size = 26, face = "plain"))


# ---- Data ------------------------------------------------------------------------------------------------------------

load_ids <- function() read_tsv_("resources/SHCS_change_ids.txt") %>% mutate(ID = as.character(ID))

# the 71 final trajectories (longitudinal assay, LoD) with the new ids; AF at an imputed timepoint is LoD / 2
load_calls <- function() {
  read_tsv_("resources/chip_calls_long_filtered.txt") %>%
    inner_join(load_ids(), by = "ID") %>%
    mutate(source = factor(source, levels = names(source_shapes)), label = short_variant(traj_id, Gene.refGene)) %>%
    arrange(ID, traj_id, Age)
}

# REML curves within the observed ages (VAF floored at LoD / 2 of the trajectory; VAF_spline before the floor)
load_pred <- function(dir = rev$lod_resources) read_tsv_(file.path(dir, "SHCS_predicted_trajectories_derivatives.txt")) %>%
  inner_join(load_ids(), by = "ID")
# curves extrapolated back to LoD / 2 of the trajectory
load_pred2 <- function(dir = rev$lod_resources) read_tsv_(file.path(dir, "SHCS_predicted_trajectories_derivatives2_LoD2_extrap.txt")) %>%
  inner_join(load_ids(), by = "ID")
load_sum <- function(dir = rev$lod_resources) read_tsv_(file.path(dir, "SHCS_predicted_trajectories_summary.txt")) %>%
  inner_join(load_ids(), by = "ID") %>%
  mutate(Gene.refGene = group_genes(Gene.refGene))
load_ages <- function(dir = rev$lod_resources) read_tsv_(file.path(dir, "SHCS_trajectories_threshold_ages.txt"))
# optional: an earlier two-cluster consensus (1 = "not-growing", 2 = "growing"), used only for the label order of k = 2
load_clusters <- function() read_tsv_("results/clustering/SHCS_DTW_consensus_k2.txt")

load_sard_calls <- function(dir = rev$lod_resources) read_tsv_(file.path(dir, "sardinia_chip_calls_long_filtered.txt")) %>%
  mutate(Chr = as.character(Chr), source = factor(source, levels = names(source_shapes)))
load_sard_pred <- function(dir = rev$lod_resources) read_tsv_(file.path(dir, "Sardinia_predicted_trajectories_derivatives.txt"))
load_sard_sum <- function(dir = rev$lod_resources) read_tsv_(file.path(dir, "Sardinia_predicted_trajectories_summary.txt")) %>%
  mutate(Gene.refGene = group_genes(Gene.refGene))
load_sard_ages <- function(dir = rev$lod_resources) read_tsv_(file.path(dir, "Sardinia_trajectories_threshold_ages.txt"))

load_fabre <- function() read.csv("resources/SardiNIA/Sardinia_growth_coefficients_and_age_at_onset_estimates_for_all_clones.csv", check.names = FALSE)

load_meta <- function() read_tsv_("resources/SHCS_Metadata_univ.txt")

# the four comparison sets of the SHCS vs SardiNIA analysis (curves at every half year, with Cohort), 5 main genes or all
set_keys <- c("unmatched", "unmatched_filtered", "matched", "matched_filtered")
set_files <- c(unmatched = "unmatched", unmatched_filtered = "unmatched_age_filtered", matched = "matched", matched_filtered = "matched_age_filtered")
set_labels <- c(unmatched = "Unmatched", unmatched_filtered = "Unmatched, age-filtered", matched = "Matched", matched_filtered = "Matched, age-filtered")
read_set <- function(key, genes = c("5", "all"), dir = rev$lod_resources) {
  genes <- match.arg(genes)
  read_tsv_(file.path(dir, sprintf("SHCS_Sardinia_%s_trajectories%s.txt", set_files[[key]], if (genes == "all") "_all_genes" else ""))) %>%
    mutate(Cohort = ifelse(Cohort == "Sardinia", "SardiNIA", "SHCS"))
}

# the interpolated blood counts of the SHCS database on the half-year grid of each individual, with the age at ART start
blood_cols <- c("LEU", "HEM", "PLA", "LYM", "CD3", "CD3P", "CD4", "CD4P", "CD8", "CD8P", "CD4_CD8_ratio")
blood_labels <- c(LEU = "LEU", HEM = "HEM", PLA = "PLT", LYM = "LYM", CD3 = "CD3", CD3P = "CD3P", CD4 = "CD4", CD4P = "CD4P",
                  CD8 = "CD8", CD8P = "CD8P", CD4_CD8_ratio = "CD4:CD8")
load_blood_grid <- function() read_tsv_("resources/SHCS_CHIP_blood_trajectories_merged.txt") %>%
  distinct(ID, Age, AgeArtStart, across(all_of(blood_cols)))


# ---- Uncertainty of the per-trajectory statistics -------------------------------------------------------------------
#
# Read-count bootstrap: at every measured timepoint
# the variant reads are redrawn, alt ~ Binomial(depth, VAF), the imputed timepoints keep LoD / 2, and the trajectory is
# refitted with the main spline. Per trajectory: 95% percentile intervals of the statistics and the share of replicates
# with the sign of the main analysis (a sign is "robust" when it is at least 95%).

bootstrap_stats <- function(calls, B = 200, seed = 2026, spec = spec_main) {
  set.seed(seed)
  calls %>%
    group_split(traj_id) %>%
    map_dfr(function(obs) {
      obs <- arrange(obs, Age)
      meas <- obs$source != "imputed"
      reps <- t(replicate(B, {
        af <- obs$AF
        af[meas] <- rbinom(sum(meas), obs$lod_depth[meas], obs$AF[meas]) / obs$lod_depth[meas]
        traj_stats(obs$Age, af, spec, obs$traj_lod_mean[1])
      }))
      as_tibble(reps) %>% mutate(traj_id = obs$traj_id[1], b = row_number(), .before = 1)
    })
}

bootstrap_ci <- function(boot, sum_main, stats = c("mean_VAF_speed", "mean_VAF_acceleration", "max_VAF_speed", "max_VAF", "Age_0.02", "Age_LoD2")) {
  boot %>%
    select(traj_id, b, all_of(stats)) %>%
    pivot_longer(all_of(stats), names_to = "stat", values_to = "value") %>%
    inner_join(sum_main %>% select(traj_id, all_of(stats)) %>% pivot_longer(all_of(stats), names_to = "stat", values_to = "main"),
               by = c("traj_id", "stat")) %>%
    group_by(traj_id, stat) %>%
    summarise(main = first(main), lo = quantile(value, 0.025, na.rm = TRUE), hi = quantile(value, 0.975, na.rm = TRUE),
              same_sign = mean(sign(value) == sign(first(main)), na.rm = TRUE), estimated = mean(!is.na(value)), .groups = "drop") %>%
    mutate(robust = !is.na(main) & same_sign >= 0.95)
}


# ---- Figures in the style of the paper --------------------------------------------------------------------------------

# Figure 2a-e: one statistic by gene (and all genes), boxplots in the gene colours
box_by_gene <- function(sum_df, stat, ylab, hline = 0, size = 17) {
  d <- bind_rows(sum_df %>% mutate(Gene.refGene = as.character(group_genes(Gene.refGene))),
                 sum_df %>% mutate(Gene.refGene = "All genes")) %>%
    mutate(Gene.refGene = factor(Gene.refGene, levels = c("All genes", gene_levels))) %>%
    filter(!is.na(.data[[stat]]))
  p <- ggplot(d, aes(y = .data[[stat]], x = Gene.refGene)) +
    geom_boxplot(aes(fill = Gene.refGene), outlier.size = 0.8) +
    scale_fill_manual(values = gene_colors) +
    theme_classic() + xlab("") + ylab(ylab) +
    theme(axis.text.x = element_text(size = size - 3, angle = 90, hjust = 1, vjust = 0.5), legend.position = "none",
          axis.title.y = element_text(size = size - 1), axis.text.y = element_text(size = size))
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed")
  p
}

# Spearman correlation heatmap of the summary statistics (Figure 2f), lower triangle
corr_heatmap <- function(df, col_order, labels = NULL, size = 16) {
  rc <- Hmisc::rcorr(as.matrix(df[, col_order]), type = "spearman")
  if (is.null(labels)) labels <- setNames(col_order, col_order)
  plot_df <- as.data.frame(rc$r) %>%
    tibble::rownames_to_column("Var1") %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "r") %>%
    mutate(p = as.vector(rc$P), stars = stars_of(p),
           i = match(Var1, col_order), j = match(Var2, col_order),
           Var1 = factor(Var1, levels = rev(col_order), labels = rev(labels[col_order])),
           Var2 = factor(Var2, levels = col_order, labels = labels[col_order]),
           r = ifelse(i < j, NA_real_, r), stars = ifelse(i < j, "", stars))
  ggplot() +
    geom_tile(data = subset(plot_df, !is.na(r)), aes(x = Var2, y = Var1, fill = r), color = "white", linewidth = 0.3) +
    geom_text(data = subset(plot_df, !is.na(r)), aes(x = Var2, y = Var1, label = stars), size = 5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), oob = squish, name = "Correlation") +
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(), text = element_text(size = size),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = size - 2), plot.title = element_text(face = "bold"))
}

stat_labels <- c(N_CHIP_more1 = "N CHIP > 1", mean_VAF = "mean VAF", sd_VAF = "sd VAF", min_VAF = "min VAF", max_VAF = "max VAF",
                 mean_VAF_speed = "mean VAF speed", sd_VAF_speed = "sd VAF speed", min_VAF_speed = "min VAF speed", max_VAF_speed = "max VAF speed",
                 mean_VAF_acceleration = "mean VAF acceleration", sd_VAF_acceleration = "sd VAF acceleration",
                 min_VAF_acceleration = "min VAF acceleration", max_VAF_acceleration = "max VAF acceleration",
                 mean_Age = "mean Age", min_Age = "min Age", max_Age = "max Age", Time = "Follow-up time",
                 Age_zero = "Age at 0.1% VAF", Age_LoD2 = "Age at LoD/2", Age_LoD = "Age at LoD", Age_0.02 = "Age at 2% VAF")

# Figure 2h: the curves by cluster, coloured by gene
cluster_curves_plot <- function(pred, clusters, size = 15, nrow = 1) {
  pred %>%
    inner_join(clusters %>% select(traj_id, cluster), by = "traj_id") %>%
    mutate(cluster = cluster_factor_n(cluster, clusters), Gene = group_genes(Gene.refGene)) %>%
    ggplot(aes(x = Age, y = VAF, group = traj_id, color = Gene)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~cluster, nrow = nrow) +
    scale_color_manual(values = gene_colors) +
    theme_classic(base_size = size) + theme(strip.background = element_blank(), legend.position = "none")
}

# Figure 2j: gene composition of the clusters
cluster_genes_plot <- function(clusters, gene_tbl, size = 15) {
  clusters %>%
    inner_join(gene_tbl, by = "traj_id") %>%
    mutate(Gene = factor(as.character(group_genes(Gene)), levels = rev(gene_levels)),
           cluster = cluster_factor_n(cluster, clusters, sep = "\n")) %>%
    count(cluster, Gene) %>%
    group_by(cluster) %>% mutate(Prop = n / sum(n)) %>% ungroup() %>%
    ggplot(aes(x = cluster, y = Prop, fill = Gene)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = gene_colors) +
    labs(x = NULL, y = "Proportion") +
    theme_bw(base_size = size) + theme(legend.position = "none")
}

# Figure 2i and Figure S7: one statistic by cluster with the p value of the test;
# the test is Wilcoxon for two clusters and Kruskal-Wallis for more; `clusters` gives the labels
cluster_box <- function(sum_df, stat, ylab, hline = 0, size = 15, clusters = NULL) {
  d <- sum_df %>% filter(!is.na(.data[[stat]]))
  d$cluster <- if (is.null(clusters)) factor(cluster_labels[as.character(d$cluster)], levels = cluster_labels) else cluster_factor(d$cluster, clusters)
  k <- nlevels(d$cluster)
  p <- if (k == 2) wilcox.test(d[[stat]] ~ d$cluster, exact = FALSE)$p.value else kruskal.test(d[[stat]] ~ d$cluster)$p.value
  fills <- setNames(cluster_fill_k[seq_len(k)], levels(d$cluster))
  n_lev <- table(d$cluster); lab_n <- paste0(levels(d$cluster), "\nn = ", n_lev[levels(d$cluster)])
  names(fills) <- lab_n; d$cluster <- factor(lab_n[as.integer(d$cluster)], levels = lab_n)
  pl <- ggplot(d, aes(x = cluster, y = .data[[stat]])) +
    geom_boxplot(aes(fill = cluster), outlier.size = 0.8) +
    scale_fill_manual(values = fills) +
    annotate("text", x = (k + 1) / 2, y = Inf, vjust = 1.5, label = p_label(p), size = 4.2) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = size) + theme(legend.position = "none", axis.text.x = element_text(size = size - 3))
  if (!is.null(hline)) pl <- pl + geom_hline(yintercept = hline, linetype = "dashed")
  pl
}

# Figure 4 and Figures S15-S18: one statistic by cohort with the Wilcoxon p
cohort_box <- function(s, stat, ylab, hline = NULL, size = 14) {
  d <- s %>% filter(!is.na(.data[[stat]])) %>% mutate(Cohort = factor(as.character(Cohort), levels = c("SardiNIA", "SHCS")))
  p <- wilcox.test(d[[stat]] ~ d$Cohort, exact = FALSE)$p.value
  pl <- ggplot(d, aes(y = .data[[stat]], x = Cohort, fill = Cohort)) +
    geom_boxplot(outlier.size = 0.6) +
    annotate("text", x = 1.5, y = Inf, vjust = 1.5, label = p_label(p), size = 3.8) +
    scale_fill_manual(values = cohort_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  if (!is.null(hline)) pl <- pl + geom_hline(yintercept = hline, linetype = "dashed")
  pl
}

# summary statistics of the curves of a comparison set (traj_summary() on the curves of the set, the ages at the
# thresholds from the full extrapolated curves)
set_summary <- function(d, ages) {
  traj_summary(d %>% rename(Gene.refGene = Gene), ages) %>%
    rename(Gene = Gene.refGene) %>%
    left_join(d %>% distinct(traj_id, Cohort, traj_lod), by = "traj_id")
}


# ---- Mixed models of log10 VAF on age and cohort (Figure 4) -----------------------------------------------------------

age_centre <- 65

# the curve at the measured ages (data_point == 1); group_filter: the rows of the group; zero: no age x cohort interaction
fit_cohort_model <- function(data, group_filter, random_slope = TRUE, zero = TRUE) {
  df <- data %>% filter({{ group_filter }}) %>%
    mutate(Age_c = Age - age_centre, Cohort = factor(as.character(Cohort), levels = c("SardiNIA", "SHCS")))
  fixed <- if (zero) "log10(VAF) ~ Age_c + Cohort" else "log10(VAF) ~ Age_c * Cohort"
  random <- if (random_slope) "(Age_c | traj_id)" else "(1 | traj_id)"
  model <- tryCatch(lmer(as.formula(paste(fixed, "+", random)), data = df, control = lmerControl(optimizer = "bobyqa")), error = function(e) e)
  if (inherits(model, "error")) return(list(error = conditionMessage(model), data = df))
  list(model = model, data = df)
}

cohort_coef_table <- function(m) {
  tidy(m$model, effects = "fixed", conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    transmute(Term = factor(gsub("CohortSHCS", "HIV", gsub("Age_c", "Age", term)), levels = c("Age", "HIV", "Age:HIV")),
              Estimate = estimate, `CI Lower` = conf.low, `CI Upper` = conf.high, `p-value` = p.value, Sig = stars_of(p.value)) %>%
    arrange(Term) %>% mutate(Term = as.character(Term))
}

all_cohort_coefficients <- function(results) {
  imap_dfr(results, function(m, label) {
    if (!is.null(m$error)) return(NULL)
    cohort_coef_table(m) %>%
      mutate(label = label, trajectories = n_distinct(m$data$traj_id), measurements = nrow(m$data), singular = isSingular(m$model),
             converged = length(m$model@optinfo$conv$lme4$messages) == 0, REML = REMLcrit(m$model), .before = 1)
  })
}

# Figure 4b-e and Figures S12-S14: the measured VAF of every trajectory (points and thin lines) on a log
# scale and the population lines of the model, with the estimates printed in the panel
cohort_model_panel <- function(m, title = NULL, ylim = c(1e-3, 0.7), size = 14, show_points = TRUE, text_size = 3.4) {
  fe <- fixef(m$model)
  ln <- tibble(Age = seq(40, 90, by = 1)) %>%
    crossing(Cohort = c("SardiNIA", "SHCS")) %>%
    mutate(VAF = 10^(fe["(Intercept)"] + fe["Age_c"] * (Age - age_centre) + fe["CohortSHCS"] * (Cohort == "SHCS") +
                       ifelse("Age_c:CohortSHCS" %in% names(fe), fe["Age_c:CohortSHCS"], 0) * (Age - age_centre) * (Cohort == "SHCS")))
  est <- cohort_coef_table(m) %>%
    transmute(txt = sprintf("%s: %.3f (%.3f to %.3f), %s", Term, Estimate, `CI Lower`, `CI Upper`, p_label(`p-value`)))
  n_coh <- m$data %>% distinct(traj_id, Cohort) %>% count(Cohort)
  est <- bind_rows(est, tibble(txt = sprintf("SHCS N=%d, SardiNIA N=%d", sum(n_coh$n[n_coh$Cohort == "SHCS"]), sum(n_coh$n[n_coh$Cohort == "SardiNIA"]))))
  p <- ggplot()
  if (show_points) p <- p +
    geom_line(data = m$data, aes(x = Age, y = VAF, group = traj_id, color = Cohort), alpha = 0.3, linewidth = 0.3) +
    geom_point(data = m$data, aes(x = Age, y = VAF, color = Cohort), size = 0.7, alpha = 0.6)
  p +
    geom_line(data = ln, aes(x = Age, y = VAF, color = Cohort), linewidth = 1.5) +
    annotate("label", x = 40, y = ylim[2], hjust = 0, vjust = 1, size = text_size, label = paste(est$txt, collapse = "\n"),
             label.size = 0, fill = alpha("white", 0.7)) +
    scale_y_log10(labels = function(x) paste0(100 * x, "%"), limits = ylim) +
    scale_color_manual(values = cohort_colors, name = "Cohort") +
    labs(title = title, y = "VAF", x = "Age") +
    theme_classic(base_size = size) + theme(legend.position = "bottom", plot.title = element_text(size = size))
}


# ---- Blood parameters (Figure 3) --------------------------------------------------------------------------------------

outcome_labels <- c(VAF = "VAF", VAF_speed = "VAF speed", VAF_acceleration = "VAF acceleration")

# the curves of one trajectory per person (the largest clone by mean VAF) joined to the blood counts, after ART start
blood_data <- function(pred, blood) {
  pred %>%
    select(ID, SHCS_ID_new, traj_id, Age, VAF, VAF_speed, VAF_acceleration, data_point, Gene.refGene) %>%
    inner_join(blood, by = c("ID", "Age")) %>%
    mutate(StartedART = (Age - AgeArtStart) > 0.5) %>%
    filter(StartedART) %>%
    group_by(traj_id) %>% mutate(mv = mean(VAF)) %>% group_by(ID) %>% filter(mv == max(mv)) %>% ungroup() %>%
    select(-mv)
}

# outcome ~ age (z) x log10(parameter) [+ log10 leukocytes + log10 lymphocytes]; on the half-year grid: random intercept and
# slope per trajectory; on the measured timepoints: random intercept per person. ML fit; the coefficient of the parameter
fit_blood_model <- function(d, outcome, blood, adjust = FALSE, measured = FALSE) {
  d <- d %>% mutate(y = if (outcome == "VAF") log10(VAF + 1e-4) else .data[[outcome]],
                    Age_z = as.numeric(scale(Age)), x = log10(.data[[blood]] + 1e-4),
                    LEU_l = log10(LEU + 1e-4), LYM_l = log10(LYM + 1e-4)) %>%
    filter(!is.na(x), !is.na(y))
  adj <- if (!adjust) "" else if (blood == "LEU") "+ LYM_l" else if (blood == "LYM") "+ LEU_l" else "+ LEU_l + LYM_l"
  if (adjust) d <- d %>% filter(!is.na(LEU_l), !is.na(LYM_l))
  re <- if (measured) "(1 | ID)" else "(Age_z | traj_id)"
  m <- tryCatch(lmer(as.formula(paste("y ~ Age_z * x", adj, "+", re)), data = d, REML = FALSE, control = lmerControl(optimizer = "bobyqa")),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  tidy(m, effects = "fixed", conf.int = TRUE) %>% filter(term == "x") %>%
    transmute(outcome, parameter = blood, estimate, conf.low, conf.high, p.value, n = nrow(d), n_persons = n_distinct(d$ID))
}

blood_models <- function(d, adjust = FALSE, measured = FALSE, outcomes = names(outcome_labels)) {
  map_dfr(outcomes, function(o) map_dfr(blood_cols, function(b) fit_blood_model(d, o, b, adjust, measured))) %>%
    mutate(sig = stars_of(p.value))
}

# Figure 3 and Figures S9 and S10: one column per outcome, squares with 95% CI, stars
blood_forest <- function(res, color = "deepskyblue2", size = 15, outcomes = names(outcome_labels)) {
  d <- res %>%
    filter(outcome %in% outcomes) %>%
    mutate(parameter = factor(parameter, levels = rev(blood_cols), labels = rev(blood_labels[blood_cols])),
           outcome = factor(outcome, levels = names(outcome_labels), labels = outcome_labels))
  ggplot(d, aes(x = estimate, y = parameter, xmin = conf.low, xmax = conf.high)) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.6) +
    geom_errorbarh(height = 0.2, color = color) +
    geom_point(size = 2.6, shape = 15, color = color) +
    geom_text(aes(label = sig), vjust = -0.5, size = 4) +
    facet_wrap(~outcome, scales = "free_x") +
    scale_y_discrete(expand = expansion(add = c(0.6, 1))) +   # room for the significance mark of the top row
    scale_x_continuous(breaks = scales::breaks_extended(4)) +
    labs(x = "Estimate (95% CI)", y = NULL) +
    theme_classic(base_size = size) + theme(strip.background = element_rect(linewidth = 0.6), panel.spacing.x = unit(1.6, "lines"))
}
