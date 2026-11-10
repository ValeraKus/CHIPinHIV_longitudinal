# Spline functions shared by the notebooks in analysis/
#
# A trajectory is a data frame of measurements with columns ID, traj_id, N_CHIP, Gene.refGene, Age, AF (and LABYEAR,
# optional). Each trajectory is smoothed on its own; the curve is predicted every `step` years and the derivatives are
# finite differences of the curve over `step` years.

library(tidyverse)

step     <- 0.5     # years between two predicted points, and over which the derivatives are taken
zero_vaf <- 0.001   # the extrapolated curve starts where the spline is below this VAF (age of origin)
chip_vaf <- 0.02    # CHIP VAF threshold
eps_group <- 1e-4   # |mean speed| or |mean acceleration| below this is "stable" / "constant"

gene_levels <- c("DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "Other")
other_genes <- c("GNAS", "JAK2", "SF3B1", "SRSF2", "KRAS")
main_genes  <- c("DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53")   # the genes of the SHCS vs SardiNIA and Fabre comparisons
group_genes <- function(x) factor(ifelse(x %in% other_genes, "Other", as.character(x)), levels = gene_levels)


# ---- Spline specifications ------------------------------------------------------------------------------------------

# spar: smooth.spline() smoothing parameter; df: smooth.spline() with fixed effective degrees of freedom (at most the
# number of distinct ages); reml: penalised regression spline mgcv::gam(AF ~ s(Age, k = 10), method = "REML"), with k at
# most the number of distinct ages
spec_grid <- list(
  "spar 0.1" = list(type = "spar", value = 0.1),
  "spar 0.2" = list(type = "spar", value = 0.2),
  "spar 0.3" = list(type = "spar", value = 0.3),
  "spar 0.4" = list(type = "spar", value = 0.4),
  "spar 0.5" = list(type = "spar", value = 0.5),
  "df 2"     = list(type = "df", value = 2),
  "df 2.5"   = list(type = "df", value = 2.5),
  "df 3"     = list(type = "df", value = 3),
  "df 4"     = list(type = "df", value = 4),
  "REML"     = list(type = "reml", value = 10))
spec_original <- "spar 0.2"   # the smooth.spline() setting of an earlier version of the analysis

# Main specification (compared with the others in 4_sensitivity_analyses/01_spline_sensitivity_and_simulation.Rmd):
# the hold-out error of the specifications of the grid is the same within noise (Table S4), so the smoothness is not
# tuned globally but estimated for each trajectory by REML
spec_main_name <- "REML"
spec_main <- spec_grid[[spec_main_name]]


# ---- Fitting ----------------------------------------------------------------------------------------------------------

# a trajectory with 3 measurements gets the midpoint of its first two measurements as a 4th point (the spline needs 4
# distinct ages); fewer than 3 measurements cannot be fitted
fit_points <- function(age, af) {
  o <- order(age)
  age <- age[o]
  af <- af[o]
  if (length(age) < 3) return(NULL)
  if (length(age) < 4) {
    age <- c(age, (age[1] + age[2]) / 2)
    af <- c(af, (af[1] + af[2]) / 2)
  }
  list(x = age, y = af)
}

# returns the fitted curve as a function of age and the effective degrees of freedom
fit_traj <- function(age, af, spec) {
  p <- fit_points(age, af)
  if (is.null(p)) return(NULL)
  n_x <- length(unique(p$x))
  if (spec$type == "spar") {
    f <- smooth.spline(x = p$x, y = p$y, spar = spec$value)
    list(predict = function(a) predict(f, x = a)$y, df = f$df, model = f)
  } else if (spec$type == "df") {
    f <- smooth.spline(x = p$x, y = p$y, df = min(spec$value, n_x))
    list(predict = function(a) predict(f, x = a)$y, df = f$df, model = f)
  } else if (spec$type == "reml") {
    d <- data.frame(Age = p$x, AF = p$y)
    m <- mgcv::gam(AF ~ s(Age, k = min(spec$value, n_x)), data = d, method = "REML")
    list(predict = function(a) as.numeric(predict(m, newdata = data.frame(Age = a))), df = sum(m$edf), model = m)
  } else stop("unknown spec type ", spec$type)
}

# VAF (negative values set to 0), speed (forward difference over `step`) and acceleration (difference of consecutive
# speeds / `step`; 0 at the first age by construction) at the ages `ages`
spline_curves <- function(pfun, ages) {
  n <- length(ages)
  y <- pfun(c(ages, ages + step))
  speed <- (y[(n + 1):(2 * n)] - y[1:n]) / step
  tibble(Age = ages, VAF = pmax(y[1:n], 0), VAF_speed = speed, VAF_acceleration = c(0, diff(speed)) / step)
}

# ages of the predicted curve: every `step` years from the first to the last measurement; with extrapolate = TRUE and a
# VAF > zero_vaf at the first measurement, the curve starts at the last age before it at which the spline is < zero_vaf
# (searched down to age 0)
curve_ages <- function(age, pfun, extrapolate = FALSE) {
  ages <- seq(min(age), max(age), by = step)
  if (extrapolate && max(pfun(min(ages)), 0) > zero_vaf) {
    ages0 <- seq(0, max(age), by = step)
    below <- which(pfun(ages0) < zero_vaf)
    below <- below[below < which(ages0 == min(age))]
    if (length(below) > 0) ages <- ages0[max(below):length(ages0)]
  }
  ages
}

predict_trajectory <- function(obs, spec, extrapolate = FALSE) {
  obs <- arrange(obs, Age)
  f <- fit_traj(obs$Age, obs$AF, spec)
  if (is.null(f)) return(NULL)
  out <- spline_curves(f$predict, curve_ages(obs$Age, f$predict, extrapolate)) %>%
    mutate(ID = obs$ID[1], traj_id = obs$traj_id[1], N_CHIP = obs$N_CHIP[1], Gene.refGene = obs$Gene.refGene[1],
           data_point = as.numeric(Age %in% obs$Age),
           out_of_range = as.numeric(!(Age %in% obs$Age)),   # 1 for every age without a measurement
           df = f$df)
  if ("LABYEAR" %in% names(obs)) out$LABDATE <- out$Age + obs$LABYEAR[1] - obs$Age[1]
  out %>%
    arrange(Age) %>%
    mutate(Time = 0:(n() - 1), VAF_acceleration_sign = sign(VAF_acceleration))
}

predict_all <- function(calls, spec, extrapolate = FALSE) {
  calls %>% group_split(traj_id) %>% map(~ predict_trajectory(.x, spec, extrapolate)) %>% bind_rows()
}


# ---- Ages at VAF thresholds --------------------------------------------------------------------------------------------

# Ages on the rise of the extrapolated curve to CHIP, i.e. up to the first age at which it reaches 2% (the whole curve
# if it never does):
#   Age_zero = last age at which the curve is <= zero_vaf (origin of the clone)
#   Age_LoD2, Age_LoD = last age at which the curve crosses LoD / 2 and LoD upwards (linear interpolation between grid points)
#   Age_0.02 = the age at which it reaches 2%
# NA (left-censored) when the curve is above the threshold at every age of the rise, i.e. it starts above it. For a
# curve that rises steadily this is the first age at which the curve reaches the threshold; for a curve that starts just
# above a threshold (e.g. at an imputed LoD / 2 point), dips below it and rises again, it is the start of the final rise,
# and the four ages are always ordered.
threshold_ages_one <- function(age, vaf, lod) {
  o <- order(age)
  age <- age[o]
  vaf <- vaf[o]
  i2 <- which(vaf >= chip_vaf)[1]
  reaches <- !is.na(i2)
  if (!reaches) i2 <- length(vaf)
  a <- age[1:i2]
  v <- vaf[1:i2]
  last_cross <- function(thr) {
    if (is.na(thr) || length(v) < 2) return(NA_real_)
    i <- which(v[-1] >= thr & v[-length(v)] < thr) + 1
    if (length(i) == 0) return(NA_real_)
    i <- max(i)
    a[i - 1] + (thr - v[i - 1]) / (v[i] - v[i - 1]) * (a[i] - a[i - 1])
  }
  z <- which(v <= zero_vaf)
  c(Age_zero = if (length(z) > 0) a[max(z)] else NA_real_,
    Age_LoD2 = last_cross(lod / 2),
    Age_LoD  = last_cross(lod),
    Age_0.02 = if (reaches) last_cross(chip_vaf) else NA_real_)
}

# pred2: extrapolated curves; lod: traj_id and lod (mean LoD of the trajectory; NA when unknown)
threshold_ages <- function(pred2, lod) {
  pred2 %>%
    left_join(lod, by = "traj_id") %>%
    group_by(traj_id) %>%
    summarise(as_tibble(as.list(threshold_ages_one(Age, VAF, first(lod)))), .groups = "drop")
}


# ---- Summary statistics ------------------------------------------------------------------------------------------------

dynamics_group <- function(mean_speed, mean_acc) {
  case_when(mean_speed < -eps_group               ~ "Shrinking",
            abs(mean_speed) <= eps_group          ~ "Stable",
            abs(mean_acc) <= eps_group            ~ "Growing_constant",
            mean_acc > eps_group                  ~ "Growing_accel",
            TRUE                                  ~ "Growing_decel")
}

# one row per trajectory: mean / sd / min / max of VAF, speed, acceleration and of the sign of the acceleration over the
# predicted ages, ages, follow-up span, the ages at VAF thresholds and the dynamics group
traj_summary <- function(pred, ages) {
  pred %>%
    group_by(traj_id, ID, Gene.refGene) %>%
    summarise(N_CHIP = max(N_CHIP),
              mean_VAF_acceleration = mean(VAF_acceleration), sd_VAF_acceleration = sd(VAF_acceleration),
              max_VAF_acceleration = max(VAF_acceleration), min_VAF_acceleration = min(VAF_acceleration),
              mean_VAF_speed = mean(VAF_speed), sd_VAF_speed = sd(VAF_speed),
              max_VAF_speed = max(VAF_speed), min_VAF_speed = min(VAF_speed),
              mean_VAF_acceleration_sign = mean(VAF_acceleration_sign), sd_VAF_acceleration_sign = sd(VAF_acceleration_sign),
              max_VAF_acceleration_sign = max(VAF_acceleration_sign), min_VAF_acceleration_sign = min(VAF_acceleration_sign),
              mean_VAF = mean(VAF), sd_VAF = sd(VAF), min_VAF = min(VAF), max_VAF = max(VAF),
              VAF_more5perc = as.numeric(any(VAF >= 0.05)),
              min_Age = min(Age), max_Age = max(Age), mean_Age = mean(Age), .groups = "drop") %>%
    mutate(deltaVAF = max_VAF - min_VAF, Time = max_Age - min_Age, N_CHIP_more1 = as.numeric(N_CHIP > 1),
           Gene.refGene = group_genes(Gene.refGene)) %>%
    left_join(ages, by = "traj_id") %>%
    mutate(group = dynamics_group(mean_VAF_speed, mean_VAF_acceleration))
}

# everything derived from one set of trajectories and one specification; `trim` (optional) keeps only part of the
# curves within the observed ages before the summary, the ages at thresholds are always from the full
# extrapolated curves
run_pipeline <- function(calls, spec, lod = NULL, trim = NULL) {
  p1 <- predict_all(calls, spec, FALSE)
  p2 <- predict_all(calls, spec, TRUE)
  if (is.null(lod)) lod <- tibble(traj_id = unique(calls$traj_id), lod = NA_real_)
  ages <- threshold_ages(p2, lod)
  if (!is.null(trim)) p1 <- trim(p1, calls)
  list(pred = p1, pred2 = p2, ages = ages, sum = traj_summary(p1, ages))
}

# summary statistics of one trajectory without the tables (for the bootstrap); same values as
# traj_summary() + threshold_ages()
traj_stats <- function(age, af, spec, lod = NA_real_) {
  f <- fit_traj(age, af, spec)
  if (is.null(f)) return(NULL)
  c1 <- spline_curves(f$predict, curve_ages(age, f$predict, FALSE))
  a2 <- curve_ages(age, f$predict, TRUE)
  c(mean_VAF_speed = mean(c1$VAF_speed), max_VAF_speed = max(c1$VAF_speed), min_VAF_speed = min(c1$VAF_speed),
    sd_VAF_speed = sd(c1$VAF_speed), mean_VAF_acceleration = mean(c1$VAF_acceleration),
    max_VAF_acceleration = max(c1$VAF_acceleration), min_VAF_acceleration = min(c1$VAF_acceleration),
    max_VAF = max(c1$VAF), threshold_ages_one(a2, pmax(f$predict(a2), 0), lod), df = f$df)
}


# ---- Choice of the specification -------------------------------------------------------------------------------------

# hold out each internal measurement (not the first or the last of the trajectory), fit the others, predict it; only
# trajectories with >= 3 measurements left. err = |prediction - VAF| / range of the measured VAF of the trajectory.
# "linear" = linear interpolation between the two neighbours, as a benchmark
holdout_errors <- function(calls, specs = spec_grid) {
  calls %>%
    group_split(traj_id) %>%
    map_dfr(function(obs) {
      obs <- arrange(obs, Age)
      n <- nrow(obs)
      if (n < 4) return(NULL)
      rng <- diff(range(obs$AF))
      map_dfr(2:(n - 1), function(i) {
        rest <- obs[-i, ]
        pred <- c(map_dbl(specs, function(s) fit_traj(rest$Age, rest$AF, s)$predict(obs$Age[i])),
                  linear = approx(rest$Age, rest$AF, xout = obs$Age[i])$y)
        tibble(traj_id = obs$traj_id[1], held_out = i, n = n,
               position = case_when(i == 2 ~ "second", i == n - 1 ~ "second-to-last", TRUE ~ "middle"),
               source = obs$source[i], AF = obs$AF[i], depth = obs$lod_depth[i], range = rng,
               setting = names(pred), prediction = unname(pred))
      })
    }) %>%
    mutate(err = abs(pmax(prediction, 0) - AF) / range,
           z = ifelse(source == "imputed", NA_real_, (pmax(prediction, 0) - AF) / sqrt(AF * (1 - AF) / depth)))
}

# median effective df of each specification over the trajectories
spec_df <- function(calls, specs = spec_grid) {
  calls %>%
    group_split(traj_id) %>%
    map_dfr(function(obs) tibble(traj_id = obs$traj_id[1], setting = names(specs),
                                 df = map_dbl(specs, ~ fit_traj(obs$Age, obs$AF, .x)$df)))
}

# rule: among the specifications whose median hold-out error is within `tol` (relative) of the smallest, the smoothest
# (lowest median effective df)
select_spec <- function(calls, specs = spec_grid, tol = 0.02) {
  errors <- holdout_errors(calls, specs)
  dfs <- spec_df(calls, specs)
  crit <- errors %>%
    group_by(setting) %>%
    summarise(points = n(), median_err = median(err), mean_err = mean(err),
              median_abs_z = median(abs(z), na.rm = TRUE), .groups = "drop") %>%
    left_join(dfs %>% group_by(setting) %>% summarise(median_df = median(df), .groups = "drop"), by = "setting") %>%
    mutate(setting = factor(setting, levels = c(names(specs), "linear"))) %>%
    arrange(setting)
  cand <- crit %>% filter(setting != "linear", median_err <= (1 + tol) * min(median_err[setting != "linear"]))
  chosen <- as.character(cand$setting[which.min(cand$median_df)])
  list(errors = errors, df = dfs, criteria = crit, chosen = chosen)
}


# ---- Results the sensitivity analyses check ----------------------------------------------------------------------------

# gene-level mean of the per-trajectory mean speed, for the 5 genes
gene_speed <- function(sum) {
  sum %>%
    filter(Gene.refGene %in% main_genes) %>%
    group_by(Gene = as.character(Gene.refGene)) %>%
    summarise(n = n(), mean_VAF_speed_gene = mean(mean_VAF_speed), .groups = "drop")
}

gene_order <- function(sum) {
  g <- gene_speed(sum) %>% arrange(desc(mean_VAF_speed_gene))
  paste(g$Gene, collapse = " > ")
}

# Fabre et al. gene effect (largest GeneEffect_mean of the gene) against the SHCS gene mean speed: Pearson, one-sided
fabre_cor <- function(sum, fabre) {
  comp <- gene_speed(sum) %>%
    inner_join(fabre %>% group_by(Gene) %>% summarise(GeneEffect_mean = max(GeneEffect_mean), .groups = "drop"), by = "Gene")
  ct <- cor.test(comp$mean_VAF_speed_gene, comp$GeneEffect_mean, method = "pearson", alternative = "greater")
  list(r = unname(ct$estimate), p = ct$p.value, genes = comp)
}

adjusted_rand <- function(a, b) {
  t <- table(a, b)
  n <- sum(t)
  s <- sum(choose(t, 2))
  sa <- sum(choose(rowSums(t), 2))
  sb <- sum(choose(colSums(t), 2))
  e <- sa * sb / choose(n, 2)
  (s - e) / ((sa + sb) / 2 - e)
}

# z-scored VAF series of the curves (SHCS: ages 40-85; SardiNIA: the first 15 years of each curve); flat or
# one-point series cannot be z-scored and are left out
cluster_series <- function(pred, cohort = c("SHCS", "SardiNIA")) {
  cohort <- match.arg(cohort)
  d <- pred %>% arrange(traj_id, Age)
  d <- if (cohort == "SHCS") filter(d, Age >= 40, Age <= 85) else group_by(d, traj_id) %>% filter(Time <= 30) %>% ungroup()
  d <- d %>% group_by(traj_id) %>% filter(n() > 1, sd(VAF) > 0) %>% mutate(z = as.numeric(scale(VAF))) %>% ungroup()
  split(d$z, d$traj_id)
}

# one DTW clustering run (k = 2, DBA centroids); labels renamed to agree best with `ref` (traj_id, cluster)
cluster_dtw <- function(pred, ref = NULL) {
  series <- cluster_series(pred, "SHCS")
  pc <- dtwclust::tsclust(series, k = 2L, distance = "dtw_basic", centroid = "dba", seed = 8, norm = "L2",
                          window.size = 20L)
  cl <- tibble(traj_id = names(series), cluster = pc@cluster)
  if (!is.null(ref)) {
    j <- inner_join(cl, ref, by = "traj_id", suffix = c("", ".ref"))
    if (mean(j$cluster == j$cluster.ref) < 0.5) cl$cluster <- 3L - cl$cluster
  }
  list(model = pc, clusters = cl, swapped = !is.null(ref) && !identical(cl$cluster, pc@cluster))
}

# SardiNIA trajectories assigned to the nearest SHCS centroid (dtwclust predict); the labels follow those of cluster_dtw()
assign_sardinia <- function(cl, pred_sard) {
  series <- cluster_series(pred_sard, "SardiNIA")
  lab <- predict(cl$model, newdata = series)
  if (cl$swapped) lab <- 3L - lab
  tibble(traj_id = names(series), cluster = as.integer(lab))
}

# ---- Consensus DTW clustering ------------------------------------------------------------------------------------------
#
# One tsclust() run depends on the seed of the DBA centroids (different seeds give different partitions of the same
# curves). The consensus runs it with many seeds; M[i, j] = share of the runs in which trajectories i and j are in the
# same cluster (co-assignment), and the consensus partition is the average-linkage tree of 1 - M cut into k clusters
# (Monti et al. 2003). confidence = mean co-assignment of a trajectory with the other members of its consensus cluster.
# PAC = share of pairs with an ambiguous co-assignment (0.1-0.9); lower is more stable.

dtw_fit <- function(series, k, seed) {
  dtwclust::tsclust(series, k = k, distance = "dtw_basic", centroid = "dba", seed = seed, norm = "L2", window.size = 20L)
}

# label mapping lab -> target that maximizes the agreement (map[lab] gives the target label)
relabel_map <- function(lab, target, k) {
  as.integer(clue::solve_LSAP(table(factor(lab, 1:k), factor(target, 1:k)), maximum = TRUE))
}

consensus_series <- function(series, seeds = 1:100, k = 2L) {
  models <- lapply(seeds, function(s) dtw_fit(series, k, s))
  runs <- sapply(models, function(m) m@cluster)
  M <- matrix(0, length(series), length(series), dimnames = list(names(series), names(series)))
  for (j in seq_len(ncol(runs))) M <- M + outer(runs[, j], runs[, j], "==")
  M <- M / ncol(runs)
  cons <- cutree(hclust(as.dist(1 - M), method = "average"), k = k)
  cons <- match(cons, order(table(factor(cons, 1:k))))   # cluster 1 = the smallest cluster
  conf <- sapply(seq_along(cons), function(i) {
    same <- cons == cons[i]
    same[i] <- FALSE
    if (any(same)) mean(M[i, same]) else NA_real_
  })
  v <- M[lower.tri(M)]
  list(clusters = tibble(traj_id = names(series), cluster = cons, confidence = conf), M = M, runs = runs, models = models,
       maps = lapply(seq_len(ncol(runs)), function(j) relabel_map(runs[, j], cons, k)),
       pac = mean(v > 0.1 & v < 0.9), k = k)
}

# consensus clustering of the curves (SHCS series); labels renamed to agree best with `ref` (traj_id, cluster)
cluster_consensus <- function(pred, seeds = 1:100, k = 2L, ref = NULL) {
  cc <- consensus_series(cluster_series(pred, "SHCS"), seeds, k)
  if (!is.null(ref)) {
    j <- inner_join(cc$clusters, ref %>% select(traj_id, cluster), by = "traj_id", suffix = c("", ".ref"))
    m <- relabel_map(j$cluster, j$cluster.ref, k)
    cc$clusters$cluster <- m[cc$clusters$cluster]
    cc$maps <- lapply(cc$maps, function(x) m[x])
  }
  cc
}

# new series assigned to the nearest centroid of each run (dtwclust predict), labels mapped to the consensus,
# majority vote; share = fraction of the runs that agree with the vote
consensus_predict <- function(cc, series) {
  votes <- sapply(seq_along(cc$models), function(j) cc$maps[[j]][as.integer(predict(cc$models[[j]], newdata = series))])
  if (is.null(dim(votes))) votes <- matrix(votes, nrow = length(series))
  tibble(traj_id = names(series),
         cluster = apply(votes, 1, function(v) as.integer(names(which.max(table(v))))),
         share = apply(votes, 1, function(v) max(table(v)) / length(v)))
}

assign_sardinia_consensus <- function(cc, pred_sard) consensus_predict(cc, cluster_series(pred_sard, "SardiNIA"))


# SHCS vs SardiNIA, 5 genes: Wilcoxon test of the per-trajectory statistics, Fisher test of the cluster shares
cohort_stats <- c("mean_VAF_speed", "max_VAF_speed", "min_VAF_speed", "mean_VAF_acceleration", "max_VAF", "Age_zero", "Age_0.02")

cohort_tests <- function(sum_shcs, sum_sard, cl_shcs, cl_sard) {
  a <- sum_shcs %>% filter(Gene.refGene %in% main_genes)
  b <- sum_sard %>% filter(Gene.refGene %in% main_genes)
  tests <- map_dfr(cohort_stats, function(s) {
    x <- a[[s]][!is.na(a[[s]])]
    y <- b[[s]][!is.na(b[[s]])]
    tibble(stat = s, median_SHCS = median(x), median_SardiNIA = median(y), direction = sign(median(x) - median(y)),
           p = wilcox.test(x, y, exact = FALSE)$p.value)
  })
  ca <- cl_shcs %>% filter(traj_id %in% a$traj_id)
  cb <- cl_sard %>% filter(traj_id %in% b$traj_id)
  tab <- rbind(SHCS = table(factor(ca$cluster, 1:2)), SardiNIA = table(factor(cb$cluster, 1:2)))
  shares <- tibble(stat = "cluster 1 share", median_SHCS = tab[1, 1] / sum(tab[1, ]), median_SardiNIA = tab[2, 1] / sum(tab[2, ]),
                   direction = sign(tab[1, 1] / sum(tab[1, ]) - tab[2, 1] / sum(tab[2, ])), p = fisher.test(tab)$p.value)
  bind_rows(tests, shares)
}

# everything checked for one setting, SHCS (and SardiNIA when given): res = list(sum, clusters, sard_sum, sard_clusters)
key_stats <- c("mean_VAF_speed", "max_VAF_speed", "min_VAF_speed", "sd_VAF_speed",
               "mean_VAF_acceleration", "max_VAF_acceleration", "min_VAF_acceleration", "max_VAF",
               "Age_zero", "Age_LoD2", "Age_LoD", "Age_0.02")

setting_results <- function(res, fabre) {
  out <- list(sum = res$sum, clusters = res$clusters, gene_order = gene_order(res$sum), fabre = fabre_cor(res$sum, fabre))
  if (!is.null(res$sard_sum)) out$cohort <- cohort_tests(res$sum, res$sard_sum, res$clusters, res$sard_clusters)
  out
}

# `subset`: trajectories on which the per-trajectory agreement (rho, signs, groups) is computed (default: all shared)
compare_results <- function(x, ref, subset = NULL) {
  j <- inner_join(x$sum, ref$sum, by = "traj_id", suffix = c("", ".ref"))
  if (!is.null(subset)) j <- filter(j, traj_id %in% subset)
  rho <- map_dbl(key_stats, ~ suppressWarnings(cor(j[[.x]], j[[paste0(.x, ".ref")]], method = "spearman", use = "pairwise.complete.obs")))
  gx <- gene_speed(x$sum)
  gr <- gene_speed(ref$sum)
  gj <- inner_join(gx, gr, by = "Gene", suffix = c("", ".ref"))
  cj <- inner_join(x$clusters, ref$clusters, by = "traj_id", suffix = c("", ".ref"))
  row <- tibble(trajectories = nrow(j),
                `same sign, mean speed` = mean(sign(j$mean_VAF_speed) == sign(j$mean_VAF_speed.ref)),
                `same sign, mean acceleration` = mean(sign(j$mean_VAF_acceleration) == sign(j$mean_VAF_acceleration.ref)),
                `same dynamics group` = mean(j$group == j$group.ref),
                `gene order` = x$gene_order, `same gene order` = x$gene_order == ref$gene_order,
                `Kendall tau, gene order` = cor(gj$mean_VAF_speed_gene, gj$mean_VAF_speed_gene.ref, method = "kendall"),
                `Fabre r` = x$fabre$r, `Fabre p` = x$fabre$p,
                `clustered trajectories` = nrow(cj), `cluster ARI` = adjusted_rand(cj$cluster, cj$cluster.ref),
                `same cluster` = mean(cj$cluster == cj$cluster.ref))
  row <- bind_cols(row, as_tibble(as.list(setNames(rho, paste0("rho ", key_stats)))))
  if (!is.null(x$cohort) && !is.null(ref$cohort)) {
    cc <- inner_join(x$cohort, ref$cohort, by = "stat", suffix = c("", ".ref"))
    row$`cohort tests, same direction` <- sum(cc$direction == cc$direction.ref)
    row$`cohort tests, same significance` <- sum((cc$p < 0.05) == (cc$p.ref < 0.05))
    # same conclusion: same significance, and the same direction where the main analysis is significant
    row$`cohort tests, same conclusion` <- sum((cc$p < 0.05) == (cc$p.ref < 0.05) & (cc$p.ref >= 0.05 | cc$direction == cc$direction.ref))
    row$`cohort tests` <- nrow(cc)
  }
  row
}
