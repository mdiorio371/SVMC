# 09-rasr-phi.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- phi_from_abcd  [from ch4_functions.R:1619] ---
phi_from_abcd <- function(n11, n10, n01, n00) {
  a  <- as.numeric(n11); b <- as.numeric(n10); d01 <- as.numeric(n01); d00 <- as.numeric(n00)
  denom <- (a + b) * (d01 + d00) * (a + d01) * (b + d00)
  if (!is.finite(denom) || denom <= 0) return(0)
  (a * d00 - b * d01) / sqrt(denom)
}

# --- fisher_p_from_abcd  [from ch4_functions.R:1626] ---
fisher_p_from_abcd <- function(n11, n10, n01, n00) {
  m <- matrix(c(n11, n10, n01, n00), nrow = 2, byrow = TRUE)
  stats::fisher.test(m, alternative = "two.sided")$p.value
}

# --- score_event_pairs  [from ch4_functions.R:1632] ---
score_event_pairs <- function(P_sr, col_idx) {
  ev_ids <- colnames(P_sr)[col_idx]
  n_ev <- length(ev_ids); if (n_ev < 2) return(tibble())
  out <- vector("list", n_ev * (n_ev - 1) / 2); k <- 1L
  
  for (i in seq_len(n_ev - 1L)) {
    xi <- .col01(P_sr, col_idx[i])
    for (j in (i + 1L):n_ev) {
      yj <- .col01(P_sr, col_idx[j])
      n11 <- sum(xi == 1L & yj == 1L)
      n10 <- sum(xi == 1L & yj == 0L)
      n01 <- sum(xi == 0L & yj == 1L)
      n00 <- sum(xi == 0L & yj == 0L)
      supA <- n11 + n10; supB <- n11 + n01
      
      out[[k]] <- tibble(
        event_A   = ev_ids[i],
        event_B   = ev_ids[j],
        n11 = n11, n10 = n10, n01 = n01, n00 = n00,
        support_A = supA, support_B = supB,
        overlap   = n11,
        jaccard   = ifelse((n11 + n10 + n01) > 0, n11 / (n11 + n10 + n01), 0),
        phi       = phi_from_abcd(n11, n10, n01, n00),
        p_fisher  = fisher_p_from_abcd(n11, n10, n01, n00)
      )
      k <- k + 1L
    }
  }
  
  bind_rows(out) %>% mutate(q_fisher = p.adjust(p_fisher, "BH"))
}

# --- .safe_phi  [from ch4_functions.R:2104] ---
.safe_phi <- function(n11, n10, n01, n00){
  a <- as.numeric(n11); b <- as.numeric(n10); cc <- as.numeric(n01); d <- as.numeric(n00)
  num <- (a*d) - (b*cc); den <- sqrt((a+b)*(cc+d)*(a+cc)*(b+d))
  if (!is.finite(den) || den <= 0) return(0); num/den
}

# --- ._fisher_p  [from ch4_functions.R:2109] ---
._fisher_p <- function(n11, n10, n01, n00){
  m <- matrix(c(n11,n10,n01,n00),2,byrow=TRUE)
  if (any(rowSums(m)==0)||any(colSums(m)==0)) return(1)
  suppressWarnings(fisher.test(m)$p.value)
}

# --- classify_rasr  [from ch4_functions.R:1808] ---
classify_rasr <- function(
    P_sr, seeds,
    min_support_for_restricted = 5,
    max_outside_abs           = 2,
    max_outside_frac          = 0.05,
    min_inside_frac           = 0.80
){
  stopifnot(nrow(P_sr) > 0, ncol(P_sr) > 0)
  
  # identify seed columns A/B (if present)
  seedA <- if (length(seeds) >= 1) seeds[1] else NA_character_
  seedB <- if (length(seeds) >= 2) seeds[2] else NA_character_
  
  # 0/1 full matrix for fast colSums
  X <- .to01_matrix(P_sr)
  
  # groups from seeds (vectorized)
  a <- if (!is.na(seedA) && seedA %in% colnames(P_sr)) X[, match(seedA, colnames(P_sr))] else rep(0L, nrow(X))
  b <- if (!is.na(seedB) && seedB %in% colnames(P_sr)) X[, match(seedB, colnames(P_sr))] else rep(0L, nrow(X))
  groups <- ifelse(a==0 & b==0, "WT",
                   ifelse(a==1 & b==0, "A",
                          ifelse(a==0 & b==1, "B",
                                 "Mixed")))
  A_ids  <- which(groups=="A");  B_ids <- which(groups=="B")
  WT_ids <- which(groups=="WT");  Mx_ids <- which(groups=="Mixed")
  
  # per-event counts by branch (all vectorized via colSums)
  n_A     <- if (length(A_ids))  colSums(X[A_ids, , drop=FALSE])  else rep(0L, ncol(X))
  n_B     <- if (length(B_ids))  colSums(X[B_ids, , drop=FALSE])  else rep(0L, ncol(X))
  n_WT    <- if (length(WT_ids)) colSums(X[WT_ids,, drop=FALSE])  else rep(0L, ncol(X))
  n_Mixed <- if (length(Mx_ids)) colSums(X[Mx_ids,,drop=FALSE])   else rep(0L, ncol(X))
  n_total <- n_A + n_B + n_WT + n_Mixed
  
  # “best” branch per event and inside/outside tallies (vectorized)
  best_is_A   <- n_A >= n_B
  inside      <- ifelse(best_is_A, n_A, n_B)
  outside     <- ifelse(best_is_A, n_B, n_A) + n_WT + n_Mixed
  outside_cap <- pmax(max_outside_abs, ceiling(max_outside_frac * pmax(1L, n_total)))
  inside_frac <- ifelse(n_total > 0, inside / n_total, 0)
  
  # labels (vectorized)
  is_seed            <- colnames(P_sr) %in% seeds
  is_restricted      <- (n_total >= min_support_for_restricted) &
    (outside <= outside_cap) &
    (inside_frac >= min_inside_frac)
  label              <- ifelse(is_seed, "Branch-defining",
                               ifelse(is_restricted, "Branch-restricted", "Neutral"))
  restricted_branch  <- ifelse(label=="Branch-restricted",
                               ifelse(best_is_A, "A", "B"),
                               NA_character_)
  
  tibble(
    event_id         = colnames(P_sr),
    n_A              = as.integer(n_A),
    n_B              = as.integer(n_B),
    n_WT             = as.integer(n_WT),
    n_Mixed          = as.integer(n_Mixed),
    n_total          = as.integer(n_total),
    label            = label,
    restricted_branch= restricted_branch
  )
}

# --- identify_RASR_branches  [from ch4_functions.R:799] ---
identify_RASR_branches <- function(
    sv_rec,
    genome_len,
    out_dir = "sr_outputs",
    min_recurrence = 50,         # candidate seed recurrence
    top_k = 40,                  # candidate pool size
    max_branches  = 4,           # max seeds to return (A, B, C, …)
    exclusivity = list(          # seed exclusivity thresholds
      phi_max = -0.15,
      q_alpha = 0.05,
      overlap_frac_max = 0.03
    ),
    tier_criteria = list(        # labeling thresholds
      core_phi_min = 0.30, core_q_alpha = 0.01,
      core_p_e_s_min = 0.80, core_p_s_e_min = 0.80,
      exclto_phi_min = 0.20, exclto_q_alpha = 0.05,
      exclto_p_s_e_min = 0.90,
      exclto_other_phi_max = -0.05, exclto_other_q_min = 0.10,
      soft_phi_min = 0.15, soft_q_alpha = 0.10
    ),
    strata = NULL,               # optional named vector: names=sample ids
    plot_spokes = TRUE,
    save_files = TRUE,
    verbose = TRUE
){
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(purrr); library(Matrix)
    library(readr); library(stringr); library(ggplot2)
  })
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---------- helpers ---------------------------------------------------------
  phi_from_abcd <- function(a,b,c,d){
    a <- as.numeric(a); b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
    denom <- (a+b)*(c+d)*(a+c)*(b+d)
    if (!is.finite(denom) || denom <= 0) return(0)
    (a*d - b*c)/sqrt(denom)
  }
  fisher_p_from_abcd <- function(a,b,c,d){
    stats::fisher.test(matrix(c(a,b,c,d), 2, byrow=TRUE), alternative="two.sided")$p.value
  }
  score_event_pairs <- function(P_mat, col_idx){
    n_ev <- length(col_idx); if (n_ev < 2) return(tibble())
    ev_ids <- colnames(P_mat)[col_idx]
    out <- vector("list", n_ev*(n_ev-1)/2); k <- 1L
    for (i in seq_len(n_ev-1)) {
      xi <- as.integer(P_mat[, col_idx[i]] != 0)
      for (j in (i+1):n_ev) {
        yj <- as.integer(P_mat[, col_idx[j]] != 0)
        a <- sum(xi==1 & yj==1); b <- sum(xi==1 & yj==0)
        c <- sum(xi==0 & yj==1); d <- sum(xi==0 & yj==0)
        out[[k]] <- tibble(
          event_A = ev_ids[i], event_B = ev_ids[j],
          a=a,b=b,c=c,d=d,
          support_A = a+b, support_B = a+c,
          overlap = a,
          jaccard = ifelse((a+b+c)>0, a/(a+b+c), 0),
          phi = phi_from_abcd(a,b,c,d),
          p_fisher = fisher_p_from_abcd(a,b,c,d)
        )
        k <- k + 1L
      }
    }
    bind_rows(out) %>% mutate(q_fisher = p.adjust(p_fisher, "BH"))
  }

  compute_conditionals <- function(assoc_long){
    assoc_long %>%
      dplyr::mutate(
        p_event_given_seed = ifelse((a + c) > 0, a / (a + c), NA_real_),  # P(E|S)
        p_seed_given_event = ifelse((a + b) > 0, a / (a + b), NA_real_)   # P(S|E)
      )
  }
  
  assoc_vs_seeds <- function(P_sr, seeds){
    if (length(seeds) == 0) return(tibble())
    ev_ids <- colnames(P_sr)
    out <- purrr::map_dfr(ev_ids, function(eid){
      xi <- as.integer(P_sr[, eid] != 0)
      purrr::map_dfr(seeds, function(sid){
        yj <- as.integer(P_sr[, sid] != 0)
        a <- sum(xi==1 & yj==1); b <- sum(xi==1 & yj==0)
        c <- sum(xi==0 & yj==1); d <- sum(xi==0 & yj==0)
        tibble(
          event_id = eid, seed_id = sid,
          a=a,b=b,c=c,d=d,
          support_event = a+c, support_seed = a+b,
          overlap = a,
          phi = phi_from_abcd(a,b,c,d),
          p_fisher = fisher_p_from_abcd(a,b,c,d)
        )
      })
    })
    # BH within seed
    out %>%
      group_by(seed_id) %>% mutate(q_fisher = p.adjust(p_fisher, "BH")) %>% ungroup()
  }
  # seed discovery (agnostic, 0/1/N seeds)
  discover_branches <- function(P_sr, sr_meta, top_k, min_recurrence, max_branches, exclusivity, out_dir, verbose=TRUE){
    cand <- sr_meta %>% filter(n_samples >= min_recurrence) %>%
      arrange(desc(n_samples), desc(width_median)) %>% slice_head(n = top_k)
    keep_ids <- intersect(colnames(P_sr), cand$event_id)
    if (!length(keep_ids)) {
      msg <- sprintf("No SR meets min_recurrence ≥ %d.", min_recurrence)
      if (verbose) message(msg)
      return(list(seeds=character(0), candidates=cand, pair_scores=tibble(), seed_summary=msg))
    }
    ps <- score_event_pairs(P_sr, match(keep_ids, colnames(P_sr)))
    if (!nrow(ps)) {
      solo <- cand %>% filter(event_id %in% keep_ids) %>% slice_head(n=1)
      seeds <- if (nrow(solo)) solo$event_id else character(0)
      msg <- if (length(seeds)) paste0("Single robust SR: ", seeds) else "No viable seeds."
      if (verbose) message(msg)
      return(list(seeds=seeds, candidates=cand, pair_scores=ps, seed_summary=msg))
    }
    exclusivity_ok <- function(overlap, phi, q, supA, supB){
      rel_ok <- overlap <= exclusivity$overlap_frac_max * min(supA, supB)
      (phi <= exclusivity$phi_max) && (q <= exclusivity$q_alpha) && rel_ok
    }
    ps_excl <- ps %>%
      mutate(rel_overlap = overlap / pmax(pmin(support_A, support_B), 1),
             excl = mapply(exclusivity_ok, overlap, phi, q_fisher, support_A, support_B))
    
    deg <- ps_excl %>%
      dplyr::filter(excl) %>%
      dplyr::select(event_A, event_B) %>%
      tidyr::pivot_longer(cols = dplyr::everything(), values_to = "event") %>%
      dplyr::count(event, name = "excl_deg")
    
    ranker <- cand %>%
      filter(event_id %in% keep_ids) %>%
      left_join(deg, by = c("event_id"="event")) %>%
      mutate(excl_deg = tidyr::replace_na(excl_deg, 0L)) %>%
      arrange(desc(n_samples), desc(excl_deg))
    seeds <- character(0)
    for (ev in ranker$event_id) {
      if (!length(seeds)) {
        has_excl <- any(ps_excl$excl & (ps_excl$event_A==ev | ps_excl$event_B==ev))
        if (has_excl || (cand %>% filter(event_id==ev) %>% pull(n_samples)) >= min_recurrence*1.5) {
          seeds <- c(seeds, ev)
        }
      } else {
        ok_all <- TRUE
        for (s in seeds) {
          row <- ps_excl %>% filter((event_A==ev & event_B==s) | (event_A==s & event_B==ev))
          if (!nrow(row) || !row$excl[1]) { ok_all <- FALSE; break }
        }
        if (ok_all) seeds <- c(seeds, ev)
      }
      if (length(seeds) >= max_branches) break
    }
    msg <- if (!length(seeds)) "No mutually-exclusive seeds found."
    else if (length(seeds)==1) paste0("Found 1 seed: ", seeds)
    else paste0("Found ", length(seeds), " mutually-exclusive seeds:\n  - ", paste(seeds, collapse="\n  - "))
    if (verbose) message(msg)
    list(seeds=seeds, candidates=cand, pair_scores=ps_excl, seed_summary=msg)
  }
  # v2 tiered classifier (CoreCo / ExclusiveTo / Pref / Neutral)
  classify_tiered_v2 <- function(sr_meta, assoc_long, seeds, crit = tier_criteria){
    if (!length(seeds)) {
      return(list(
        tier_calls = sr_meta %>%
          transmute(event_id, n_samples, width_median, width_pct, label_v2 = "Neutral"),
        assoc_cond = tibble()
      ))
    }
    AL <- compute_conditionals(assoc_long) %>%
      filter(seed_id %in% seeds) %>%
      select(event_id, seed_id, phi, q_fisher, overlap, support_seed, support_event,
             p_event_given_seed, p_seed_given_event)
    
    # Core (symmetric, near-deterministic)
    core_flags <- AL %>%
      mutate(CoreCo = (phi >= crit$core_phi_min) &
               (q_fisher <= crit$core_q_alpha) &
               (p_event_given_seed >= crit$core_p_e_s_min) &
               (p_seed_given_event >= crit$core_p_s_e_min)) %>%
      filter(CoreCo) %>%
      group_by(event_id) %>% summarise(core_seeds = list(seed_id), n_core = n(), .groups="drop")
    
    # ExclusiveTo (asymmetric: inside a branch only)
    exclto_list <- lapply(seeds, function(s) {
      this <- AL %>% filter(seed_id==s) %>%
        mutate(ok_self = (phi >= crit$exclto_phi_min) &
                 (q_fisher <= crit$exclto_q_alpha) &
                 (p_seed_given_event >= crit$exclto_p_s_e_min)) %>%
        select(event_id, ok_self)
      others <- AL %>% filter(seed_id != s) %>%
        group_by(event_id) %>%
        summarise(ok_vs_others = all( (phi <= crit$exclto_other_phi_max) | (q_fisher >= crit$exclto_other_q_min) ),
                  .groups="drop")
      cand <- this %>% left_join(others, by="event_id") %>%
        mutate(ExclusiveTo_seed = s,
               ExclusiveTo = (tidyr::replace_na(ok_self, FALSE)) &
                 (tidyr::replace_na(ok_vs_others, TRUE)))
      cand %>% filter(ExclusiveTo) %>% select(event_id, ExclusiveTo_seed)
    })
    exclto_tbl <- bind_rows(exclto_list) %>%
      group_by(event_id) %>% summarise(exclto_seeds = list(unique(ExclusiveTo_seed)), n_exclto = n(), .groups="drop")
    
    # soft preference fallback
    soft_best <- AL %>%
      filter(phi >= crit$soft_phi_min, q_fisher <= crit$soft_q_alpha) %>%
      group_by(event_id) %>% slice_max(order_by = phi, n = 1, with_ties = FALSE) %>%
      ungroup() %>% select(event_id, soft_seed = seed_id, soft_phi = phi, soft_q = q_fisher)
    
    out <- sr_meta %>%
      select(event_id, n_samples, width_median, width_pct) %>%
      left_join(core_flags, by="event_id") %>%
      left_join(exclto_tbl, by="event_id") %>%
      left_join(soft_best, by="event_id") %>%
      mutate(n_core = tidyr::replace_na(n_core, 0L),
             n_exclto = tidyr::replace_na(n_exclto, 0L)) %>%
      rowwise() %>%
      mutate(
        label_v2 = dplyr::case_when(
          n_core >= 1 ~ paste0("CoreCo:", (core_seeds %||% list(NA_character_))[[1]][1]),
          n_exclto == 1 ~ paste0("ExclusiveTo:", (exclto_seeds %||% list(NA_character_))[[1]][1]),
          !is.na(soft_seed) ~ paste0("Pref:", soft_seed),
          TRUE ~ "Neutral"
        )
      ) %>% ungroup()
    
    list(tier_calls = out, assoc_cond = AL)
  }
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # readable seed aliases like A(0.857Mb), B(0.547Mb)
  make_seed_alias <- function(seeds, sr_meta){
    sizes <- sr_meta %>% filter(event_id %in% seeds) %>%
      mutate(sizeMb = width_median/1e6) %>%
      select(event_id, sizeMb) %>% arrange(match(event_id, seeds))
    alias <- sprintf("%s(%.3fMb)", LETTERS[seq_along(seeds)], sizes$sizeMb)
    names(alias) <- sizes$event_id
    alias
  }
  regex_escape <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  
  make_readable_seed_table <- function(tier_calls_v2, assoc_long, seeds, sr_meta, min_n = 0, out_csv = NULL){
    if (!length(seeds)) return(tier_calls_v2 %>% mutate(event_label=sprintf("%.3f Mb (n=%d)", width_median/1e6, n_samples)))
    seed_alias <- make_seed_alias(seeds, sr_meta)
    ev <- tier_calls_v2 %>%
      filter(n_samples >= min_n | event_id %in% seeds) %>%
      mutate(size_Mb = width_median/1e6,
             event_label = paste0(sprintf("%.3f", size_Mb), " Mb (n=", n_samples, ")")) %>%
      select(event_id, event_label, n_samples, size_Mb, label_v2)
    assoc_wide <- compute_conditionals(assoc_long) %>%
      filter(seed_id %in% seeds, event_id %in% ev$event_id, event_id != seed_id) %>%
      transmute(event_id, seed_id,
                phi = round(phi, 3),
                q   = signif(q_fisher, 3),
                P_E_given_S = round(p_event_given_seed, 3),
                P_S_given_E = round(p_seed_given_event, 3)) %>%
      pivot_wider(names_from = seed_id, values_from = c(phi, q, P_E_given_S, P_S_given_E), names_sep = "__")
    if (nrow(assoc_wide)) {
      cn <- colnames(assoc_wide)
      for (s in names(seed_alias)) {
        pat <- paste0("__", regex_escape(s), "$")
        repl <- paste0("__", seed_alias[[s]])
        hit <- grepl(pat, cn, perl = TRUE); cn[hit] <- sub(pat, repl, cn[hit], perl = TRUE)
      }
      colnames(assoc_wide) <- cn
    }
    out <- ev %>% left_join(assoc_wide, by="event_id") %>% arrange(desc(n_samples), desc(size_Mb))
    if (!is.null(out_csv)) { write_csv(out, out_csv); if (verbose) message("Wrote readable table -> ", out_csv) }
    out
  }
  
  # --- Drop-in replacement: preview on screen, optional saving -----------------
  plot_seed_spokes <- function(tier_calls_v2, seeds, P_sr, sr_meta,
                               out_png = NULL, out_pdf = NULL, save = FALSE) {
    stopifnot(is.data.frame(tier_calls_v2), is.character(seeds) || length(seeds)==0)
    
    # vectorized helper to extract the primary seed from label_v2
    parse_primary <- function(lbl_vec) {
      out <- rep(NA_character_, length(lbl_vec))
      ok  <- !is.na(lbl_vec) & grepl("^(CoreCo|ExclusiveTo|Pref):", lbl_vec)
      out[ok] <- sub("^(CoreCo|ExclusiveTo|Pref):", "", lbl_vec[ok])
      out
    }
    
    WT_n <- sum(Matrix::rowSums(P_sr) == 0)
    
    # Build edges from event -> its primary seed
    edges <- tier_calls_v2 %>%
      dplyr::mutate(primary_seed = parse_primary(label_v2),
                    tier = sub(":.*$", "", label_v2)) %>%
      dplyr::filter(!is.na(primary_seed), primary_seed %in% seeds,
                    event_id != primary_seed) %>%
      dplyr::select(event_id, n_samples, width_median, primary_seed, tier) %>%
      dplyr::mutate(size_Mb = width_median/1e6,
                    node_lab = sprintf("%d x %.3f Mb", n_samples, size_Mb)) %>%
      dplyr::filter(!is.na(primary_seed), primary_seed %in% seeds,
                    event_id != primary_seed)   
    
    # Seed labels
    seed_df <- sr_meta %>%
      dplyr::filter(event_id %in% seeds) %>%
      dplyr::arrange(match(event_id, seeds)) %>%
      dplyr::mutate(seed_lab = paste0(LETTERS[seq_along(seeds)],
                                      " (", sprintf("%.3f Mb", width_median/1e6), ")"))
    
    xs <- seq_along(seeds)
    nodes_seed <- data.frame(name = seeds, x = xs, y = 0, type = "seed",
                             label = seed_df$seed_lab, stringsAsFactors = FALSE)
    nodes_wt <- data.frame(name = "WT", x = ifelse(length(xs)>0, mean(xs), 0), y = -1.5,
                           type = "WT", label = paste0("WT: ", WT_n), stringsAsFactors = FALSE)
    
    # Events around each seed
    nodes_evt <- edges %>%
      dplyr::group_by(primary_seed) %>%
      dplyr::mutate(theta = seq(20, 340, length.out = dplyr::n()), r = 1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(x = xs[match(primary_seed, seeds)] + r*cospi(theta/180),
                    y = 0 + r*sinpi(theta/180),
                    name = event_id, type = tier, label = node_lab) %>%
      dplyr::select(name, x, y, type, label)
    
    nodes <- dplyr::bind_rows(nodes_seed, nodes_evt, nodes_wt)
    
    # Edge segments seed -> event
    e <- edges %>%
      dplyr::transmute(x = xs[match(primary_seed, seeds)], y = 0,
                       xend = nodes_evt$x[match(event_id, nodes_evt$name)],
                       yend = nodes_evt$y[match(event_id, nodes_evt$name)],
                       tier = tier)
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = e,
                            ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linetype = tier),
                            alpha = 0.8) +
      ggplot2::geom_point(data = subset(nodes, type %in% c("seed","WT")),
                          ggplot2::aes(x = x, y = y, color = type), size = 4) +
      ggplot2::geom_point(data = subset(nodes, !(type %in% c("seed","WT"))),
                          ggplot2::aes(x = x, y = y), size = 2, color = "grey30") +
      ggplot2::geom_text(data = nodes, ggplot2::aes(x = x, y = y, label = label),
                         vjust = -0.6, size = 3) +
      ggplot2::scale_color_manual(values = c(seed = "#1f78b4", WT = "#585858")) +
      ggplot2::scale_linetype_manual(values = c(CoreCo="solid", ExclusiveTo="solid", Pref="dotted", Neutral="blank")) +
      ggplot2::coord_equal() + ggplot2::theme_void() +
      ggplot2::ggtitle("RASR seed-centric view (solid=Core/Restricted; dotted=Associated; WT shown)")
    
    # Only save if asked to
    if (isTRUE(save)) {
      if (!is.null(out_png)) ggplot2::ggsave(out_png, p, width = 10, height = 6, dpi = 300)
      if (!is.null(out_pdf)) ggplot2::ggsave(out_pdf, p, width = 10, height = 6)
    }
    
    if (interactive()) print(p)     # show immediately in RStudio/interactive
    invisible(p)
  }
  
  
  # ---------- 1) Build presence, SR slice, metadata ---------------------------
  sv_tbl <- sv_rec %>%
    mutate(event_id = paste(variant, sv_group, sep="|"))
  samples <- sort(unique(sv_tbl$qid))
  events  <- sort(unique(sv_tbl$event_id))
  
  # presence_long (distinct pairs)
  presence_long <- sv_tbl %>%
    distinct(sample_id = qid, event_id) %>% mutate(present = 1L)
  
  # sparse presence matrix
  i <- match(presence_long$sample_id, samples)
  j <- match(presence_long$event_id,  events)
  P <- sparseMatrix(i=i, j=j, x=1L, dims=c(length(samples), length(events)),
                    dimnames = list(samples, events))
  
  # SR-only
  sr_ids <- sv_tbl %>% filter(variant == "Structural_rearrangement") %>% pull(event_id) %>% unique()
  P_sr <- if (length(sr_ids)) P[, intersect(colnames(P), sr_ids), drop=FALSE] else P[, 0, drop=FALSE]
  
  # sr_meta
  sr_meta <- sv_tbl %>%
    filter(event_id %in% colnames(P_sr)) %>%
    group_by(event_id, variant, sv_group) %>%
    summarise(n_samples = n_distinct(qid),
              width_median = median(width, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(width_pct = 100 * width_median / genome_len)
  
  if (verbose) {
    message(sprintf("Genome length used for width_pct: %s bp", scales::comma(genome_len)))
    if (nrow(sr_meta)) {
      message("Top Structural_rearrangement events by n_samples:")
      print(sr_meta %>% arrange(desc(n_samples)) %>% slice_head(n=20))
    }
  }
  
  # ---------- 2) Discover seeds ----------------------------------------------
  disc <- discover_branches(P_sr, sr_meta, top_k, min_recurrence, max_branches, exclusivity, out_dir, verbose)
  seeds <- disc$seeds
  
  # ---------- 3) Associations vs seeds + tier calls ---------------------------
  assoc_long <- assoc_vs_seeds(P_sr, seeds)
  cls <- classify_tiered_v2(sr_meta, assoc_long, seeds, tier_criteria)
  tier_calls_v2 <- cls$tier_calls
  
  # RASR labels (your 3-class scheme)
  tier_calls_v2 <- tier_calls_v2 %>%
    mutate(
      rasr_label = case_when(
        event_id %in% seeds ~ "Branch-defining",
        grepl("^CoreCo:", label_v2) & !(event_id %in% seeds) ~ "Branch-restricted",
        grepl("^ExclusiveTo:", label_v2) ~ "Branch-restricted",
        TRUE ~ "Neutral"
      )
    )
  
  # ---------- 4) Branch sample sets & summary --------------------------------
  WT <- rownames(P_sr)[rowSums(P_sr) == 0]
  seed_cols <- if (length(seeds)) P_sr[, seeds, drop=FALSE] else Matrix(0, nrow(P_sr), 0)
  by_seed <- lapply(seq_along(seeds), function(k){
    s <- seeds[k]
    mine <- rownames(P_sr)[which(seed_cols[, s] != 0)]
    others <- if (ncol(seed_cols) > 1) rownames(P_sr)[which(rowSums(seed_cols[, setdiff(seeds, s), drop=FALSE]) != 0)] else character(0)
    setdiff(mine, others)
  })
  names(by_seed) <- seeds
  Mixed <- if (length(seeds) >= 2) rownames(P_sr)[rowSums(seed_cols) >= 2] else character(0)
  
  sets <- c(list(WT=WT), by_seed, list(Mixed=Mixed))
  total_n <- nrow(P_sr)
  summary_tbl <- tibble(
    set = c("WT", seeds, "Mixed"),
    n = c(length(WT), sapply(by_seed, length), length(Mixed)),
    pct = round(100 * n / total_n, 1)
  )
  
  # ---------- 5) Readable table ----------------------------------------------
  readable_tbl <- make_readable_seed_table(tier_calls_v2, assoc_long, seeds, sr_meta, min_n = 0,
                                           out_csv = if (save_files) file.path(out_dir,"sr_readable_seed_table.csv") else NULL)
  
  # ---------- 6) Exclusion summary among seeds -------------------------------
  seed_pairs <- if (length(seeds) >= 2) {
    ps <- score_event_pairs(P_sr, match(seeds, colnames(P_sr))) %>%
      filter((event_A %in% seeds) & (event_B %in% seeds)) %>%
      arrange(phi, q_fisher)
    if (save_files) write_csv(ps, file.path(out_dir,"seed_exclusivity_pairs.csv"))
    ps
  } else tibble()
  
  # ---------- 7) Plot: seed-centric spokes with WT ---------------------------
  plot_paths <- list()
  if (plot_spokes) {
    p_png <- file.path(out_dir, "rasr_seed_spokes.png")
    p_pdf <- file.path(out_dir, "rasr_seed_spokes.pdf")
    plot_seed_spokes(tier_calls_v2, seeds, P_sr, sr_meta, p_png, p_pdf)
    plot_paths <- list(spokes_png = p_png, spokes_pdf = p_pdf)
  }
  
  # ---------- 8) Save core files ---------------------------------------------
  if (save_files) {
    saveRDS(P,      file.path(out_dir, "P.rds"))
    saveRDS(P_sr,   file.path(out_dir, "P_sr.rds"))
    saveRDS(presence_long, file.path(out_dir, "presence_long.rds"))
    write_csv(sr_meta, file.path(out_dir, "sr_meta.csv"))
    write_lines(disc$seed_summary, file.path(out_dir, "seed_summary.txt"))
    write_csv(tier_calls_v2 %>% select(event_id, n_samples, width_median, width_pct, label_v2, rasr_label),
              file.path(out_dir, "branch_calls_rasr_labels.csv"))
    if (nrow(assoc_long)) write_csv(assoc_long, file.path(out_dir, "sr_seed_associations_long.csv"))
    write_csv(readable_tbl, file.path(out_dir, "sr_readable_seed_table.csv"))
    write_csv(summary_tbl,  file.path(out_dir, "branch_samples.csv"))
  }
  
  
  labels3 <- classify_rasr(P_sr, seeds)   
  
  
  readable_table <- sr_meta %>%
    dplyr::transmute(event_id, n_samples, size_Mb = width_median / 1e6) %>%
    dplyr::left_join(labels3 %>% dplyr::rename(rasr_label = label), by = "event_id")
  
  
  if (isTRUE(save_files)) {
    readr::write_csv(readable_table, file.path(out_dir, "sr_readable_seed_table.csv"))
    readr::write_csv(labels3,       file.path(out_dir, "sr_labels_threeway.csv"))
  }
  return(list(
    seeds = seeds,
    P_sr = P_sr,
    sr_meta = sr_meta,
    readable_table = readable_table,
    labels_threeway = labels3,
    seed_exclusivity = disc$pair_scores,
    seed_summary = disc$seed_summary,
    P = P, P_sr = P_sr, presence_long = presence_long,
    sr_meta = sr_meta,
    assoc_long = assoc_long,
    tier_calls = tier_calls_v2,          # includes label_v2 and rasr_label
    readable_table = readable_tbl,
    branch_sets = list(sets = sets, summary = summary_tbl),
    plot_paths = plot_paths,
    params = list(
      min_recurrence=min_recurrence, top_k=top_k, max_branches=max_branches,
      exclusivity=exclusivity, tier_criteria=tier_criteria
    )
  ))

}

# --- rasr_seed_discovery  [from ch4_functions.R:1677] ---
rasr_seed_discovery <- function(
    P_sr, sr_meta,
    out_dir            = NULL,     # e.g., "sr_outputs"
    top_k              = 40,       # candidate pool size
    min_recurrence     = 50,       # minimum n_samples to consider as candidates
    max_branches       = 4,        # cap on number of mutually-exclusive seeds to return
    phi_max            = -0.15,    # exclusivity strength (more negative = stricter)
    q_alpha            = 0.05,     # BH-adjusted p-value threshold
    overlap_frac_max   = 0.03,     # allowed overlap as fraction of min(supports)
    verbose            = TRUE
) {
  if (!is.null(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1) Candidate RASRs by recurrence and present in matrix
  cand <- sr_meta %>%
    filter(n_samples >= min_recurrence, event_id %in% colnames(P_sr)) %>%
    arrange(desc(n_samples), desc(width_median)) %>%
    slice_head(n = top_k)
  
  if (!nrow(cand)) {
    msg <- sprintf("No SRs meet min_recurrence ≥ %d. Returning 0 seeds.", min_recurrence)
    if (verbose) message(msg)
    if (!is.null(out_dir)) writeLines(msg, file.path(out_dir, "seed_summary.txt"))
    groups <- setNames(rep("WT", nrow(P_sr)), rownames(P_sr))
    return(list(seeds = character(0), pair_scores = tibble(), groups = groups,
                group_counts = tibble(group="WT", n=length(groups), pct=100), summary = msg))
  }
  
  keep_idx <- match(cand$event_id, colnames(P_sr))
  keep_idx <- keep_idx[!is.na(keep_idx)]
  if (!length(keep_idx)) {
    msg <- "Candidate IDs not present in P_sr. Returning 0 seeds."
    if (verbose) message(msg)
    if (!is.null(out_dir)) writeLines(msg, file.path(out_dir, "seed_summary.txt"))
    groups <- setNames(rep("WT", nrow(P_sr)), rownames(P_sr))
    return(list(seeds = character(0), pair_scores = tibble(), groups = groups,
                group_counts = tibble(group="WT", n=length(groups), pct=100), summary = msg))
  }
  
  # 2) Pairwise exclusivity within candidates
  ps <- score_event_pairs(P_sr, keep_idx)
  
  if (!nrow(ps)) {
    # Single strong candidate path
    seeds <- cand$event_id[1]
    msg <- paste0("Only one candidate available; returning 1 seed: ", seeds)
    if (verbose) message(msg)
    if (!is.null(out_dir)) writeLines(msg, file.path(out_dir, "seed_summary.txt"))
    groups <- branch_groups_from_seeds(P_sr, seeds[1], NA_character_)
    tab <- table(groups)
    group_counts <- tibble(group = names(tab), n = as.integer(tab)) %>%
      mutate(pct = round(100 * n / sum(n), 1)) %>%
      arrange(match(group, c("WT","A","B","Mixed")))
    return(list(seeds = seeds, pair_scores = ps, groups = groups,
                group_counts = group_counts, summary = msg))
  }
  
  # exclusivity condition
  exclusivity_ok <- function(overlap, phi, q, supA, supB) {
    rel_ok <- overlap <= overlap_frac_max * min(supA, supB)
    (phi <= phi_max) && (q <= q_alpha) && rel_ok
  }
  
  ps_excl <- ps %>%
    mutate(rel_overlap = overlap / pmax(pmin(support_A, support_B), 1),
           excl = mapply(exclusivity_ok, overlap, phi, q_fisher, support_A, support_B))
  
  # 3) Greedy selection of a mutually-exclusive set of seeds
  deg <- ps_excl %>% filter(excl) %>%
    select(event_A, event_B) %>%
    pivot_longer(cols = everything(), values_to = "event") %>%
    count(event, name = "excl_deg")
  
  ranker <- cand %>%
    left_join(deg, by = c("event_id" = "event")) %>%
    mutate(excl_deg = tidyr::replace_na(excl_deg, 0L)) %>%
    arrange(desc(n_samples), desc(excl_deg))
  
  seeds <- character(0)
  for (ev in ranker$event_id) {
    if (!length(seeds)) {
      has_excl <- any(ps_excl$excl & (ps_excl$event_A == ev | ps_excl$event_B == ev))
      if (has_excl || (cand %>% filter(event_id == ev) %>% pull(n_samples) >= min_recurrence * 1.5)) {
        seeds <- c(seeds, ev)
      }
    } else {
      ok_all <- TRUE
      for (s in seeds) {
        row <- ps_excl %>% filter((event_A == ev & event_B == s) | (event_A == s & event_B == ev))
        if (!nrow(row) || !row$excl[1]) { ok_all <- FALSE; break }
      }
      if (ok_all) seeds <- c(seeds, ev)
    }
    if (length(seeds) >= max_branches) break
  }
  
  # 4) Summary + branch labels (first two seeds → A/B)
  msg <- if (!length(seeds)) {
    "No mutually-exclusive SR seeds found under thresholds. Returning 0 seeds."
  } else if (length(seeds) == 1) {
    paste0("Found a single robust SR seed (no exclusive partner under thresholds):\n  - ", seeds)
  } else {
    paste0("Found ", length(seeds), " mutually-exclusive seeds:\n  - ", paste(seeds, collapse = "\n  - "))
  }
  if (verbose) message(msg)
  if (!is.null(out_dir)) writeLines(msg, file.path(out_dir, "seed_summary.txt"))
  
  seedA <- if (length(seeds) >= 1) seeds[1] else NA_character_
  seedB <- if (length(seeds) >= 2) seeds[2] else NA_character_
  groups <- branch_groups_from_seeds(P_sr, seedA, seedB)
  
  tab <- table(groups)
  group_counts <- tibble(group = names(tab), n = as.integer(tab)) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    arrange(match(group, c("WT","A","B","Mixed")))
  
  list(
    seeds        = seeds,
    pair_scores  = ps_excl,
    groups       = groups,
    group_counts = group_counts,
    summary      = msg,
    params       = list(top_k=top_k, min_recurrence=min_recurrence, max_branches=max_branches,
                        phi_max=phi_max, q_alpha=q_alpha, overlap_frac_max=overlap_frac_max)
  )
}

# --- rasr_pairwise_simple  [from ch4_functions.R:1890] ---
rasr_pairwise_simple <- function(P_sr, events = NULL, psr_ids = NULL, q_alpha = 0.05) {
  stopifnot(is.matrix(P_sr) || inherits(P_sr, "Matrix"))
  if (is.null(colnames(P_sr))) stop("P_sr must have colnames=event_ids.")
  if (is.null(rownames(P_sr))) stop("P_sr must have rownames=sample_ids.")
  
  ids <- if (!is.null(events)) as.character(events) else if (!is.null(psr_ids)) as.character(psr_ids) else colnames(P_sr)
  ids <- intersect(ids, colnames(P_sr))
  k <- length(ids)
  
  empty_pairs <- tibble(
    event_A = character(), event_B = character(),
    n11 = integer(), n10 = integer(), n01 = integer(), n00 = integer(),
    phi = numeric(), p = numeric(), q = numeric(), jaccard = numeric(),
    pair_label = character()
  )
  empty_cooc <- matrix(0L, nrow = k, ncol = k, dimnames = list(ids, ids))
  
  if (k == 0) {
    per_event <- tibble(event_id = character(), label = character(), n_pairs = integer())
    return(list(pairs = empty_pairs, per_event = per_event, cooc_matrix = empty_cooc))
  }
  if (k == 1) {
    per_event <- tibble(event_id = ids, label = "Neutral", n_pairs = 0L)
    return(list(pairs = empty_pairs, per_event = per_event, cooc_matrix = empty_cooc))
  }
  
  X <- (P_sr[, ids, drop = FALSE] != 0)
  n <- nrow(X)
  n1 <- Matrix::colSums(X)
  n11_mat <- Matrix::crossprod(X)
  
  idx <- which(upper.tri(matrix(0, k, k)), arr.ind = TRUE)
  i_idx <- idx[,1]; j_idx <- idx[,2]
  
  pairs <- purrr::map_dfr(seq_along(i_idx), function(t) {
    i <- i_idx[t]; j <- j_idx[t]
    a <- as.numeric(n11_mat[i, j])
    b <- as.numeric(n1[i] - a)
    cc <- as.numeric(n1[j] - a)
    d <- as.numeric(n - a - b - cc)
    tibble(
      event_A = ids[i], event_B = ids[j],
      n11 = a, n10 = b, n01 = cc, n00 = d,
      phi = .safe_phi(a, b, cc, d),
      p   = ._fisher_p(a, b, cc, d),
      jaccard = if ((a + b + cc) > 0) a / (a + b + cc) else 0
    )
  }) %>% mutate(q = p.adjust(p, method = "BH")) %>%
    mutate(pair_label = case_when(
      (n11 == 0) ~ "Exclusive",
      (phi > 0 & q <= q_alpha) ~ "Co-occurrant",
      (phi < 0 & q <= q_alpha) ~ "Exclusive",
      TRUE ~ "Neutral"
    ))
  
  per_event <- bind_rows(
    pairs %>% transmute(event_id = event_A, pair_label),
    pairs %>% transmute(event_id = event_B, pair_label)
  ) %>%
    group_by(event_id) %>%
    summarise(
      n_pairs = n(),
      n_cooc  = sum(pair_label == "Co-occurrant"),
      n_excl  = sum(pair_label == "Exclusive"),
      label = case_when(
        n_cooc > 0 ~ "Co-occurrant",
        n_cooc == 0 & n_excl > 0 ~ "Exclusive",
        TRUE ~ "Neutral"
      ),
      .groups = "drop"
    )
  
  cooc_matrix <- matrix(0L, nrow = k, ncol = k, dimnames = list(ids, ids))
  cooc_matrix[upper.tri(cooc_matrix)] <- as.integer(pairs$n11 > 0)
  cooc_matrix <- cooc_matrix + t(cooc_matrix); diag(cooc_matrix) <- 0L
  
  list(pairs = pairs, per_event = per_event, cooc_matrix = cooc_matrix)
}

# --- rasr_scoreboard  [from ch4_functions.R:2117] ---
rasr_scoreboard <- function(P_sr_group, sr_meta, events, q_alpha=0.05,
                            strict_exclusive=TRUE, allow_partial=TRUE) {
  events <- trimws(as.character(events))
  cols   <- colnames(P_sr_group)
  missing <- setdiff(events, cols)
  if (length(missing)) {
    msg <- paste0("These event_id(s) not found in P_sr_group and will be dropped: ",
                  paste(missing, collapse=", "))
    if (allow_partial) warning(msg) else stop(msg)
  }
  ids <- intersect(events, cols)
  if (length(ids) < 2L) stop("Fewer than 2 requested events are present after matching.")
  
  X <- (P_sr_group[, ids, drop=FALSE] != 0)
  n <- nrow(X); k <- ncol(X)
  
  # label map
  lab <- sr_meta %>%
    transmute(event_id = trimws(as.character(event_id)),
              RASR_id  = ifelse(is.na(width_median), NA_character_,
                                paste0("RASR_", sprintf("%.2f", width_median/1e6), "Mbp"))) %>%
    filter(event_id %in% ids) %>% { setNames(.$RASR_id, .$event_id) }
  
  # counts
  n1  <- Matrix::colSums(X)
  n11 <- Matrix::crossprod(X)  # k x k
  
  # pairwise stats
  pairs <- tibble::tibble()
  for (i in 1:(k-1)) for (j in (i+1):k) {
    a <- as.numeric(n11[i,j]); b <- as.numeric(n1[i]-a); cc <- as.numeric(n1[j]-a)
    d <- as.numeric(n - a - b - cc)
    pairs <- dplyr::bind_rows(pairs, tibble::tibble(
      event_A = ids[i], event_B = ids[j],
      n11=a, n10=b, n01=cc, n00=d,
      phi = {
        num <- (a*d) - (b*cc); den <- sqrt((a+b)*(cc+d)*(a+cc)*(b+d))
        if (!is.finite(den) || den <= 0) 0 else num/den
      },
      p = {
        M <- matrix(c(a,b,cc,d),2,byrow=TRUE)
        if (any(rowSums(M)==0) || any(colSums(M)==0)) 1 else suppressWarnings(fisher.test(M)$p.value)
      }
    ))
  }
  if (nrow(pairs)) pairs$q <- p.adjust(pairs$p, "BH")
  
  # logical adjacency: observed, exclusive, significant co-occurrent
  OBS <- (n11 > 0)
  E   <- matrix(FALSE, k, k, dimnames=list(ids, ids))
  C   <- matrix(FALSE, k, k, dimnames=list(ids, ids))
  
  if (nrow(pairs)) {
    e_idx <- with(pairs, if (strict_exclusive) which(n11==0)
                  else which((n11==0) | (phi<0 & q<=q_alpha)))
    if (length(e_idx)) {
      A <- cbind(match(pairs$event_A[e_idx], ids), match(pairs$event_B[e_idx], ids))
      E[A] <- E[A[,2:1]] <- TRUE
    }
    c_idx <- with(pairs, which(phi>0 & q<=q_alpha & n11>0))
    if (length(c_idx)) {
      A <- cbind(match(pairs$event_A[c_idx], ids), match(pairs$event_B[c_idx], ids))
      C[A] <- C[A[,2:1]] <- TRUE
    }
  }
  diag(OBS) <- diag(E) <- diag(C) <- FALSE
  
  # occurs-alone (within this selected set)
  occurs_alone <- sapply(seq_len(k), function(i) any(X[,i] & rowSums(X[,-i, drop=FALSE])==0))
  names(occurs_alone) <- ids
  
  # list-of-elements form
  out_list <- lapply(ids, function(e) list(
    event_id = e,
    RASR_id  = lab[[e]],
    occurs_alone     = occurs_alone[[e]],
    observed_with    = unname(lab[ids[OBS[e, ]]]),            # NEW
    exclusive_to     = unname(lab[ids[E[e, ]]]),
    cooccurrent_with = unname(lab[ids[C[e, ]]])
  ))
  names(out_list) <- unname(lab[ids])
  
  # tidy table form
  out_tbl <- tibble::tibble(
    event_id = ids,
    RASR_id  = unname(lab[ids]),
    occurs_alone = occurs_alone,
    observed_with    = sapply(ids, function(e) paste(unname(lab[ids[OBS[e, ]]]), collapse=", ")), # NEW
    exclusive_to     = sapply(ids, function(e) paste(unname(lab[ids[E[e, ]]]),   collapse=", ")),
    cooccurrent_with = sapply(ids, function(e) paste(unname(lab[ids[C[e, ]]]),   collapse=", "))
  )
  
  list(list = out_list, table = out_tbl, pairs = pairs)
}

# --- rasr_presence_table  [from ch4_functions.R:2253] ---
rasr_presence_table <- function(P_sr_group, sr_meta, events, use_labels = TRUE, digits = 2) {
  if (is.null(colnames(P_sr_group)) || is.null(rownames(P_sr_group))) {
    stop("P_sr_group must have both rownames (accessions) and colnames (event_id).")
  }
  events <- as.character(events)
  present <- intersect(events, colnames(P_sr_group))
  if (length(present) == 0) stop("None of the requested events are columns in P_sr_group.")
  if (length(present) < length(events)) {
    warning("Dropping missing events: ", paste(setdiff(events, present), collapse = ", "))
  }
  
  # label map: event_id -> RASR label
  lab_map <- sr_meta %>%
    filter(event_id %in% present) %>%
    transmute(event_id, label = rasr_make_label(width_median, digits)) %>%
    { setNames(.$label, .$event_id) }
  # fallback to event_id when label is NA
  for (e in present) if (is.na(lab_map[[e]])) lab_map[[e]] <- e
  
  X <- (P_sr_group[, present, drop = FALSE] != 0)
  df <- as_tibble(as.matrix(X)) %>%
    mutate(accession = rownames(X), .before = 1) %>%
    mutate(across(-accession, ~ as.integer(.x)))
  
  if (use_labels) {
    new_names <- c("accession", unname(lab_map[present]))
    colnames(df) <- new_names
  } else {
    colnames(df) <- c("accession", present)
  }
  df %>% arrange(accession)
}

# --- build_P_sr_group  [from ch4_functions.R:2215] ---
build_P_sr_group <- function(sv_rec, variant = "Structural_rearrangement") {
  stopifnot(all(c("variant","sv_group","qid") %in% names(sv_rec)))
  df <- sv_rec %>%
    mutate(
      variant_chr   = trimws(as.character(.data$variant)),
      sv_group_chr  = trimws(as.character(.data$sv_group)),
      sample        = as.character(.data$qid)
    ) %>%
    filter(variant_chr == !!variant) %>%
    transmute(event_id = paste(variant_chr, sv_group_chr, sep="|"),
              sample) %>%
    distinct()
  
  samp <- sort(unique(df$sample))
  evts <- sort(unique(df$event_id))
  sparseMatrix(i = match(df$sample, samp),
               j = match(df$event_id, evts),
               x = 1L,
               dims = c(length(samp), length(evts)),
               dimnames = list(samp, evts))
}

# --- scan_one_feature  [from SVCM_result1_helpers.R:335] ---
scan_one_feature <- function(feature_id, M_bin, Dmat, cluster_raw,
                             nperm = 199L) {
  carriers_idx     <- which(M_bin[, feature_id] != 0)
  noncarriers_idx  <- which(M_bin[, feature_id] == 0)
  n_car <- length(carriers_idx)
  n_non <- length(noncarriers_idx)

  empty_row <- tibble::tibble(
    feature_id                       = feature_id,
    carriers_overlap                 = as.integer(n_car),
    mean_within                      = NA_real_,
    mean_between                     = NA_real_,
    delta                            = NA_real_,
    delta_z                          = NA_real_,
    pct_dom_nonnoise                 = NA_real_,
    outside_dom_carriers_nonnoise    = NA_integer_
  )
  if (n_car < 3L || n_non < 3L) return(empty_row)

  mw <- mean_within_dist(carriers_idx, Dmat)
  mb <- mean_between_dist(carriers_idx, noncarriers_idx, Dmat)
  delta_obs <- mb - mw

  N <- nrow(M_bin)
  perm_delta <- vapply(seq_len(nperm), function(.k) {
    sampled <- sample.int(N, n_car)
    mw_p <- mean_within_dist(sampled, Dmat)
    mb_p <- mean_between_dist(sampled, setdiff(seq_len(N), sampled), Dmat)
    mb_p - mw_p
  }, FUN.VALUE = numeric(1))

  perm_mean <- mean(perm_delta, na.rm = TRUE)
  perm_sd   <- stats::sd(perm_delta, na.rm = TRUE)
  delta_z <- if (is.finite(perm_sd) && perm_sd > 0) {
    (delta_obs - perm_mean) / perm_sd
  } else {
    NA_real_
  }

  # Cluster purity of carriers across non-noise clusters
  carrier_clusters <- as.character(cluster_raw[carriers_idx])
  non_noise <- carrier_clusters[
    !is.na(carrier_clusters) &
    !carrier_clusters %in% c("0", "Noise", "")
  ]
  if (length(non_noise) > 0L) {
    tab <- sort(table(non_noise), decreasing = TRUE)
    pct_dom     <- as.numeric(tab[1]) / length(non_noise)
    outside_dom <- length(non_noise) - as.integer(tab[1])
  } else {
    pct_dom     <- NA_real_
    outside_dom <- NA_integer_
  }

  tibble::tibble(
    feature_id                    = feature_id,
    carriers_overlap              = as.integer(n_car),
    mean_within                   = mw,
    mean_between                  = mb,
    delta                         = delta_obs,
    delta_z                       = delta_z,
    pct_dom_nonnoise              = pct_dom,
    outside_dom_carriers_nonnoise = outside_dom
  )
}

# --- assign_states  [from SVCM_result1_helpers.R:405] ---
assign_states <- function(M_mod, major_ids, label_map = NULL) {
  if (is.null(label_map)) {
    label_map <- stats::setNames(major_ids, major_ids)
  }
  if (!all(major_ids %in% colnames(M_mod))) {
    missing <- setdiff(major_ids, colnames(M_mod))
    stop("major_ids not in M_mod columns: ", paste(missing, collapse = ", "))
  }
  M_sub <- M_mod[, major_ids, drop = FALSE]
  presence <- rowSums(M_sub != 0)
  qids <- rownames(M_mod)

  state <- vapply(seq_along(qids), function(i) {
    np <- presence[i]
    if (np == 0L) return("None")
    if (np > 1L)  return("Multiple")
    present_idx <- which(M_sub[i, ] != 0)[1]
    unname(label_map[[major_ids[present_idx]]])
  }, FUN.VALUE = character(1))

  n_multiple <- sum(state == "Multiple")
  if (n_multiple > 0L) {
    warning(sprintf(
      "%d genomes carry multiple major motifs ('Multiple' state). Selected modules may not be fully mutually exclusive.",
      n_multiple
    ))
  }

  tibble::tibble(qid = qids, state = state)
}

# --- .col01  [from ch4_functions.R:1610] ---
.col01 <- function(M, colname_or_index) {
  j <- if (is.numeric(colname_or_index)) colname_or_index else match(colname_or_index, colnames(M))
  if (is.na(j)) return(rep(0L, nrow(M)))
  v <- if (inherits(M, "sparseMatrix")) Matrix::as.matrix(M[, j, drop = FALSE])[, 1] else M[, j]
  v[is.na(v)] <- 0L
  as.integer(v != 0L)
}

# --- .branch_groups_from_seeds  [from ch4_functions.R:737] ---
.branch_groups_from_seeds <- function(P, seedA, seedB) {
  a <- if (!is.na(seedA) && seedA %in% colnames(P)) .col01(P, seedA) else rep(0L, nrow(P))
  b <- if (!is.na(seedB) && seedB %in% colnames(P)) .col01(P, seedB) else rep(0L, nrow(P))
  dplyr::case_when(
    a == 0 & b == 0 ~ "WT",
    a == 1 & b == 0 ~ "A",
    a == 0 & b == 1 ~ "B",
    a == 1 & b == 1 ~ "Mixed",
    TRUE            ~ "WT"
  )
}

# --- branch_groups_from_seeds  [from ch4_functions.R:1665] ---
branch_groups_from_seeds <- function(P, seedA, seedB) {
  if (is.null(rownames(P))) rownames(P) <- paste0("sample_", seq_len(nrow(P)))
  a <- if (!is.na(seedA) && seedA %in% colnames(P)) .col01(P, seedA) else rep(0L, nrow(P))
  b <- if (!is.na(seedB) && seedB %in% colnames(P)) .col01(P, seedB) else rep(0L, nrow(P))
  grp <- ifelse(a==0 & b==0, "WT",
                ifelse(a==1 & b==0, "A",
                       ifelse(a==0 & b==1, "B", "Mixed")))
  names(grp) <- rownames(P)
  grp
}

# --- rasr_make_label  [from ch4_functions.R:2244] ---
rasr_make_label <- function(width_bp, digits = 2) {
  ifelse(is.na(width_bp),
         NA_character_,
         paste0("RASR_", sprintf(paste0("%.", digits, "f"), width_bp / 1e6), "Mbp"))
}

# --- rasr_prelim_inputs  [from ch4_functions.R:1336] ---
rasr_prelim_inputs <- function(
    sv_rec,
    genome_len,                        # REQUIRED (e.g., reflen)
    variant_name = "Structural_rearrangement",
    out_dir     = NULL,                # e.g., "sr_outputs" or NULL to skip saving
    save_files  = FALSE,
    verbose     = TRUE
){
  # 1) minimally required columns
  need <- c("variant","sv_group","qid")
  miss <- setdiff(need, names(sv_rec))
  if (length(miss)) stop("sv_rec is missing required columns: ", paste(miss, collapse=", "))
  
  # 2) event_id and width
  df <- sv_rec %>%
    mutate(
      variant   = as.character(variant),
      sv_group  = as.character(sv_group),
      qid       = as.character(qid),
      event_id  = paste0(variant, "|", sv_group),
      width     = .compute_width(cur_data())
    )
  
  # 3) presence_long & full P
  presence_long <- df %>%
    distinct(qid, event_id) %>%
    mutate(present = 1L)
  
  P <- .make_presence_matrix(presence_long, sample_col = "qid", event_col = "event_id")
  
  # 4) SR-only slice (by event_id prefix to keep it agnostic)
  sr_cols <- grepl(paste0("^", stringr::fixed(variant_name), "\\|"), colnames(P))
  P_sr <- if (any(sr_cols)) P[, sr_cols, drop = FALSE] else {
    if (verbose) message("No columns matched variant_name='", variant_name, "'. Returning 0-column P_sr.")
    P[, 0]
  }
  
  # 5) per-event SR metadata
  sr_meta <- df %>%
    filter(variant == !!variant_name) %>%
    group_by(event_id, variant, sv_group) %>%
    summarise(
      n_samples    = n_distinct(qid),
      width_median = median(as.numeric(width), na.rm = TRUE),
      width_pct    = round(100 * width_median / as.numeric(genome_len), 3),
      .groups = "drop"
    ) %>%
    arrange(desc(n_samples), desc(width_median))
  
  if (verbose) {
    cat("Genome length used for width_pct:", format(genome_len, big.mark = ","), "bp\n")
    if (nrow(sr_meta)) {
      cat("Top", min(20, nrow(sr_meta)), variant_name, "events by n_samples:\n")
      print(sr_meta %>% slice_head(n = 20), n = min(20, nrow(sr_meta)))
    } else {
      cat("No", variant_name, "events found in sv_rec.\n")
    }
  }
  
  # 6) optional saves
  if (isTRUE(save_files) && !is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(P,      file.path(out_dir, "P.rds"))
    saveRDS(P_sr,   file.path(out_dir, "P_sr.rds"))
    saveRDS(presence_long, file.path(out_dir, "presence_long.rds"))
    readr::write_csv(sr_meta, file.path(out_dir, "sr_meta.csv"))
  }
  
  # return clean bundle
  list(
    P = P,
    P_sr = P_sr,
    presence_long = presence_long,
    sr_meta = sr_meta,
    params = list(genome_len = genome_len, variant_name = variant_name)
  )
}

# --- rasr_basic_lists  [from ch4_functions.R:2017] ---
rasr_basic_lists <- function(sv_rec, res_pre, min_recur = 10, q_alpha = 0.05,
                             include_negative_exclusive = TRUE) {
  # 1) select the recurrent RASRs and get labels
  ev_tbl <- res_pre$sr_meta %>%
    filter(variant == "Structural_rearrangement", n_samples >= min_recur) %>%
    select(event_id, n_samples, width_median) %>%
    mutate(RASR_id = rasr_make_label(width_median))
  
  # 2) Build group-level presence matrix & restrict to these events
  Pg <- build_P_sr_group(sv_rec)
  keep <- intersect(ev_tbl$event_id, colnames(Pg))
  if (length(keep) < 2L) stop("Fewer than 2 recurrent RASRs found in P matrix.")
  X <- (Pg[, keep, drop = FALSE] != 0)
  k <- ncol(X); ids <- colnames(X)
  
  # 3) pairwise counts/stats
  n <- nrow(X)
  n1 <- Matrix::colSums(X)
  n11 <- Matrix::crossprod(X)  # k x k
  
  pairs <- tibble()
  for (i in 1:(k-1)) for (j in (i+1):k) {
    a <- as.numeric(n11[i, j])                  # n11
    b <- as.numeric(n1[i] - a)                  # n10
    cc <- as.numeric(n1[j] - a)                 # n01
    d <- as.numeric(n - a - b - cc)             # n00
    pairs <- bind_rows(pairs, tibble(
      event_A = ids[i], event_B = ids[j],
      n11 = a, n10 = b, n01 = cc, n00 = d,
      phi = .safe_phi(a,b,cc,d),
      p   = ._fisher_p(a,b,cc,d)
    ))
  }
  pairs <- pairs %>% mutate(q = p.adjust(p, method = "BH"))
  
  # 4) exclusivity & co-occurrence rules
  is_exclusive_pair <- function(a_row) {
    (a_row$n11 == 0) ||
      (include_negative_exclusive && (a_row$phi < 0 && a_row$q <= q_alpha))
  }
  is_cooc_pair <- function(a_row) (a_row$phi > 0 && a_row$q <= q_alpha && a_row$n11 > 0)
  
  # map to quick lookup lists
  excl_to <- lapply(ids, function(id) {
    withA <- pairs %>% filter(event_A == id)
    withB <- pairs %>% filter(event_B == id) %>% rename(event_A = event_B, event_B = event_A)
    pb <- bind_rows(withA, withB)
    pb %>% filter(purrr::pmap_lgl(., is_exclusive_pair)) %>% pull(event_B)
  })
  names(excl_to) <- ids
  
  cooc_with <- lapply(ids, function(id) {
    withA <- pairs %>% filter(event_A == id)
    withB <- pairs %>% filter(event_B == id) %>% rename(event_A = event_B, event_B = event_A)
    pb <- bind_rows(withA, withB)
    pb %>% filter(purrr::pmap_lgl(., is_cooc_pair)) %>% pull(event_B)
  })
  names(cooc_with) <- ids
  
  # 5) occurs_alone: any sample has this RASR and none of the others
  occ_alone <- sapply(seq_len(k), function(i) {
    any(X[, i] & rowSums(X[, -i, drop = FALSE]) == 0)
  })
  names(occ_alone) <- ids
  
  # 6) assemble a BASIC named list using the pretty RASR labels
  lab_map <- ev_tbl %>% filter(event_id %in% ids) %>% select(event_id, RASR_id)
  lab <- setNames(lab_map$RASR_id, lab_map$event_id)
  
  out <- lapply(ids, function(e) list(
    event_id = e,
    RASR_id = lab[[e]],
    occurs_alone = occ_alone[[e]],
    exclusive_to = unname(lab[excl_to[[e]]]),
    cooccurrent_with = unname(lab[cooc_with[[e]]])
  ))
  names(out) <- unname(lab[ids])
  out
}

# --- rasr_add_combos  [from ch4_functions.R:2290] ---
rasr_add_combos <- function(pres_wide, include_none = TRUE) {
  stopifnot("accession" %in% names(pres_wide))
  rasr_cols <- setdiff(names(pres_wide), "accession")
  M <- as.matrix(pres_wide[, rasr_cols, drop = FALSE])
  
  combo <- apply(M, 1, function(row) {
    idx <- which(row == 1L)
    if (length(idx) == 0) {
      if (include_none) "None" else NA_character_
    } else if (length(idx) == 1) {
      paste0(colnames(M)[idx], "_only")
    } else {
      paste(sort(colnames(M)[idx]), collapse = " + ")
    }
  })
  pres_wide %>% mutate(combo = combo)
}

# --- rasr_sample_by_combo  [from ch4_functions.R:2310] ---
rasr_sample_by_combo <- function(pres_with_combo, n_per = 25, seed = 1, drop_na_combo = TRUE) {
  set.seed(seed)
  df <- pres_with_combo
  if (drop_na_combo) df <- df %>% filter(!is.na(combo))
  sampled <- df %>%
    group_by(combo) %>%
    slice_sample(n = n_per) %>%
    ungroup()
  
  summary <- df %>%
    count(combo, name = "n_total") %>%
    left_join(sampled %>% count(combo, name = "n_sampled"), by = "combo") %>%
    mutate(n_sampled = coalesce(n_sampled, 0L)) %>%
    arrange(desc(n_total))
  
  list(sampled = sampled, summary = summary)
}

# --- .rasr_branch_labels  [from ch4_functions.R:2333] ---
.rasr_branch_labels <- function(sc_table, min_branch_exclusives = 2, force_branches = NULL) {
  split_csv <- function(x) {
    x <- ifelse(is.na(x) | nchar(trimws(x)) == 0, "", x)
    lapply(strsplit(x, "\\s*,\\s*"), function(v) unique(v[nchar(v) > 0]))
  }
  ex_sets <- split_csv(sc_table$exclusive_to); names(ex_sets) <- sc_table$RASR_id
  excl_counts <- vapply(ex_sets, length, integer(1))
  branch_labels <- names(excl_counts)[excl_counts >= min_branch_exclusives]
  occ_map <- setNames(sc_table$occurs_alone, sc_table$RASR_id)
  branch_labels <- branch_labels[isTRUE(occ_map[branch_labels])]
  if (!is.null(force_branches)) branch_labels <- unique(c(branch_labels, force_branches))
  branch_labels
}
