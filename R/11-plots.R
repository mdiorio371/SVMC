# 11-plots.R
# Auto-extracted canonical functions (SVMC consolidation).
# Provenance noted per function.

# --- plot_rasr_pathway_letters  [from ch4_functions.R:2347] ---
plot_rasr_pathway_letters <- function(
    P_sr_group, sc_table, sc_pairs,
    force_branches = NULL,
    
    node_size_base = 16, node_size_scale = 8.0,
    label_size = 8.0,   # A–F
    count_size = 6.8,   # sample counts
    edge_width_range = c(1.6, 4.8)
) {
  ev_ids <- sc_table$event_id
  id2lab <- setNames(sc_table$RASR_id, sc_table$event_id)
  lab2id <- setNames(sc_table$event_id, sc_table$RASR_id)
  
  X <- (P_sr_group[, ev_ids, drop = FALSE] != 0)
  ids <- colnames(X); k <- ncol(X)
  n11 <- if (k >= 2) as.matrix(Matrix::crossprod(X)) else matrix(0, k, k, dimnames = list(ids, ids))
  
  branch_labels <- .rasr_branch_labels(sc_table, force_branches = force_branches)
  branch_ids <- unname(lab2id[branch_labels]); branch_ids <- branch_ids[!is.na(branch_ids)]
  
  nodes <- tibble(
    node  = c("WT", id2lab[ev_ids]),
    core  = c("WT", rep("event", length(ev_ids))),
    label = c("WT", id2lab[ev_ids]),
    count = c(sum(Matrix::rowSums(X) == 0L), as.integer(Matrix::colSums(X)[ev_ids]))
  )
  # stable A–F mapping using sc$table row order
  rasr_labels <- id2lab[ev_ids]
  letters_map <- setNames(LETTERS[seq_along(rasr_labels)], rasr_labels)
  nodes$label_print <- ifelse(nodes$label == "WT", "WT", letters_map[nodes$label])
  nodes$size <- node_size_base + node_size_scale * log1p(pmax(1L, nodes$count))
  
  # Exclusive edges: WT->branches + parent→child for not-alone events (pick parent with max n11 among branches)
  edges_excl <- tibble()
  if (length(branch_ids)) {
    edges_excl <- tibble(
      from = "WT", to = id2lab[branch_ids], kind = "Exclusive",
      weight = as.integer(Matrix::colSums(X[, branch_ids, drop = FALSE])[branch_ids])
    )
  }
  child_ids <- sc_table %>% dplyr::filter(!occurs_alone) %>% pull(event_id)
  if (length(child_ids) && k >= 2) {
    for (ch in child_ids) {
      parent_pool <- intersect(branch_ids, ids); if (!length(parent_pool)) parent_pool <- ids
      a <- n11[match(ch, ids), match(parent_pool, ids)]; a[is.na(a)] <- 0L
      if (length(a) && max(a) > 0) {
        par_id <- parent_pool[which.max(a)]
        edges_excl <- dplyr::bind_rows(edges_excl,
                                       tibble(from = id2lab[par_id], to = id2lab[ch], kind = "Exclusive", weight = as.integer(max(a))))
      }
    }
  }
  
  # neutral edges (n11>0) excluding any already in Exclusive
  edges_neu <- tibble()
  if (k >= 2) {
    for (i in 1:(k-1)) for (j in (i+1):k) {
      a <- as.integer(n11[i, j]); if (a <= 0) next
      from_lab <- id2lab[ids[i]]; to_lab <- id2lab[ids[j]]
      already_excl <- nrow(dplyr::filter(edges_excl,
                                         (from == from_lab & to == to_lab) | (from == to_lab & to == from_lab))) > 0
      if (!already_excl) edges_neu <- dplyr::bind_rows(edges_neu,
                                                       tibble(from = from_lab, to = to_lab, kind = "neutral", weight = a))
    }
  }
  edges <- dplyr::bind_rows(edges_excl, edges_neu)
  
  pal <- c("#999999", scales::hue_pal()(length(rasr_labels)))
  base_cols <- setNames(pal[-1], rasr_labels)
  node_cols <- base_cols[nodes$label]; node_cols[is.na(node_cols)] <- "#999999"
  
  g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  
  ggraph::ggraph(g, layout = "sugiyama") +
    # edges first
    ggraph::geom_edge_link(aes(width = weight, linetype = kind), color = "black", alpha = 0.9) +
    ggraph::scale_edge_width(range = edge_width_range, guide = "none") +
    ggraph::scale_edge_linetype_manual(values = c(Exclusive = "solid", neutral = "dotted"),
                                       name = "Relationship") +
    # big nodes
    ggraph::geom_node_point(aes(size = size), shape = 21, fill = node_cols, color = "grey30", stroke = 0.9) +
    geom_node_point(aes(size = size), shape = 21, fill = node_cols, color = "grey30", stroke = 0.9) +
    scale_size_continuous(range = c(12, 30), guide = "none")  +
    # label (A–F) on a white tag above node to avoid any edge conflicts
    ggraph::geom_node_label(aes(label = label_print),
                            vjust = -1.15, size = label_size, fontface = "bold",
                            label.size = 0, label.padding = unit(c(0.1, 0.2, 0.1, 0.2), "lines"),
                            fill = "white") +
    # count below node, also on white
    # ggraph::geom_node_label(aes(label = count),
    #                         vjust = 1.9, size = count_size,
    #                         label.size = 0, label.padding = unit(c(0.05, 0.15, 0.05, 0.15), "lines"),
    #                         fill = "white") +
    ggplot2::guides(size = "none") +
    ggplot2::theme_void(base_size = 13) +
    coord_cartesian(clip = "off") +                             # allow labels to extend past panel
    theme(plot.margin = margin(t = 28, r = 16, b = 16, l = 16))
}

# --- synteny_panel  [from ch4_functions.R:2450] ---
synteny_panel <- function(colinear = list(), inversions = list(),
                          title = "", xlab = "WT", ylab = "") {
  seg_df <- function(lst, col) {
    if (length(lst) == 0) return(NULL)
    do.call(rbind, lapply(lst, function(v) {
      data.frame(x = v[1], xend = v[2], y = v[3], yend = v[4], col = col)
    }))
  }
  d_all <- rbind(seg_df(colinear, "blue"), seg_df(inversions, "red"))
  ggplot(d_all) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend, color = col),
                 linewidth = 1.9, lineend = "round") +
    scale_color_identity() +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

# --- panel1b_synteny_grid  [from ch4_functions.R:2491] ---
panel1b_synteny_grid <- function(sc_table, reflen) {
  # order as they appear in sc$table
  rasr_labels <- sc_table$RASR_id
  # A–F letters
  letters_map <- setNames(LETTERS[seq_along(rasr_labels)], rasr_labels)
  
  # parse Mbp from label (e.g., "RASR_0.86Mbp")
  len_mbp <- as.numeric(sub("Mbp$", "", sub("^RASR_", "", rasr_labels)))
  len_bp  <- len_mbp * 1e6
  
  plots <- vector("list", length(rasr_labels))
  for (i in seq_along(rasr_labels)) {
    segs <- synteny_from_length(reflen, len_bp[i])
    title <- paste0(format(len_mbp[i], nsmall = 2), " Mbp (", letters_map[rasr_labels[i]], ")")
    plots[[i]] <- synteny_panel(segs$colinear, segs$inversions,
                                title = title, xlab = "WT", ylab = letters_map[rasr_labels[i]])
  }
  
  # 2x3 grid
  gridExtra::grid.arrange(grobs = plots, ncol = 3)
}

# --- rasr_synteny_segments  [from ch4_functions.R:2550] ---
rasr_synteny_segments <- function(reflen, inv_bp, center = 0.5) {
  stopifnot(is.finite(reflen), reflen > 0)
  inv_bp <- max(0, inv_bp)
  f <- min(inv_bp / reflen, 0.98)   # keep some flanks visible
  if (f == 0) {
    return(data.frame(x=0, y=0, xend=1, yend=1, col="blue", seg="colinear"))
  }
  a <- min(max(center - f/2, 0), 1 - f); b <- a + f
  out <- list()
  if (a > 1e-6) out[[length(out)+1]] <- data.frame(x=0, y=0, xend=a, yend=a, col="blue", seg="colinear_left")
  out[[length(out)+1]] <- data.frame(x=a, y=b, xend=b, yend=a, col="red", seg="inversion")
  if ((1-b) > 1e-6) out[[length(out)+1]] <- data.frame(x=b, y=b, xend=1, yend=1, col="blue", seg="colinear_right")
  do.call(rbind, out)
}

# --- synteny_panel_classic  [from ch4_functions.R:2659] ---
synteny_panel_classic <- function(reflen, inv_bp, length_mbp, n,
                                  xlab = "WT", ylab = "",
                                  colinear_col = "#999999",   # WT color
                                  inversion_col = "#E64B35")  # event color
{
  d_all <- rasr_synteny_segments(reflen, inv_bp)
  ann   <- paste0("Length = ", format(length_mbp, nsmall = 2), " Mbp\nn = ", n)
  
  ggplot() +
    geom_segment(data = subset(d_all, seg != "inversion"),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = 2.4, lineend = "round", color = colinear_col) +
    geom_segment(data = subset(d_all, seg == "inversion"),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = 2.4, lineend = "round", color = inversion_col) +
    coord_fixed(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
    annotate("text", x = 0.3, y = 0.98, label = ann, hjust = 0, vjust = 1,
             fontface = "bold", size = 5.2) +
    labs(title = NULL, x = xlab, y = ylab) +
    theme_classic(base_size = 14) +
    theme(panel.grid = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 18, face = "bold"))
}

# --- rasr_letter_colors  [from ch4_functions.R:2620] ---
rasr_letter_colors <- function(sc_table, wt_col = "#999999", palette = NULL) {
  rasr_labels <- sc_table$RASR_id
  if (is.null(palette)) palette <- scales::hue_pal()(length(rasr_labels))
  base_cols <- setNames(palette, rasr_labels)             # event color by RASR_id
  letters_map <- setNames(LETTERS[seq_along(rasr_labels)], rasr_labels)
  list(
    wt      = wt_col,
    by_id   = base_cols,                                  # RASR_id -> color
    by_letter = setNames(unname(base_cols), letters_map[rasr_labels])  # A/B/.. -> color
  )
}

# --- svcm_panelA_prop  [from SVCM6.R:1348] ---
svcm_panelA_prop <- function(sv_ro_tbl) {
  sv_ro_tbl %>%
    filter(!is.na(ro_bin)) %>%
    count(species, ro_bin) %>%
    group_by(species) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
}

# --- synteny_from_length  [from ch4_functions.R:2474] ---
synteny_from_length <- function(reflen, inv_bp, y_top = 0.84, y_mid = 0.5, y_bot = 0.16) {
  frac <- pmin(pmax(inv_bp / reflen, 0), 1)
  half <- frac / 2
  x0 <- 0.5 - half; x1 <- 0.5 + half
  x0 <- pmax(0, x0); x1 <- pmin(1, x1)
  
  colinear <- list(
    c(0,   x0,  y_mid, y_mid),  # left flank
    c(x1,  1,   y_mid, y_mid)   # right flank
  )
  inversions <- list(
    c(x0, x1, y_top, y_bot)     # flipped block (diagonal red)
  )
  list(colinear = colinear, inversions = inversions)
}
