#############################################################
# Simulation study:
# APA, Pearson correlation, and ARI under hierarchical and
# cross-cutting pile-sort organisations
############################################################
#
# The script compares three respondent-by-respondent agreement
# measures for unconstrained pile sorts:
#
#   1. APA: adjusted positive agreement
#   2. Pearson correlation on binary co-placement vectors
#   3. Adjusted Rand Index
#
# Two conditions are simulated:
#
#   single_culture_scenario:
#     All respondents sort according to the same taxonomic hierarchy H.
#     Respondents differ only in classificatory resolution.
#
#   multiple_cultures_scenario:
#     75% of respondents sort according to H.
#     25% sort according to a cross-cutting ecological organisation H*.
#
# The main diagnostics are:
#
#   - mean respondent-by-respondent agreement
#   - eigenvalue spectra of the respondent-by-respondent matrices
#   - eta-squared recovery of known grouping variables from eigenvectors
#
########################################################


##########################################################
# 0. Setup
#######################################################

set.seed(123)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
  library(ggraph)
  library(grid)
  library(patchwork)
})

N_SIMS <- 1000
R <- 40

OUT_DIR <- "."
SAVE_FIGURES <- TRUE
SAVE_TABLES <- TRUE


# Shared visual settings for eigenanalysis figures
eigen_base_size <- 14
eigen_strip_size <- 14
eigen_axis_title_size <- 15
eigen_axis_text_size <- 13


# Visual settings for the single-culture combined figure
single_base_size <- 11
single_strip_size <- 10
single_axis_title_size <- 12
single_axis_text_size <- 10
single_tile_text_size <- 3.6


save_plot <- function(filename, plot, width, height, dpi = 300) {
  if (!SAVE_FIGURES) return(invisible(NULL))
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}


############################################################
# 1. Latent item organisations
############################################################
#
# H is the taxonomic hierarchy. It defines a nested organisation of
# animal items. Respondents in the taxonomic condition are simulated
# as cutting this same tree at different levels of resolution.
############################################################

animal_hierarchy <- list(
  Mammals = list(
    Carnivores = list(
      Canines = c("dog", "wolf", "fox"),
      Felines = c("cat", "lion", "tiger")
    ),
    Hoofed_mammals = list(
      Equines = c("horse", "zebra", "donkey"),
      Ruminants = c("cow", "sheep", "goat")
    ),
    Primates = c("human", "chimpanzee", "gorilla")
  ),
  Birds = list(
    Birds_of_prey = c("eagle", "hawk", "owl"),
    Water_birds = c("duck", "goose", "swan"),
    Small_birds = c("sparrow", "robin", "pigeon")
  ),
  Reptiles_Amphibians = list(
    Reptiles = list(
      Scaled_reptiles = c("snake", "lizard", "gecko"),
      Armoured_reptiles = c("turtle", "crocodile", "alligator")
    ),
    Amphibians = c("frog", "toad", "salamander")
  ),
  Fish_Marine = list(
    Fish = list(
      Common_fish = c("salmon", "cod", "trout"),
      Large_fish = c("shark", "tuna", "eel")
    ),
    Marine_mammals = c("whale", "dolphin", "seal")
  ),
  Invertebrates = list(
    Insects = list(
      Social_insects = c("ant", "bee", "wasp"),
      Flying_insects = c("butterfly", "fly", "mosquito")
    )
  )
)


# Tree utilities ------------------------------------------------------
#
# The hierarchy is stored as nested lists. These helpers move between
# item names and integer item ids which are easier to use when
# constructing matrices.

get_leaves <- function(x) {
  if (is.character(x)) return(x)
  unlist(lapply(x, get_leaves), use.names = FALSE)
}

animals <- get_leaves(animal_hierarchy)
n_items <- length(animals)

items_to_ids <- function(items) {
  match(items, animals)
}

node_items <- function(node) {
  items_to_ids(get_leaves(node))
}

cat("\nNumber of animal items:", n_items, "\n")


# H*: cross-cutting ecological organisation --------------------------
#
# H* is not a hierarchy. It assigns the same animals to size and
# habitat dimensions. This gives a second classificatory logic that
# cuts across the taxonomic tree.

animal_traits <- data.frame(
  animal = animals,
  stringsAsFactors = FALSE
) %>%
  mutate(
    habitat = case_when(
      animal %in% c(
        "eagle", "hawk", "owl",
        "sparrow", "robin", "pigeon",
        "bee", "wasp", "butterfly", "fly", "mosquito"
      ) ~ "air",

      animal %in% c(
        "duck", "goose", "swan",
        "turtle", "crocodile", "alligator",
        "frog", "toad", "salamander",
        "salmon", "cod", "trout", "shark", "tuna", "eel",
        "whale", "dolphin", "seal"
      ) ~ "water",

      TRUE ~ "land"
    ),

    size = case_when(
      animal %in% c(
        "dog", "wolf", "lion", "tiger",
        "horse", "zebra", "donkey",
        "cow", "sheep", "goat",
        "human", "chimpanzee", "gorilla",
        "eagle", "hawk", "owl",
        "duck", "goose", "swan",
        "crocodile", "alligator",
        "shark", "tuna", "whale", "dolphin", "seal"
      ) ~ "big",

      TRUE ~ "small"
    )
  )

stopifnot(setequal(animal_traits$animal, animals))


############################################################
# 2. Visualise H and H*
############################################################

make_hstar_plot <- function(traits) {
  set.seed(123)

  plot_data <- traits %>%
    mutate(
      habitat = factor(habitat, levels = c("land", "water", "air")),
      size = factor(size, levels = c("small", "big")),
      x = as.numeric(habitat),
      y = as.numeric(size),
      x_jit = x + runif(n(), -0.15, 0.15),
      y_jit = y + runif(n(), -0.08, 0.08)
    )

  ggplot(plot_data) +
    geom_point(
      aes(x = x_jit, y = y_jit),
      size = 2.2,
      alpha = 0.85
    ) +
    geom_text_repel(
      aes(x = x_jit, y = y_jit, label = animal),
      size = 4.5,
      segment.color = NA,
      box.padding = 0.15,
      point.padding = 0.08,
      max.overlaps = Inf
    ) +
    scale_x_continuous(
      breaks = 1:3,
      labels = c("land", "water", "air"),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_y_continuous(
      breaks = 1:2,
      labels = c("small", "big"),
      expand = expansion(mult = c(0.08, 0.08))
    ) +
    coord_cartesian(
      xlim = c(0.7, 3.3),
      ylim = c(0.75, 2.25)
    ) +
    labs(
      x = "Habitat",
      y = "Size"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 17),
      plot.margin = margin(5, 5, 5, 5)
    )
}


hierarchy_to_edges <- function(x, parent = "Animals") {
  if (is.character(x)) {
    return(data.frame(from = parent, to = x))
  }

  edges <- data.frame(from = character(), to = character())

  for (nm in names(x)) {
    # Invertebrates are flattened by one level so that the plot is
    # not visually dominated by a single intermediate category.
    if (nm == "Invertebrates") {
      child <- x[[nm]]

      for (child_nm in names(child)) {
        edges <- rbind(
          edges,
          data.frame(from = parent, to = child_nm),
          hierarchy_to_edges(child[[child_nm]], parent = child_nm)
        )
      }
    } else {
      edges <- rbind(
        edges,
        data.frame(from = parent, to = nm),
        hierarchy_to_edges(x[[nm]], parent = nm)
      )
    }
  }

  edges
}


make_taxonomy_plot <- function(tree) {
  edges <- hierarchy_to_edges(tree)
  graph <- graph_from_data_frame(edges)

  V(graph)$label <- gsub("_", " ", V(graph)$name)

  ggraph(graph, layout = "dendrogram", circular = TRUE) +
    geom_edge_diagonal(
      alpha = 0.65,
      linewidth = 0.35
    ) +
    geom_node_point(
      size = 1.8,
      alpha = 0.85
    ) +
    geom_node_label(
      aes(label = label),
      size = 3.8,
      label.size = 0,
      label.padding = unit(0.10, "lines"),
      fill = "white",
      alpha = 0.95
    ) +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(2, 2, 2, 2)
    )
}


p_Hstar_grid <- make_hstar_plot(animal_traits)
p_circle <- make_taxonomy_plot(animal_hierarchy)

print(p_Hstar_grid)
print(p_circle)

save_plot("Hstar_size_habitat_grid.png", p_Hstar_grid, width = 8, height = 5, dpi = 300)
save_plot("H_taxonomic_hierarchy_circular.png", p_circle, width = 10, height = 10, dpi = 600)


p_H <- p_circle +
  labs(title = "") +
  coord_fixed(clip = "off") +
  theme(
    plot.margin = margin(0, -45, 0, -30)
  )

p_Hstar <- p_Hstar_grid +
  labs(title = "") +
  theme(
    plot.margin = margin(5, 0, 5, -20)
  )

p_latent_organisations <- p_H + p_Hstar +
  plot_layout(widths = c(0.95, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "",
    tag_suffix = "",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 22)
    )
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 22),
    plot.tag.position = c(0.07, 0.98)
  )

print(p_latent_organisations)

save_plot(
  "Figure_2_latent_organisations_H_and_Hstar.png",
  p_latent_organisations,
  width = 13,
  height = 7,
  dpi = 600
)


############################################################
# 3. Generate respondent partitions
#
# Respondents differ in resolution:
#   low_resolution          5-8 piles
#   intermediate_resolution 9-11 piles
#   high_resolution         12-16 piles
#
# The probabilities create a bimodal lumper-splitter distribution,
# with fewer respondents in the intermediate category.
############################################################

draw_resolution <- function() {
  sample(
    x = c("low_resolution", "intermediate_resolution", "high_resolution"),
    size = 1,
    prob = c(0.40, 0.20, 0.40)
  )
}


# Partitions from the taxonomic hierarchy ----------------------------

get_broad_partition <- function(tree) {
  lapply(tree, node_items)
}

get_intermediate_partition <- function(tree) {
  out <- list()

  for (broad in names(tree)) {
    for (mid in names(tree[[broad]])) {
      out[[paste(broad, mid, sep = "/")]] <- node_items(tree[[broad]][[mid]])
    }
  }

  out
}

get_fine_partition <- function(tree) {
  out <- list()

  for (broad in names(tree)) {
    for (mid in names(tree[[broad]])) {
      if (is.character(tree[[broad]][[mid]])) {
        out[[paste(broad, mid, sep = "/")]] <- node_items(tree[[broad]][[mid]])
      } else {
        for (fine in names(tree[[broad]][[mid]])) {
          out[[paste(broad, mid, fine, sep = "/")]] <-
            node_items(tree[[broad]][[mid]][[fine]])
        }
      }
    }
  }

  out
}


make_taxonomic_partition <- function(resolution, tree = animal_hierarchy) {
  if (resolution == "low_resolution") {

    # Start with broad groups and occasionally refine them.
    # This produces coarse cuts through the hierarchy.
    part <- get_broad_partition(tree)
    n_refine <- sample(0:2, 1)

    if (n_refine > 0) {
      refine_names <- sample(names(tree), n_refine)

      for (broad in refine_names) {
        part[[broad]] <- NULL

        for (mid in names(tree[[broad]])) {
          part[[paste(broad, mid, sep = "/")]] <-
            node_items(tree[[broad]][[mid]])
        }
      }
    }

    if (length(part) > 8) {
      return(make_taxonomic_partition("low_resolution", tree))
    }

  } else if (resolution == "intermediate_resolution") {

    # Intermediate respondents start from the mid-level taxonomy and
    # may collapse one or two broad groups upward.
    part <- get_intermediate_partition(tree)

    intermediate_move <- sample(
      x = c("collapse_two", "collapse_one", "no_change"),
      size = 1
    )

    if (intermediate_move %in% c("collapse_one", "collapse_two")) {
      n_collapse <- ifelse(intermediate_move == "collapse_one", 1, 2)
      collapse_broads <- sample(names(tree), size = n_collapse)

      for (broad in collapse_broads) {
        for (mid in names(tree[[broad]])) {
          mid_name <- paste(broad, mid, sep = "/")
          part[[mid_name]] <- NULL

          if (!is.character(tree[[broad]][[mid]])) {
            for (fine in names(tree[[broad]][[mid]])) {
              fine_name <- paste(broad, mid, fine, sep = "/")
              part[[fine_name]] <- NULL
            }
          }
        }

        part[[broad]] <- node_items(tree[[broad]])
      }
    }

    if (length(part) < 9 || length(part) > 11) {
      return(make_taxonomic_partition("intermediate_resolution", tree))
    }

  } else if (resolution == "high_resolution") {

    # Start from the fine taxonomy and occasionally collapse fine
    # branches. This keeps high-resolution respondents detailed but
    # not identical to one another.
    part <- get_fine_partition(tree)

    collapsible_mid <- data.frame()

    for (broad in names(tree)) {
      for (mid in names(tree[[broad]])) {
        if (!is.character(tree[[broad]][[mid]])) {
          collapsible_mid <- rbind(
            collapsible_mid,
            data.frame(
              broad = broad,
              mid = mid,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }

    n_collapse <- sample(0:4, 1)

    if (n_collapse > 0) {
      collapse_rows <- collapsible_mid[
        sample(seq_len(nrow(collapsible_mid)), n_collapse),
      ]

      for (i in seq_len(nrow(collapse_rows))) {
        broad <- collapse_rows$broad[i]
        mid <- collapse_rows$mid[i]

        for (fine in names(tree[[broad]][[mid]])) {
          fine_name <- paste(broad, mid, fine, sep = "/")
          part[[fine_name]] <- NULL
        }

        part[[paste(broad, mid, sep = "/")]] <-
          node_items(tree[[broad]][[mid]])
      }
    }

    if (length(part) < 12 || length(part) > 16) {
      return(make_taxonomic_partition("high_resolution", tree))
    }

  } else {
    stop("Unknown resolution: ", resolution)
  }

  part
}


# Partitions from the ecological organisation ------------------------
#
# H* respondents sort by habitat, size, or the size-by-habitat
# combination. The resolution label controls the granularity of the
# cross-cutting partition.

make_ecological_partition <- function(resolution, traits = animal_traits) {
  if (resolution == "low_resolution") {

    principle <- sample(
      x = c("habitat", "size"),
      size = 1,
      prob = c(0.70, 0.30)
    )
    labels <- traits[[principle]]

  } else if (resolution == "intermediate_resolution") {

    # Mostly habitat, but with size distinguished within one habitat.
    target_habitat <- sample(unique(traits$habitat), 1)
    labels <- traits$habitat

    labels[traits$habitat == target_habitat] <- paste(
      traits$habitat[traits$habitat == target_habitat],
      traits$size[traits$habitat == target_habitat],
      sep = "_"
    )

  } else if (resolution == "high_resolution") {

    labels <- paste(traits$habitat, traits$size, sep = "_")

  } else {
    stop("Unknown resolution: ", resolution)
  }

  split(seq_len(nrow(traits)), labels)
}


############################################################
# 4. Agreement measures
############################################################
#
# All three measures compare respondents pairwise, but they do so from
# different representations of the same sort.
#
#   APA and Pearson use binary co-placement vectors.
#   ARI uses item labels.
#
# Co-placement vectors contain one entry for each item pair:
#   1 = the respondent placed the pair in the same pile
#   0 = otherwise.
############################################################

part_to_A <- function(part, n) {
  mat <- matrix(FALSE, n, n)

  for (pile in part) {
    mat[pile, pile] <- TRUE
  }

  as.integer(mat[upper.tri(mat, diag = FALSE)])
}

part_to_labels <- function(part, n) {
  labels <- rep(NA_integer_, n)

  for (p in seq_along(part)) {
    labels[part[[p]]] <- p
  }

  labels
}

sorting_density <- function(part) {
  sum(choose(vapply(part, length, integer(1)), 2))
}

upper_vals <- function(M) {
  M[upper.tri(M)]
}


APA <- function(A_i, A_j) {
  observed <- sum(A_i * A_j)
  s_i <- sum(A_i)
  s_j <- sum(A_j)
  n_pairs <- length(A_i)

  expected <- (s_i * s_j) / n_pairs
  maximum <- min(s_i, s_j)
  denominator <- maximum - expected

  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }

  value <- (observed - expected) / denominator

  max(0, min(1, value))
}

Pearson <- function(A_i, A_j) {
  suppressWarnings(cor(A_i, A_j))
}

Adjusted_Rand <- function(labels_i, labels_j) {
  tab <- table(labels_i, labels_j)

  nij <- as.vector(tab)
  ai <- rowSums(tab)
  bj <- colSums(tab)
  n <- sum(tab)

  sum_nij <- sum(choose(nij, 2))
  sum_ai <- sum(choose(ai, 2))
  sum_bj <- sum(choose(bj, 2))
  total_pairs <- choose(n, 2)

  expected <- (sum_ai * sum_bj) / total_pairs
  maximum <- 0.5 * (sum_ai + sum_bj)
  denominator <- maximum - expected

  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }

  (sum_nij - expected) / denominator
}


############################################################
# 5. One simulation replicate
#
# A replicate returns:
#   - respondent metadata
#   - pairwise respondent results
#   - dataset-level summary values
#   - the three respondent-by-respondent matrices
############################################################

run_one_simulation <- function(sim_id,
                               condition,
                               p_ecological,
                               R = 40) {
  respondents <- vector("list", R)

  meta <- data.frame(
    sim = integer(),
    condition = character(),
    respondent = integer(),
    resolution = character(),
    latent_organisation = character(),
    k = integer(),
    density = numeric(),
    stringsAsFactors = FALSE
  )

  for (r in seq_len(R)) {
    resolution <- draw_resolution()

    latent_organisation <- if (runif(1) < p_ecological) {
      "ecological_Hstar"
    } else {
      "taxonomic_H"
    }

    part_r <- if (latent_organisation == "taxonomic_H") {
      make_taxonomic_partition(resolution)
    } else {
      make_ecological_partition(resolution)
    }

    respondents[[r]] <- part_r

    meta <- rbind(
      meta,
      data.frame(
        sim = sim_id,
        condition = condition,
        respondent = r,
        resolution = resolution,
        latent_organisation = latent_organisation,
        k = length(part_r),
        density = sorting_density(part_r),
        stringsAsFactors = FALSE
      )
    )
  }

  A_list <- lapply(respondents, part_to_A, n = n_items)
  label_list <- lapply(respondents, part_to_labels, n = n_items)

  APA_mat <- matrix(NA_real_, R, R)
  Pearson_mat <- matrix(NA_real_, R, R)
  ARI_mat <- matrix(NA_real_, R, R)

  for (i in seq_len(R)) {
    for (j in seq_len(R)) {
      APA_mat[i, j] <- if (i == j) 1 else APA(A_list[[i]], A_list[[j]])
      Pearson_mat[i, j] <- if (i == j) 1 else Pearson(A_list[[i]], A_list[[j]])
      ARI_mat[i, j] <- if (i == j) 1 else Adjusted_Rand(label_list[[i]], label_list[[j]])
    }
  }

  pair_df <- expand.grid(i = seq_len(R), j = seq_len(R)) %>%
    filter(i < j) %>%
    left_join(
      meta %>%
        select(
          i = respondent,
          resolution_i = resolution,
          latent_i = latent_organisation
        ),
      by = "i"
    ) %>%
    left_join(
      meta %>%
        select(
          j = respondent,
          resolution_j = resolution,
          latent_j = latent_organisation
        ),
      by = "j"
    ) %>%
    mutate(
      sim = sim_id,
      condition = condition,
      resolution_pair = paste(resolution_i, resolution_j, sep = "-"),
      latent_pair = paste(latent_i, latent_j, sep = "-"),
      APA = mapply(function(i, j) APA_mat[i, j], i, j),
      Pearson = mapply(function(i, j) Pearson_mat[i, j], i, j),
      ARI = mapply(function(i, j) ARI_mat[i, j], i, j)
    )

  summary_df <- data.frame(
    sim = sim_id,
    condition = condition,
    p_ecological = p_ecological,

    APA_mean = mean(upper_vals(APA_mat), na.rm = TRUE),
    Pearson_mean = mean(upper_vals(Pearson_mat), na.rm = TRUE),
    ARI_mean = mean(upper_vals(ARI_mat), na.rm = TRUE),

    APA_sd = sd(upper_vals(APA_mat), na.rm = TRUE),
    Pearson_sd = sd(upper_vals(Pearson_mat), na.rm = TRUE),
    ARI_sd = sd(upper_vals(ARI_mat), na.rm = TRUE)
  )

  list(
    meta = meta,
    pairs = pair_df,
    summary = summary_df,
    APA_mat = APA_mat,
    Pearson_mat = Pearson_mat,
    ARI_mat = ARI_mat
  )
}


############################################################
# 6. Run all simulations
############################################################

condition_grid <- data.frame(
  condition = c(
    "single_culture_scenario",
    "multiple_cultures_scenario"
  ),
  p_ecological = c(0.00, 0.25),
  stringsAsFactors = FALSE
)

sim_grid <- expand.grid(
  sim = seq_len(N_SIMS),
  row_id = seq_len(nrow(condition_grid))
) %>%
  left_join(
    condition_grid %>% mutate(row_id = row_number()),
    by = "row_id"
  )

sim_list <- vector("list", nrow(sim_grid))

for (s in seq_len(nrow(sim_grid))) {
  sim_list[[s]] <- run_one_simulation(
    sim_id = sim_grid$sim[s],
    condition = sim_grid$condition[s],
    p_ecological = sim_grid$p_ecological[s],
    R = R
  )
}

meta_all <- bind_rows(lapply(sim_list, `[[`, "meta"))
pair_all <- bind_rows(lapply(sim_list, `[[`, "pairs"))
summary_all <- bind_rows(lapply(sim_list, `[[`, "summary"))


summary_print <- summary_all %>%
  group_by(condition, p_ecological) %>%
  summarise(
    APA_mean = mean(APA_mean),
    APA_sd_across_sims = sd(APA_mean),

    Pearson_mean = mean(Pearson_mean),
    Pearson_sd_across_sims = sd(Pearson_mean),

    ARI_mean = mean(ARI_mean),
    ARI_sd_across_sims = sd(ARI_mean),

    .groups = "drop"
  )

print(summary_print)


############################################################
# 7. Single-culture agreement heatmap
############################################################
#
# This figure checks whether measures preserve agreement across
# respondents who share the same hierarchy but cut it at different
# levels of resolution.
############################################################

pair_heat <- pair_all %>%
  filter(condition == "single_culture_scenario") %>%
  mutate(
    resolution_i = factor(
      resolution_i,
      levels = c("low_resolution", "intermediate_resolution", "high_resolution"),
      labels = c("Lumper", "Intermediate", "Splitter")
    ),
    resolution_j = factor(
      resolution_j,
      levels = c("low_resolution", "intermediate_resolution", "high_resolution"),
      labels = c("Lumper", "Intermediate", "Splitter")
    )
  ) %>%
  group_by(resolution_i, resolution_j) %>%
  summarise(
    APA = mean(APA, na.rm = TRUE),
    Pearson = mean(Pearson, na.rm = TRUE),
    ARI = mean(ARI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(APA, Pearson, ARI),
    names_to = "measure",
    values_to = "mean_agreement"
  ) %>%
  mutate(
    measure = recode(
      measure,
      APA = "APA",
      Pearson = "Pearson correlation",
      ARI = "Adjusted Rand Index"
    ),
    measure = factor(
      measure,
      levels = c("APA", "Pearson correlation", "Adjusted Rand Index")
    )
  )

p_heat <- ggplot(
  pair_heat,
  aes(x = resolution_i, y = resolution_j, fill = mean_agreement)
) +
  geom_tile(
    color = "white",
    linewidth = 1
  ) +
  geom_text(
    aes(label = round(mean_agreement, 2)),
    size = 5
  ) +
  facet_wrap(~ measure, nrow = 1) +
  coord_equal() +
  scale_fill_gradient(
    limits = c(0, 1),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    panel.grid = element_blank(),
    legend.position = "none"
  )

print(p_heat)
save_plot("sim_heatmap_single_culture_APA_Pearson_ARI.png", p_heat, width = 11, height = 4.2)


############################################################
# 8. Structure-pair diagnostic
############################################################
#
# A compact check of mean agreement by latent pairing. This is most
# useful for inspecting the multiple-culture condition.
############################################################

structure_summary <- pair_all %>%
  group_by(condition, latent_pair) %>%
  summarise(
    n_pairs = n(),
    APA_mean = mean(APA, na.rm = TRUE),
    Pearson_mean = mean(Pearson, na.rm = TRUE),
    ARI_mean = mean(ARI, na.rm = TRUE),
    .groups = "drop"
  )

print(structure_summary, n = Inf)


############################################################
# 9. Eigenanalysis of respondent-by-respondent matrices
############################################################
#
# The eigenvectors are used descriptively. Respondent coordinates on
# the leading eigenvectors are compared with known simulated grouping
# variables: resolution and latent organisation.
############################################################

analyse_eigen <- function(M, meta, measure, condition, sim_id) {
  M <- (M + t(M)) / 2
  eig <- eigen(M, symmetric = TRUE)

  values <- eig$values
  vectors <- eig$vectors

  eigenvalues <- data.frame(
    sim = sim_id,
    condition = condition,
    measure = measure,
    component = seq_along(values),
    eigenvalue = values,
    abs_eigenvalue = abs(values),
    prop_abs_eigenvalue = abs(values) / sum(abs(values)),
    negative = values < 0
  )

  eigenvectors <- data.frame(
    sim = sim_id,
    condition = condition,
    measure = measure,
    respondent = meta$respondent,
    resolution = meta$resolution,
    latent_organisation = meta$latent_organisation,
    EV1 = vectors[, 1],
    EV2 = vectors[, 2],
    EV3 = vectors[, 3],
    EV4 = vectors[, 4]
  )

  list(
    eigenvalues = eigenvalues,
    eigenvectors = eigenvectors
  )
}


eigen_list <- list()
counter <- 1

for (idx in seq_along(sim_list)) {
  sim_obj <- sim_list[[idx]]
  meta_i <- sim_obj$meta

  sim_id_i <- unique(meta_i$sim)
  condition_i <- unique(meta_i$condition)

  eigen_list[[counter]] <- analyse_eigen(
    M = sim_obj$APA_mat,
    meta = meta_i,
    measure = "APA",
    condition = condition_i,
    sim_id = sim_id_i
  )
  counter <- counter + 1

  eigen_list[[counter]] <- analyse_eigen(
    M = sim_obj$Pearson_mat,
    meta = meta_i,
    measure = "Pearson correlation",
    condition = condition_i,
    sim_id = sim_id_i
  )
  counter <- counter + 1

  eigen_list[[counter]] <- analyse_eigen(
    M = sim_obj$ARI_mat,
    meta = meta_i,
    measure = "Adjusted Rand Index",
    condition = condition_i,
    sim_id = sim_id_i
  )
  counter <- counter + 1
}

eigenvalues_all <- bind_rows(lapply(eigen_list, `[[`, "eigenvalues"))
eigenvectors_all <- bind_rows(lapply(eigen_list, `[[`, "eigenvectors"))


# Eta-squared recovery statistic -------------------------------------
#
# eta^2 = SS_between / SS_total
#
# Here x is an eigenvector coordinate and group is a known simulated
# grouping variable. High eta^2 means that the eigenvector separates
# respondents according to that variable.

eta_squared <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]
  group <- as.factor(group[ok])

  if (length(unique(group)) < 2) {
    return(NA_real_)
  }

  grand_mean <- mean(x)
  group_means <- tapply(x, group, mean)
  group_ns <- table(group)

  ss_between <- sum(group_ns * (group_means - grand_mean)^2)
  ss_total <- sum((x - grand_mean)^2)

  if (ss_total <= 0) {
    return(NA_real_)
  }

  ss_between / ss_total
}


ev_separation <- eigenvectors_all %>%
  group_by(condition, measure, sim) %>%
  summarise(
    EV1_eta_resolution = eta_squared(EV1, resolution),
    EV2_eta_resolution = eta_squared(EV2, resolution),
    EV3_eta_resolution = eta_squared(EV3, resolution),
    EV4_eta_resolution = eta_squared(EV4, resolution),

    EV1_eta_latent = eta_squared(EV1, latent_organisation),
    EV2_eta_latent = eta_squared(EV2, latent_organisation),
    EV3_eta_latent = eta_squared(EV3, latent_organisation),
    EV4_eta_latent = eta_squared(EV4, latent_organisation),

    .groups = "drop"
  )

ev_separation_summary <- ev_separation %>%
  group_by(condition, measure) %>%
  summarise(
    EV1_eta_resolution = mean(EV1_eta_resolution, na.rm = TRUE),
    EV2_eta_resolution = mean(EV2_eta_resolution, na.rm = TRUE),
    EV3_eta_resolution = mean(EV3_eta_resolution, na.rm = TRUE),
    EV4_eta_resolution = mean(EV4_eta_resolution, na.rm = TRUE),

    EV1_eta_latent = mean(EV1_eta_latent, na.rm = TRUE),
    EV2_eta_latent = mean(EV2_eta_latent, na.rm = TRUE),
    EV3_eta_latent = mean(EV3_eta_latent, na.rm = TRUE),
    EV4_eta_latent = mean(EV4_eta_latent, na.rm = TRUE),

    .groups = "drop"
  )

print(ev_separation_summary, n = Inf)


# Eigenvalue ratios ---------------------------------------------------
#
# The EV1/EV2 ratio summarises how dominant the first eigenvalue is
# relative to the second in the multiple-culture condition.

eigenvalue_ratio_multiple <- eigenvalues_all %>%
  filter(
    condition == "multiple_cultures_scenario",
    component %in% c(1, 2)
  ) %>%
  select(
    sim,
    measure,
    component,
    eigenvalue
  ) %>%
  pivot_wider(
    names_from = component,
    values_from = eigenvalue,
    names_prefix = "EV"
  ) %>%
  mutate(
    EV1_EV2_ratio = EV1 / EV2
  )

eigenvalue_ratio_multiple_summary <- eigenvalue_ratio_multiple %>%
  group_by(measure) %>%
  summarise(
    mean_EV1 = mean(EV1, na.rm = TRUE),
    mean_EV2 = mean(EV2, na.rm = TRUE),
    mean_EV1_EV2_ratio = mean(EV1_EV2_ratio, na.rm = TRUE),
    sd_EV1_EV2_ratio = sd(EV1_EV2_ratio, na.rm = TRUE),
    .groups = "drop"
  )

print(eigenvalue_ratio_multiple_summary, n = Inf)


############################################################
# 10. Single-culture combined figure
#
# Panel a: mean agreement by resolution pairing.
# Panel b: eta-squared recovery of resolution from EV1 and EV2.
############################################################

p_heat_single <- ggplot(
  pair_heat,
  aes(x = resolution_i, y = resolution_j, fill = mean_agreement)
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = round(mean_agreement, 2)),
    size = single_tile_text_size
  ) +
  facet_wrap(~ measure, nrow = 1) +
  coord_equal() +
  scale_fill_gradient(
    limits = c(0, 1),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = single_base_size) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = single_strip_size
    ),
    axis.text.x = element_text(
      size = single_axis_text_size,
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      margin = margin(t = 2)
    ),
    axis.text.y = element_text(
      size = single_axis_text_size
    ),
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.1, "lines"),
    legend.position = "none"
  )

print(p_heat_single)


ev_sep_single_resolution <- ev_separation_summary %>%
  filter(condition == "single_culture_scenario") %>%
  select(
    condition,
    measure,
    EV1_eta_resolution,
    EV2_eta_resolution
  ) %>%
  pivot_longer(
    cols = c(EV1_eta_resolution, EV2_eta_resolution),
    names_to = "diagnostic",
    values_to = "eta_squared"
  ) %>%
  separate(
    diagnostic,
    into = c("component", "eta", "grouping"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    measure = factor(
      measure,
      levels = c("APA", "Pearson correlation", "Adjusted Rand Index")
    ),
    component = factor(
      component,
      levels = c("EV1", "EV2")
    )
  )

p_ev_sep_single_resolution <- ggplot(
  ev_sep_single_resolution,
  aes(x = component, y = eta_squared)
) +
  geom_col(
    fill = "grey40",
    colour = "black",
    linewidth = 0.2,
    width = 0.55
  ) +
  facet_wrap(~ measure, nrow = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Eigenvector",
    y = expression(eta^2)
  ) +
  theme_minimal(base_size = single_base_size) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = single_strip_size
    ),
    axis.title = element_text(
      size = single_axis_title_size
    ),
    axis.text = element_text(
      size = single_axis_text_size
    ),
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing.x = unit(1.1, "lines")
  )

print(p_ev_sep_single_resolution)

save_plot(
  "sim_eigenvector_separation_single_culture_APA_Pearson_ARI.png",
  p_ev_sep_single_resolution,
  width = 12,
  height = 3.2
)


p_heat_single_tagged <- p_heat_single +
  labs(title = "") +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

p_ev_sep_single_resolution_tagged <- p_ev_sep_single_resolution +
  labs(title = "") +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

p_single_culture_combined <- p_heat_single_tagged / p_ev_sep_single_resolution_tagged +
  plot_layout(
    heights = c(1.05, 1.00)
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "",
    tag_suffix = "",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 13
      )
    )
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 13
    ),
    plot.tag.position = c(0.01, 0.98)
  )

print(p_single_culture_combined)

save_plot(
  "Figure_single_culture_heatmap_and_resolution_separation.png",
  p_single_culture_combined,
  width = 12,
  height = 8.8,
  dpi = 600
)


############################################################
# 11. Multiple-culture eigenvalue spectra
############################################################

eigen_spectrum_plot_data <- eigenvalues_all %>%
  filter(
    component <= 10,
    condition == "multiple_cultures_scenario"
  ) %>%
  mutate(
    measure = factor(
      measure,
      levels = c("APA", "Pearson correlation", "Adjusted Rand Index")
    )
  ) %>%
  group_by(measure, component) %>%
  summarise(
    mean_eigenvalue = mean(eigenvalue, na.rm = TRUE),
    .groups = "drop"
  )

p_spectrum <- ggplot(
  eigen_spectrum_plot_data,
  aes(x = component, y = mean_eigenvalue, group = measure)
) +
  geom_line(linewidth = 0.45, colour = "black") +
  geom_point(size = 2.0, colour = "black") +
  facet_wrap(~ measure, nrow = 1) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "Eigenvalue rank",
    y = "Mean eigenvalue"
  ) +
  theme_minimal(base_size = eigen_base_size) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = eigen_strip_size
    ),
    axis.title = element_text(
      size = eigen_axis_title_size
    ),
    axis.text = element_text(
      size = eigen_axis_text_size
    ),
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_spectrum)

save_plot(
  "sim_eigenvalue_spectra_multiple_cultures_APA_Pearson_ARI.png",
  p_spectrum,
  width = 12,
  height = 3.8
)


############################################################
# 12. Multiple-culture eigenvector recovery
############################################################

ev_sep_long <- ev_separation_summary %>%
  filter(condition == "multiple_cultures_scenario") %>%
  pivot_longer(
    cols = starts_with("EV"),
    names_to = "diagnostic",
    values_to = "eta_squared"
  ) %>%
  separate(
    diagnostic,
    into = c("component", "eta", "grouping"),
    sep = "_",
    extra = "drop",
    remove = FALSE
  ) %>%
  mutate(
    measure = factor(
      measure,
      levels = c("APA", "Pearson correlation", "Adjusted Rand Index")
    ),
    grouping = recode(
      grouping,
      resolution = "Lumping/splitting",
      latent = "Latent organisation"
    ),
    grouping = factor(
      grouping,
      levels = c("Lumping/splitting", "Latent organisation")
    ),
    component = factor(
      component,
      levels = c("EV1", "EV2", "EV3")
    )
  ) %>%
  filter(component != "EV4")

p_ev_sep <- ggplot(
  ev_sep_long,
  aes(x = component, y = eta_squared)
) +
  geom_col(
    fill = "grey40",
    colour = "black",
    linewidth = 0.2,
    width = 0.55
  ) +
  facet_grid(grouping ~ measure) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Eigenvector",
    y = expression(eta^2)
  ) +
  theme_minimal(base_size = eigen_base_size) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = eigen_strip_size
    ),
    axis.title = element_text(
      size = eigen_axis_title_size
    ),
    axis.text = element_text(
      size = eigen_axis_text_size
    ),
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_ev_sep)

save_plot(
  "sim_eigenvector_separation_multiple_cultures_APA_Pearson_ARI.png",
  p_ev_sep,
  width = 12,
  height = 5.0
)


p_spectrum_tagged <- p_spectrum +
  labs(title = "") +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

p_ev_sep_tagged <- p_ev_sep +
  labs(title = "") +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

p_eigen_combined <- p_spectrum_tagged / p_ev_sep_tagged +
  plot_layout(
    heights = c(0.85, 1.25)
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "",
    tag_suffix = "",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 22
      )
    )
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 22
    ),
    plot.tag.position = c(0.01, 0.98)
  )

print(p_eigen_combined)

save_plot(
  "Figure_eigenvalue_spectra_and_eigenvector_separation.png",
  p_eigen_combined,
  width = 12,
  height = 8.8,
  dpi = 600
)


############################################################
# 13. Save tables
############################################################

if (SAVE_TABLES) {
  write.csv(summary_all, file.path(OUT_DIR, "sim_summary.csv"), row.names = FALSE)
  write.csv(pair_all, file.path(OUT_DIR, "sim_pairwise.csv"), row.names = FALSE)
  write.csv(meta_all, file.path(OUT_DIR, "sim_metadata.csv"), row.names = FALSE)
  write.csv(structure_summary, file.path(OUT_DIR, "sim_structure_pair_summary.csv"), row.names = FALSE)

  write.csv(eigenvalues_all, file.path(OUT_DIR, "sim_eigenvalues_APA_Pearson_ARI.csv"), row.names = FALSE)
  write.csv(eigenvectors_all, file.path(OUT_DIR, "sim_eigenvectors_APA_Pearson_ARI.csv"), row.names = FALSE)
  write.csv(ev_separation, file.path(OUT_DIR, "sim_eigenvector_separation_by_sim.csv"), row.names = FALSE)
  write.csv(ev_separation_summary, file.path(OUT_DIR, "sim_eigenvector_separation_summary.csv"), row.names = FALSE)
  write.csv(eigenvalue_ratio_multiple_summary, file.path(OUT_DIR, "sim_eigenvalue_ratio_multiple_summary.csv"), row.names = FALSE)
}