#################################################
## 06_figures_and_tables.R
##
## Figures and tables for the main text and the
## supplement, from the model summaries (step 05)
## and the per-observation effect sizes (step 02).
## Files are named by their place in the paper
## (Figure_2_..., Figure_S13_..., Table_S7_...).
##
## Every estimate is shown as an X-fold change in
## species richness (exp of the log response ratio)
## with its 80% (thick) and 95% (thin) credible
## intervals. "N = x (y)" labels give the number of
## observations (studies) each estimate is based on.
##
## Output: output/figures/*.png, output/tables/*.csv,
## and the data-source, language and off-land-site
## summaries in output/summaries/
#################################################

library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)

paths_file <- Find(file.exists, c("R/paths.R", file.path(Sys.getenv("LBMS_ANALYSIS_DIR"), "R", "paths.R")))
if(is.null(paths_file)) stop("Run from the repository root (the folder that contains run_analysis.R).")
source(paths_file)
source(file.path(code_dir, "R", "helpers.R"))
source(file.path(code_dir, "R", "models.R"))
source(file.path(code_dir, "R", "plotting.R"))

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

fold_change_label <- "X-fold change in species richness"

## Figure 2 order, top to bottom: nature-based
## strategies, then the others. Figures that facet by
## intervention follow it; y axes run bottom to top,
## so they use rev().
intervention_display_order <- c(
  "Reforestation",
  "Grassland restoration",
  "Inland wetland restoration",
  "Coastal wetland restoration",
  "Agroforestry",
  "Afforestation",
  "Bioenergy cropping",
  "Wind",
  "Solar"
)

## Factor for facets: Figure 2 order, then any other
## label (e.g. "Forestation").
order_interventions <- function(x){
  x <- as.character(x)
  factor(x, levels = c(intervention_display_order,
                       sort(setdiff(unique(x), intervention_display_order))))
}

#################################################
## Model summaries (step 05)
#################################################

read_summary <- function(name){
  df <- read.csv(file.path(summary_dir, paste0(name, ".csv")), stringsAsFactors = FALSE)
  if("intervention" %in% names(df)) df$Intervention <- plot_labels[df$intervention]
  df
}

main_effects <- read_summary("main_effects")
taxa_effects <- read_summary("all_taxa_effects")
zone_effects <- read_summary("all_zone_effects")
biome_effects <- read_summary("all_biome_effects")
structure_effects <- read_summary("all_structure_effects")
human_use_effects <- read_summary("all_human_use_effects")
experiment_type_effects <- read_summary("all_experiment_type_effects")
forestation_method_effects <- read_summary("all_forestation_method_effects")
aggregation_effects <- read_summary("all_aggregation_effects")
weighted_effects <- read_summary("all_weighted_effects")
unweighted_effects <- read_summary("all_unweighted_effects")
pet_effects <- read_summary("all_pet_effects")
peese_effects <- read_summary("all_peese_effects")
high_n_excl_effects <- read_summary("all_high_n_excl_effects")
one_per_study_effects <- read_summary("all_one_per_study_effects")
prior_sensitivity_effects <- read_summary("all_prior_sensitivity_effects")
area_curves <- read_summary("all_area_curves")
time_curves <- read_summary("all_time_curves")
all_model_coefficients <- read_summary("all_model_coefficients")
model_counts <- read_summary("all_model_counts")

## A model's slope: the fixed effect that is neither
## the intercept nor S_aggregation
is_moderator_coefficient <- function(parameter){
  parameter != "b_Intercept" & !grepl("^b_S_aggregation", parameter)
}

#################################################
## Sample-size labels
##
## n_label(): "N = observations (studies)" for one
## model and level, from all_model_counts.csv.
## add_n_labels(): those labels on a forest plot,
## left-aligned against the y axis and dodged like
## the points. n_strip_labeller(): in facet strips.
#################################################

## Same size as the axis tick labels (14 pt x 0.95)
n_label_size <- 14 * 0.95 / ggplot2::.pt

n_label <- function(intervention, model, variable = "overall", level = "overall"){
  lookup <- setNames(
    paste0("N = ", model_counts$n_observations, " (", model_counts$n_studies, ")"),
    paste(model_counts$intervention, model_counts$model, model_counts$variable,
          model_counts$level, sep = "\r")
  )
  unname(lookup[paste(intervention, model, variable, level, sep = "\r")])
}

## x = 0 is the left edge of a log2 axis (log2(0) =
## -Inf); on a linear axis pass x = -Inf. With
## colour_by_group the text takes the colour of its
## estimate.
add_n_labels <- function(df, y = "Intervention", group = NULL,
                         dodge_width = 0.4, size = n_label_size, colour_by_group = FALSE,
                         x = 0){
  df <- df[!is.na(df$n_text), , drop = FALSE]
  df$n_text <- paste0("   ", df$n_text)
  df$.n_x <- x
  aes_list <- aes(x = .n_x, y = .data[[y]], label = n_text)
  if(!is.null(group)) aes_list$group <- rlang::quo(.data[[group]])
  if(colour_by_group && !is.null(group)) aes_list$colour <- rlang::quo(.data[[group]])
  args <- list(
    data = df, mapping = aes_list, inherit.aes = FALSE,
    hjust = 0, size = size, show.legend = FALSE,
    position = if(is.null(group)) "identity" else position_dodge(width = dodge_width)
  )
  if(!(colour_by_group && !is.null(group))) args$colour <- "grey20"
  do.call(geom_text, args)
}

n_subtitle <- function(...){
  x <- n_label(...)
  if(length(x) == 0 || is.na(x)) NULL else x
}

n_strip_labeller <- function(model, variable = "overall", level = "overall"){
  ivs <- c(interventions, "forestation")
  labs <- setNames(
    paste0(plot_labels[ivs], ifelse(is.na(n_label(ivs, model, variable, level)), "",
                                    paste0("\n", n_label(ivs, model, variable, level)))),
    plot_labels[ivs]
  )
  labeller(Intervention = labs)
}

#################################################
## Per-observation data (step 02)
#################################################

## The observations the models can use
## (modelled_rows(), R/models.R). Chao, rarefied and
## other estimates, and summed richness at unequal
## sampling effort, are left out here too, so Table 1, the
## maps, the design table and the scatter points count
## the same rows as the models.
read_effects <- function(intervention){
  df <- read.csv(file.path(effects_dir, paste0(intervention, ".csv")), stringsAsFactors = FALSE)
  df[modelled_rows(df), ]
}

raw_data <- bind_rows(lapply(interventions, function(intervention){
  df <- read_effects(intervention)
  df$intervention <- intervention
  df[, intersect(c("intervention", "study_id", "observation_id", "site_latitude_dd",
                   "site_longitude_dd", "country", "taxa_grouped", "zone", "lrr_richness"),
                 names(df))]
}))

raw_data$Intervention <- plot_labels[raw_data$intervention]

#################################################
## Table helper
##
## Identifier columns followed by the 2.5, 5, 10, 50,
## 90, 95 and 97.5% posterior quantiles as fold
## changes.
#################################################

make_ci_table <- function(df, id_cols){

  df %>%
    mutate(
      `0.025` = exp(q025),
      `0.05`  = exp(q05),
      `0.1`   = exp(q10),
      `0.5`   = exp(median),
      `0.9`   = exp(q90),
      `0.95`  = exp(q95),
      `0.975` = exp(q975)
    ) %>%
    select(all_of(id_cols), `0.025`, `0.05`, `0.1`, `0.5`, `0.9`, `0.95`, `0.975`)

}

write_table <- function(df, name, dir = table_dir){
  write.csv(df, file.path(dir, paste0(name, ".csv")), row.names = FALSE)
}

#################################################
## Study maps (Figure 1, Supplementary Figures 4-12)
##
## Robinson projection; Antarctica dropped (no
## sites). Polygons crossing the antimeridian are
## split so they don't smear across the map.
#################################################

world_map <-
  ne_countries(scale = 110, returnclass = "sf") %>%
  filter(continent != "Antarctica") %>%
  sf::st_break_antimeridian(lon_0 = 0)

## Height / width of the projected map
world_map_bbox <- sf::st_bbox(sf::st_transform(world_map, "+proj=robin"))

map_aspect <- unname(
  (world_map_bbox["ymax"] - world_map_bbox["ymin"]) /
    (world_map_bbox["xmax"] - world_map_bbox["xmin"])
)

make_study_map <- function(map_points, color_values, legend_title = "Intervention"){

  ggplot() +
    geom_sf(data = world_map, fill = "grey88", colour = "black", linewidth = 0.15) +
    geom_point(
      data = map_points,
      aes(x = site_longitude_dd, y = site_latitude_dd, colour = Intervention),
      size = 1.3, alpha = 0.85
    ) +
    scale_colour_manual(values = color_values, name = legend_title) +
    coord_sf(crs = "+proj=robin", default_crs = sf::st_crs(4326), expand = FALSE) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )

}

#################################################
## Figure 1: study locations
##
## One point per intervention x study x site. Sites
## more than 50 km from land in the 1:50m Natural
## Earth layer (which includes small islands) are
## left off the maps and listed in
## map_sites_off_land.csv; they are almost all
## coordinate errors. The models use every
## observation.
#################################################

map_points_all <-
  raw_data %>%
  filter(!is.na(site_latitude_dd), !is.na(site_longitude_dd)) %>%
  distinct(intervention, study_id, site_latitude_dd, site_longitude_dd, .keep_all = TRUE)

map_land_buffer_km <- 50

map_land <-
  ne_countries(scale = 50, returnclass = "sf") %>%
  sf::st_break_antimeridian(lon_0 = 0) %>%
  sf::st_transform("+proj=robin") %>%
  sf::st_make_valid() %>%
  sf::st_union()

map_points_all$km_from_land <-
  as.numeric(
    sf::st_distance(
      sf::st_as_sf(map_points_all, coords = c("site_longitude_dd", "site_latitude_dd"), crs = 4326) %>%
        sf::st_transform("+proj=robin"),
      map_land
    )
  ) / 1000

map_points_off_land <-
  map_points_all %>%
  filter(km_from_land > map_land_buffer_km) %>%
  arrange(desc(km_from_land))

write_table(
  map_points_off_land %>%
    transmute(intervention, study_id, site_latitude_dd, site_longitude_dd,
              km_from_land = round(km_from_land)),
  "map_sites_off_land", summary_dir
)

message(nrow(map_points_off_land), " map site(s) more than ", map_land_buffer_km,
        " km from land left off the maps (see map_sites_off_land.csv).")

map_points_all <- map_points_all %>% filter(km_from_land <= map_land_buffer_km)

map_points_all$Intervention <-
  factor(plot_labels[map_points_all$intervention], levels = intervention_display_order)

figure1 <-
  make_study_map(map_points_all, intervention_colors) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.01, 0.02),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", colour = NA)
  )

save_lbms_plot(figure1, file.path(figure_dir, "Figure_1_study_locations.png"),
               width = 12, height = 12 * map_aspect)

#################################################
## Figure 2: main effects
#################################################

figure2_counts <- data.frame(
  Intervention = unname(plot_labels[interventions]),
  label = paste0("   ", n_label(interventions, "main_model"))
)

figure2 <-
  ggplot(main_effects, aes(y = Intervention, x = fold_change)) +
  geom_errorbarh(aes(xmin = exp(q025), xmax = exp(q975)),
                 height = 0, linewidth = 0.4, color = "dodgerblue4") +
  geom_errorbarh(aes(xmin = exp(q10), xmax = exp(q90)),
                 height = 0, linewidth = 1.5, color = "dodgerblue4") +
  geom_point(size = 3, color = "dodgerblue4") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  ## N labels at the left edge of the panel
  geom_text(data = figure2_counts, aes(x = 0, y = Intervention, label = label),
            inherit.aes = FALSE, hjust = 0, size = n_label_size) +
  scale_y_discrete(limits = rev(intervention_display_order)) +
  scale_fold("x", L = 2) +
  theme_lbms() +
  labs(x = fold_change_label, y = NULL, colour = NULL)

save_lbms_plot(figure2, file.path(figure_dir, "Figure_2_main_effects.png"), width = 9, height = 6)

#################################################
## Moderator forest plot
##
## One row per intervention, one estimate per
## moderator level. For two-level moderators,
## annotate_difference marks the difference between
## the levels: "*" strong support, "." moderate
## support. For multi-level moderators,
## show_cld_letters adds compact letters (levels that
## share a letter do not differ with strong support);
## they are omitted for a row whose levels all share
## one letter.
#################################################

make_category_forest_plot <- function(df, filename, colors, y_limits,
                                      annotate_difference = FALSE,
                                      show_cld_letters = FALSE,
                                      count_model = NULL,
                                      count_variable = "overall"){

  plot_df <- df %>% filter(level != "Difference", Intervention %in% y_limits)

  ## Wider dodge when rows hold several levels
  n_levels <- max(1, n_distinct(plot_df$level))
  dodge_w <- if(n_levels > 2) 0.75 else 0.5

  p <-
    ggplot(plot_df, aes(y = Intervention, x = fold_change, colour = level)) +
    geom_errorbarh(aes(xmin = exp(q025), xmax = exp(q975)),
                   height = 0, linewidth = 0.4, position = position_dodge(width = dodge_w)) +
    geom_errorbarh(aes(xmin = exp(q10), xmax = exp(q90)),
                   height = 0, linewidth = 1.4, position = position_dodge(width = dodge_w)) +
    geom_point(size = 2.5, position = position_dodge(width = dodge_w)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_y_discrete(limits = y_limits) +
    scale_fold("x") +
    scale_color_manual(values = colors) +
    theme_lbms()

  if(!is.null(count_model)){
    plot_df$n_text <- n_label(
      plot_df$intervention, count_model, count_variable,
      if(count_variable == "overall") "overall" else as.character(plot_df$level)
    )
    p <- p + add_n_labels(plot_df, group = "level", dodge_width = dodge_w, colour_by_group = TRUE)
  }

  if(annotate_difference){

    diff_df <-
      df %>%
      filter(level == "Difference", Intervention %in% y_limits) %>%
      mutate(annotation = ifelse(strong_support, "*", ifelse(moderate_support, ".", ""))) %>%
      filter(annotation != "")

    if(nrow(diff_df) > 0){
      ## At the right edge, so the marks don't widen the axis
      p <- p + geom_text(data = diff_df, aes(y = Intervention, x = Inf, label = annotation),
                         inherit.aes = FALSE, colour = "black", size = 5, hjust = 1.5)
    }

  }

  if(show_cld_letters){

    ## One column of letters per row, right of the row's
    ## widest 95% CI
    row_max_upper <-
      plot_df %>%
      group_by(Intervention) %>%
      summarise(row_max_upper = max(exp(q975), na.rm = TRUE), .groups = "drop")

    cld_df <-
      plot_df %>%
      filter(!is.na(cld_letter) & cld_letter != "") %>%
      group_by(Intervention) %>%
      filter(n_distinct(cld_letter) > 1) %>%
      ungroup() %>%
      left_join(row_max_upper, by = "Intervention")

    if(nrow(cld_df) > 0){
      p <- p + geom_text(data = cld_df, aes(x = row_max_upper * 1.15, label = cld_letter),
                         position = position_dodge(width = dodge_w),
                         hjust = 1, size = 3.5, fontface = "bold", show.legend = FALSE)
    }

  }

  save_lbms_plot(p, file.path(figure_dir, filename), width = 12,
                 height = max(8, 1.5 + 0.32 * n_distinct(plot_df$Intervention) * n_levels))

}

#################################################
## Figure 3
##
## A: forested vs open biome for the nature-based
##    interventions (re- and afforestation pooled as
##    Forestation).
## B: forestation method x biome, with compact
##    letters.
## C: tropical vs temperate, forestation.
## A and C mark the difference between the two levels
## with "*" (strong support) or a square (moderate
## support). B and C share one x range.
#################################################

biome_effects_fig3a <-
  biome_effects %>%
  filter(intervention %in% c("grassland_restoration", "wetland_restoration-inland",
                             "wetland_restoration-coastal", "agroforestry", "forestation"))

zone_effects_fig3c <- zone_effects %>% filter(intervention == "forestation")

## Figure 3B method labels -> forest_type values (for
## the N labels)
fig3b_method_key <- c(
  "Natural forest regrowth" = "natural forest restoration",
  "Natural\nforest\nregrowth" = "natural forest restoration",
  "Plantation" = "timber plantation"
)

fig3b_n_text <- function(method, level){
  n_label("forestation", "forestation_type_model", "forest_type | forested",
          paste(fig3b_method_key[as.character(method)], as.character(level), sep = " | "))
}

## 95% CI (thin), 80% CI (thick), median, dotted line at 1
fig3_layers <- function(dodge_width){
  list(
    geom_vline(xintercept = 1, linetype = "dotted"),
    geom_errorbarh(aes(xmin = exp(q025), xmax = exp(q975)),
                   height = 0, linewidth = 0.4, position = position_dodge(width = dodge_width)),
    geom_errorbarh(aes(xmin = exp(q10), xmax = exp(q90)),
                   height = 0, linewidth = 1.4, position = position_dodge(width = dodge_width)),
    geom_point(size = 2.5, position = position_dodge(width = dodge_width))
  )
}

## Panel letter, top left inside the panel
fig3_tag <- function(tag){
  annotate("text", x = 0, y = Inf, label = tag, hjust = -0.3, vjust = 1.3, size = 5.5, fontface = "bold")
}

fig3_support_mark <- function(strong, moderate){
  ifelse(strong, "*", ifelse(moderate, "▪", ""))
}

if(nrow(biome_effects_fig3a) > 0 && nrow(forestation_method_effects) > 0 && nrow(zone_effects_fig3c) > 0){

  ## ---- Panel A ----

  fig3a_rows <- c("Forestation", "Agroforestry", "Coastal wetland restoration",
                  "Inland wetland restoration", "Grassland restoration")

  fig3a_df <-
    biome_effects_fig3a %>%
    filter(level != "Difference") %>%
    mutate(level = factor(as.character(level), levels = c("Forested", "Open")))

  ## One column of marks, right of the widest 95% CI
  ## among the marked rows; "*" is lowered to row centre
  fig3a_marks <-
    biome_effects_fig3a %>%
    filter(level == "Difference") %>%
    mutate(mark = fig3_support_mark(strong_support, moderate_support),
           mark_vjust = ifelse(mark == "*", 0.75, 0.5)) %>%
    filter(mark != "")

  fig3a_mark_x <-
    fig3a_df %>%
    filter(Intervention %in% fig3a_marks$Intervention) %>%
    summarise(x = max(exp(q975))) %>%
    pull(x) * 1.12

  fig3a_panel <-
    ggplot(fig3a_df, aes(y = Intervention, x = fold_change, colour = level)) +
    fig3_layers(0.4) +
    geom_text(data = fig3a_marks,
              aes(y = Intervention, x = fig3a_mark_x, label = mark, vjust = mark_vjust),
              inherit.aes = FALSE, size = 8, hjust = 0) +
    fig3_tag("A.") +
    add_n_labels(
      fig3a_df %>% mutate(n_text = n_label(intervention, "biome_model", "forested", as.character(level))),
      group = "level", dodge_width = 0.4, colour_by_group = TRUE
    ) +
    scale_y_discrete(limits = fig3a_rows) +
    scale_fold("x") +
    scale_colour_manual(values = biome_colors, name = "Biome type") +
    theme_lbms() +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.99, 0.99),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4)
    ) +
    labs(x = fold_change_label, y = NULL)

  ## ---- Panels B and C ----

  fig3b_df <-
    forestation_method_effects %>%
    mutate(
      level = factor(as.character(level), levels = c("Forested", "Open")),
      method = factor(as.character(method), levels = c("Natural forest regrowth", "Plantation"),
                      labels = c("Natural\nforest\nregrowth", "Plantation"))
    )

  fig3b_letters <-
    fig3b_df %>%
    group_by(method) %>%
    mutate(letter_x = max(exp(q975)) * 1.06) %>%
    ungroup()

  fig3c_df <-
    zone_effects_fig3c %>%
    filter(level != "Difference") %>%
    mutate(level = factor(as.character(level), levels = c("Temperate", "Tropical")))

  fig3c_mark <-
    zone_effects_fig3c %>%
    filter(level == "Difference") %>%
    mutate(mark = fig3_support_mark(strong_support, moderate_support)) %>%
    filter(mark != "")

  fig3c_mark_x <- max(exp(fig3c_df$q975)) * 1.08

  ## One symmetric range for B and C, letters and marks included
  fig3bc_L <- fold_axis_limit(c(
    exp(fig3b_df$q025), exp(fig3b_df$q975), fig3b_letters$letter_x * 1.05,
    exp(fig3c_df$q025), exp(fig3c_df$q975), fig3c_mark_x * 1.05
  ))

  fig3bc_x_scale <- scale_fold("x", L = fig3bc_L)

  fig3b_panel <-
    ggplot(fig3b_df, aes(y = method, x = fold_change, colour = level)) +
    fig3_layers(0.4) +
    geom_text(data = fig3b_letters, aes(x = letter_x, label = cld_letter, group = level),
              position = position_dodge(width = 0.4),
              hjust = 0, size = 4.5, fontface = "bold", colour = "black", show.legend = FALSE) +
    fig3_tag("B.") +
    add_n_labels(fig3b_df %>% mutate(n_text = fig3b_n_text(method, level)),
                 y = "method", group = "level", dodge_width = 0.4, colour_by_group = TRUE) +
    fig3bc_x_scale +
    scale_colour_manual(values = biome_colors) +
    theme_lbms() +
    theme(legend.position = "none") +
    labs(x = fold_change_label, y = NULL)

  fig3c_panel <-
    ggplot(fig3c_df, aes(y = level, x = fold_change, colour = level)) +
    fig3_layers(0) +
    fig3_tag("C.") +
    add_n_labels(fig3c_df %>% mutate(n_text = n_label("forestation", "zone_model", "zone", as.character(level))),
                 y = "level") +
    fig3bc_x_scale +
    scale_colour_manual(values = zone_colors) +
    theme_lbms() +
    theme(legend.position = "none") +
    labs(x = fold_change_label, y = NULL)

  ## The difference mark sits between the two rows
  if(nrow(fig3c_mark) > 0){
    fig3c_panel <- fig3c_panel +
      annotate("text", x = fig3c_mark_x, y = 1.5, label = fig3c_mark$mark[1], size = 8, hjust = 0,
               vjust = ifelse(fig3c_mark$mark[1] == "*", 0.75, 0.5))
  }

  ## ---- Assemble ----

  fig3_right <-
    cowplot::plot_grid(
      cowplot::ggdraw() + cowplot::draw_label("Forestation", fontface = "bold", size = 16),
      cowplot::plot_grid(fig3b_panel, fig3c_panel, ncol = 1, align = "v", axis = "lr"),
      ncol = 1,
      rel_heights = c(0.07, 1)
    )

  figure3 <- cowplot::plot_grid(fig3a_panel, fig3_right, ncol = 2, rel_widths = c(1.45, 1))

  save_lbms_plot(figure3, file.path(figure_dir, "Figure_3_biome_and_zone.png"), width = 14, height = 7)

}

#################################################
## Figures 4 and 5: area and age
##
## Posterior median and 80% / 95% bands of the slope
## model across each intervention's observed range,
## over the observations (grey). Blue solid line:
## strong support for the slope; grey dashed: not.
## The y axis is 1/4 to 4 in every panel (points
## beyond it are drawn at the edge).
#################################################

slope_support <- function(model_name){
  all_model_coefficients %>%
    filter(model == model_name, is_moderator_coefficient(parameter)) %>%
    select(intervention, strong_support)
}

add_significance_group <- function(curves, model_name){
  curves %>%
    left_join(slope_support(model_name), by = "intervention") %>%
    mutate(significance_group = ifelse(!is.na(strong_support) & strong_support,
                                       "Significant", "Not significant"))
}

slope_curve_layers <- function(){
  list(
    geom_ribbon(aes(ymin = exp(lower95), ymax = exp(upper95), fill = significance_group), alpha = 0.15),
    geom_ribbon(aes(ymin = exp(lower80), ymax = exp(upper80), fill = significance_group), alpha = 0.30),
    geom_line(aes(colour = significance_group, linetype = significance_group), linewidth = 1),
    geom_hline(yintercept = 1, linetype = "dashed")
  )
}

slope_curve_scales <- function(){
  list(
    scale_fold("y", L = 4),
    scale_linetype_manual(name = "Slope support",
                          values = c("Significant" = "solid", "Not significant" = "dashed")),
    scale_color_manual(name = "Slope support", values = slope_support_values),
    scale_fill_manual(name = "Slope support", values = slope_support_values)
  )
}

if(nrow(area_curves) > 0){

  area_curves <- add_significance_group(area_curves, "area_model")

  area_scatter <- bind_rows(lapply(unique(area_curves$intervention), function(iv){
    d <- read_effects(iv)
    d <- d[!is.na(d$lrr_richness) & !is.na(d$mean_size_ha) & d$mean_size_ha > 0, ]
    if(nrow(d) == 0) return(NULL)
    data.frame(intervention = iv, log_area = log(d$mean_size_ha), fold_change = exp(d$lrr_richness))
  }))

  area_scatter$Intervention <- plot_labels[area_scatter$intervention]

  ## Layer order: scatter, 95% band, 80% band, line, y = 1
  figure4 <-
    ggplot(area_curves, aes(x = exp(predictor), y = exp(median))) +
    geom_point(data = area_scatter, aes(x = exp(log_area), y = fold_change),
               inherit.aes = FALSE, colour = "grey60", alpha = 0.25, size = 0.7) +
    slope_curve_layers() +
    scale_x_log10() +
    slope_curve_scales() +
    facet_wrap(vars(Intervention = order_interventions(Intervention)), scales = "free_x", ncol = 3,
               labeller = n_strip_labeller("area_model")) +
    theme_lbms() +
    theme(legend.position = "bottom") +
    labs(x = "Intervention area (ha)", y = "Biodiversity fold-change")

  save_lbms_plot(figure4, file.path(figure_dir, "Figure_4_intervention_area.png"), width = 12, height = 8)

}

if(nrow(time_curves) > 0){

  time_curves <- add_significance_group(time_curves, "time_model")

  time_scatter <- bind_rows(lapply(unique(time_curves$intervention), function(iv){
    d <- read_effects(iv)
    d <- d[!is.na(d$lrr_richness) & !is.na(d$intervention_length), ]
    if(nrow(d) == 0) return(NULL)
    data.frame(intervention = iv, years = d$intervention_length, fold_change = exp(d$lrr_richness))
  }))

  time_scatter$Intervention <- plot_labels[time_scatter$intervention]

  ## Each panel runs from 0 to 1.1 x its oldest site,
  ## capped at 100 years; points and curve beyond the
  ## cap are dropped from the panel.
  time_axis_caps <-
    time_scatter %>%
    group_by(Intervention) %>%
    summarise(cap = pmin(1.1 * max(years, na.rm = TRUE), 100), .groups = "drop")

  time_scatter <- time_scatter %>% left_join(time_axis_caps, by = "Intervention") %>% filter(years <= cap)
  time_curves <- time_curves %>% left_join(time_axis_caps, by = "Intervention") %>% filter(predictor <= cap)

  time_axis_anchors <- time_scatter %>% distinct(Intervention) %>% mutate(x = 0)

  figure5 <-
    ggplot(time_curves, aes(x = predictor, y = exp(median))) +
    geom_point(data = time_scatter, aes(x = years, y = fold_change),
               inherit.aes = FALSE, colour = "grey60", alpha = 0.25, size = 0.7) +
    geom_blank(data = time_axis_anchors, aes(x = x), inherit.aes = FALSE) +
    slope_curve_layers() +
    slope_curve_scales() +
    ## x axis starts exactly at 0
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    facet_wrap(vars(Intervention = order_interventions(Intervention)), scales = "free_x", ncol = 3,
               labeller = n_strip_labeller("time_model")) +
    theme_lbms() +
    theme(legend.position = "bottom") +
    labs(x = "Years since intervention", y = "Biodiversity fold-change")

  save_lbms_plot(figure5, file.path(figure_dir, "Figure_5_intervention_age.png"), width = 12, height = 8)

}

#################################################
## Moderator figures (supplement)
#################################################

all_rows <- c("Forestation", rev(intervention_display_order))

if(nrow(taxa_effects) > 0){

  taxa_effects$level <- factor(taxa_effects$level, levels = taxa_groups)

  make_category_forest_plot(
    taxa_effects, "Figure_S16_taxonomic_group.png", taxa_colors,
    y_limits = rev(intervention_display_order),
    show_cld_letters = TRUE,
    count_model = "taxa_model", count_variable = "taxa_grouped"
  )

}

zone_effects$level <- factor(zone_effects$level, levels = c("Temperate", "Tropical", "Difference"))

make_category_forest_plot(
  zone_effects, "Figure_S17_climate_zone.png", zone_colors,
  y_limits = all_rows,
  annotate_difference = TRUE,
  count_model = "zone_model", count_variable = "zone"
)

structure_effects$level <- factor(structure_effects$level, levels = c("Woody", "Open", "Difference"))

make_category_forest_plot(
  structure_effects, "Figure_S14_reference_structure.png", structure_colors,
  y_limits = all_rows,
  annotate_difference = TRUE,
  count_model = "structure_model", count_variable = "reference_habitat_structure_axis"
)

if(nrow(human_use_effects) > 0){

  human_use_effects$level <-
    factor(human_use_effects$level, levels = c("Intact natural", "Degraded natural", "Human use"))

  make_category_forest_plot(
    human_use_effects, "Figure_S15_reference_human_use.png", human_use_colors,
    y_limits = all_rows,
    show_cld_letters = TRUE,
    count_model = "human_use_model", count_variable = "reference_habitat_human_use_axis"
  )

}

## Solar and wind are entirely observational, so
## they have a Natural estimate only
experiment_type_effects$level <-
  factor(experiment_type_effects$level, levels = c("Natural", "Manipulative", "Difference"))

make_category_forest_plot(
  experiment_type_effects, "Figure_S18_study_design.png", experiment_type_colors,
  y_limits = all_rows,
  annotate_difference = TRUE,
  count_model = "experiment_type_model", count_variable = "experiment_type"
)

#################################################
## Sensitivity comparisons (supplement)
##
## Main-effect intercepts under alternative analyses,
## side by side for each intervention.
#################################################

make_comparison_plot <- function(plot_df, method_levels, colors, dodge, filename, width, height){

  plot_df$Method <- factor(plot_df$Method, levels = method_levels)

  p <-
    ggplot(plot_df, aes(y = Intervention, x = fold_change, colour = Method)) +
    geom_errorbarh(aes(xmin = exp(q025), xmax = exp(q975)),
                   height = 0, linewidth = 0.4, position = position_dodge(width = dodge)) +
    geom_errorbarh(aes(xmin = exp(q10), xmax = exp(q90)),
                   height = 0, linewidth = 1.4, position = position_dodge(width = dodge)) +
    geom_point(size = 2, position = position_dodge(width = dodge)) +
    add_n_labels(plot_df %>% mutate(n_text = n_label(intervention, model)),
                 group = "Method", dodge_width = dodge, colour_by_group = TRUE) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_y_discrete(limits = rev(intervention_display_order)) +
    scale_fold("x") +
    scale_color_manual(values = colors) +
    theme_lbms() +
    labs(x = fold_change_label, y = NULL, colour = "Method")

  save_lbms_plot(p, file.path(figure_dir, filename), width = width, height = height)

}

intercepts <- function(df, method) df %>% filter(parameter == "b_Intercept") %>% mutate(Method = method)

## Aggregation method: pooled (main model) vs studies
## reporting mean or aggregated richness
if(nrow(aggregation_effects) > 0){
  make_comparison_plot(
    bind_rows(
      intercepts(main_effects, "Pooled"),
      intercepts(aggregation_effects, NA) %>%
        mutate(Method = ifelse(model == "mean_aggregation_model", "Mean", "Aggregated"))
    ),
    c("Aggregated", "Mean", "Pooled"), aggregation_colors, 0.75,
    "Figure_S20_aggregation.png", 10, 12.5
  )
}

## Variance weighting, on the observations with a
## standard error
if(nrow(weighted_effects) > 0 && nrow(unweighted_effects) > 0){
  make_comparison_plot(
    bind_rows(intercepts(weighted_effects, "Weighted"), intercepts(unweighted_effects, "Unweighted")),
    c("Unweighted", "Weighted"), weighting_colors, 0.4,
    "Figure_S21_variance_weighting.png", 9, 7
  )
}

## Publication bias: main model (all observations),
## variance-weighted (observations with a standard
## error), and the PEESE bias-adjusted intercept on
## those same observations
if(nrow(peese_effects) > 0){
  make_comparison_plot(
    bind_rows(
      intercepts(main_effects, "Main model"),
      intercepts(weighted_effects, "Variance-weighted"),
      intercepts(peese_effects, "PEESE")
    ),
    c("PEESE", "Variance-weighted", "Main model"), publication_bias_colors, 0.75,
    "Figure_S24_publication_bias.png", 10, 12.5
  )
}

if(nrow(high_n_excl_effects) > 0){
  make_comparison_plot(
    bind_rows(intercepts(main_effects, "Full sample"), intercepts(high_n_excl_effects, "Excl. top 10% studies")),
    c("Excl. top 10% studies", "Full sample"), high_n_colors, 0.4,
    "Figure_S25_high_n_exclusion.png", 9, 7
  )
}

if(nrow(one_per_study_effects) > 0){
  make_comparison_plot(
    bind_rows(intercepts(main_effects, "Full sample"), intercepts(one_per_study_effects, "One effect per study")),
    c("One effect per study", "Full sample"), one_per_study_colors, 0.4,
    "Figure_S26_one_effect_per_study.png", 9, 7
  )
}

if(nrow(prior_sensitivity_effects) > 0){
  make_comparison_plot(
    bind_rows(
      intercepts(main_effects, "Baseline prior"),
      intercepts(prior_sensitivity_effects, NA) %>%
        mutate(Method = ifelse(model == "prior_tight_model", "Tight prior", "Wide prior"))
    ),
    c("Wide prior", "Tight prior", "Baseline prior"), prior_sensitivity_colors, 0.75,
    "Figure_S19_prior_sensitivity.png", 10, 12.5
  )
}

#################################################
## PET slope (funnel-plot asymmetry)
##
## Slope of each observation's effect size on its
## standard error; blue where the 95% CI excludes
## zero.
#################################################

pet_slope_effects <-
  pet_effects %>%
  filter(parameter == "b_lrr_se_richness") %>%
  mutate(significance_group = ifelse(!is.na(strong_support) & strong_support,
                                     "Significant", "Not significant"))

if(nrow(pet_slope_effects) > 0){

  figure_pet_slope <-
    ggplot(pet_slope_effects, aes(y = Intervention, x = median, colour = significance_group)) +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0, linewidth = 0.4) +
    geom_errorbarh(aes(xmin = q10, xmax = q90), height = 0, linewidth = 1.4) +
    geom_point(size = 2.5) +
    add_n_labels(pet_slope_effects %>% mutate(n_text = n_label(intervention, "pet_model")), x = -Inf) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ## room on the left for the N labels
    scale_x_continuous(expand = expansion(mult = c(0.45, 0.05))) +
    scale_y_discrete(limits = rev(intervention_display_order)) +
    scale_color_manual(name = "Slope support", values = slope_support_values) +
    theme_lbms() +
    labs(x = "PET slope (funnel-plot asymmetry - 0 = no asymmetry)", y = NULL)

  save_lbms_plot(figure_pet_slope, file.path(figure_dir, "Figure_S23_PET_slopes.png"),
                 width = 9, height = 6)

}

#################################################
## Funnel plots
##
## Observations with a standard error (the PET-PEESE
## data): fold change against standard error
## (reversed), one panel per intervention, all on
## fold change 1/4 to 4 and SE 0 to 1 (points beyond
## are drawn at the edge as open triangles). Solid
## line: the variance-weighted estimate; shaded: the
## pseudo-95% funnel around it (+/- 1.96 SE); dashed:
## no effect.
#################################################

funnel_scatter <- bind_rows(lapply(interventions, function(iv){
  d <- read_effects(iv)
  d <- d[!is.na(d$lrr_richness) & !is.na(d$lrr_se_richness) & d$lrr_se_richness > 0, ]
  if(nrow(d) == 0) return(NULL)
  data.frame(intervention = iv, fold_change = exp(d$lrr_richness), se = d$lrr_se_richness)
}))

if(nrow(funnel_scatter) > 0){

  funnel_scatter$Intervention <- plot_labels[funnel_scatter$intervention]

  funnel_x_limit <- 4
  funnel_se_limit <- 1

  funnel_scatter <-
    funnel_scatter %>%
    mutate(
      beyond = se > funnel_se_limit | fold_change > funnel_x_limit | fold_change < 1 / funnel_x_limit,
      se_plot = pmin(se, funnel_se_limit),
      fold_plot = pmin(pmax(fold_change, 1 / funnel_x_limit), funnel_x_limit),
      Point = ifelse(beyond, "Beyond axis range (shown at edge)", "Observation")
    )

  funnel_centre <-
    weighted_effects %>%
    filter(parameter == "b_Intercept", intervention %in% unique(funnel_scatter$intervention)) %>%
    transmute(intervention, Intervention, centre_log = median)

  funnel_boundary <-
    funnel_centre %>%
    tidyr::crossing(se = seq(0, funnel_se_limit, length.out = 100)) %>%
    mutate(lower = exp(centre_log - 1.96 * se), upper = exp(centre_log + 1.96 * se))

  figure_funnel <-
    ggplot() +
    geom_ribbon(data = funnel_boundary, aes(y = se, xmin = lower, xmax = upper),
                fill = "grey85", alpha = 0.6, orientation = "y") +
    geom_path(data = funnel_boundary, aes(x = lower, y = se),
              linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    geom_path(data = funnel_boundary, aes(x = upper, y = se),
              linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_vline(data = funnel_centre, aes(xintercept = exp(centre_log)), colour = "black", linewidth = 0.5) +
    geom_point(data = funnel_scatter, aes(x = fold_plot, y = se_plot, shape = Point),
               colour = "grey25", alpha = 0.45, size = 1) +
    scale_fold("x", L = funnel_x_limit) +
    scale_y_reverse(limits = c(funnel_se_limit, 0), oob = scales::oob_keep) +
    scale_shape_manual(values = c("Observation" = 16, "Beyond axis range (shown at edge)" = 2), name = NULL) +
    facet_wrap(vars(Intervention = order_interventions(Intervention)), ncol = 3,
               labeller = n_strip_labeller("pet_model")) +
    theme_lbms() +
    theme(legend.position = "bottom") +
    labs(x = fold_change_label, y = "Standard error")

  save_lbms_plot(figure_funnel, file.path(figure_dir, "Figure_S22_funnel_plots.png"),
                 width = 12, height = 10)

  #################################################
  ## Observations outside the funnel, by side
  ##
  ## Same observations, centre and funnel as the
  ## funnel plots, on the unclipped values. Without
  ## small-study asymmetry the observations outside
  ## the funnel should split evenly between the two
  ## sides.
  #################################################

  ## An intervention whose variance-weighted model was
  ## excluded (not converged) has no funnel centre and
  ## is left out of this table.
  no_centre <- setdiff(unique(funnel_scatter$intervention), funnel_centre$intervention)
  if(length(no_centre) > 0){
    log_message("WARNING: no variance-weighted estimate (funnel centre) for ",
                paste(no_centre, collapse = ", "), " - left out of Table_S15_funnel_outside")
  }

  funnel_sides <-
    funnel_scatter %>%
    filter(intervention %in% funnel_centre$intervention) %>%
    left_join(funnel_centre %>% select(intervention, centre_log), by = "intervention") %>%
    mutate(
      z = (log(fold_change) - centre_log) / se,
      side = case_when(z < -1.96 ~ "left", z > 1.96 ~ "right", TRUE ~ "inside")
    )

  summarise_sides <- function(df){
    df %>%
      summarise(
        n_observations = n(),
        n_left = sum(side == "left"),
        n_inside = sum(side == "inside"),
        n_right = sum(side == "right"),
        .groups = "drop"
      ) %>%
      mutate(
        percent_left = round(100 * n_left / n_observations, 1),
        percent_right = round(100 * n_right / n_observations, 1),
        percent_outside = round(100 * (n_left + n_right) / n_observations, 1),
        right_share_of_outside = ifelse(n_left + n_right > 0, round(n_right / (n_left + n_right), 2), NA)
      )
  }

  write_table(
    bind_rows(
      funnel_sides %>% group_by(Intervention) %>% summarise_sides(),
      funnel_sides %>% summarise_sides() %>% mutate(Intervention = "All")
    ) %>%
      mutate(Intervention = factor(Intervention, levels = c(intervention_display_order, "All"))) %>%
      arrange(Intervention),
    "Table_S15_funnel_outside"
  )

}

#################################################
## Tropical vs temperate observation counts
#################################################

zone_counts <-
  raw_data %>%
  filter(!is.na(lrr_richness), !is.na(zone)) %>%
  count(intervention, Intervention, zone) %>%
  pivot_wider(names_from = zone, values_from = n, values_fill = 0)

figure_zone_counts <-
  ggplot(zone_counts, aes(x = Temperate, y = Tropical, colour = Intervention)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_lbms() +
  theme_lbms() +
  labs(x = "Number of temperate observations", y = "Number of tropical observations",
       colour = "Intervention")

save_lbms_plot(figure_zone_counts, file.path(figure_dir, "Figure_S2_zone_counts.png"),
               width = 8, height = 6)

#################################################
## Absolute latitude of observations
#################################################

figure_latitude <-
  raw_data %>%
  filter(!is.na(site_latitude_dd)) %>%
  mutate(abs_latitude = abs(site_latitude_dd)) %>%
  ggplot(aes(x = abs_latitude, fill = Intervention)) +
  geom_density(colour = "black", linewidth = 0.3) +
  scale_fill_lbms() +
  facet_wrap(vars(Intervention = order_interventions(Intervention)), ncol = 3) +
  theme_lbms() +
  theme(legend.position = "none") +
  labs(x = "Absolute latitude", y = "Density")

save_lbms_plot(figure_latitude, file.path(figure_dir, "Figure_S3_latitude.png"),
               width = 10, height = 8)

#################################################
## Study locations by intervention (Supplementary
## Figures 4-12, in alphabetical order of the
## intervention key)
#################################################

supp_map_numbers <- setNames(4:12, sort(interventions))

for(intervention in sort(interventions)){

  map_points_iv <- map_points_all %>% filter(intervention == !!intervention)

  if(nrow(map_points_iv) == 0) next

  label <- plot_labels[[intervention]]

  p <- make_study_map(map_points_iv, setNames(intervention_colors[label], label)) +
    theme(legend.position = "none")

  save_lbms_plot(
    p,
    file.path(figure_dir, paste0("Figure_S", supp_map_numbers[[intervention]], "_locations_",
                                 gsub("[^a-z]+", "_", tolower(label)), ".png")),
    width = 10, height = 10 * map_aspect
  )

}

#################################################
## Table 1: interventions, reference habitats and
## sample sizes (observations with an effect size)
#################################################

reference_habitats <- c(
  reforestation = "Agricultural field, pasture, or timber clear cut",
  grassland_restoration = "Degraded grassland or prairie, agricultural field, or pasture",
  `wetland_restoration-inland` = "Degraded peatland, bog, swamp, or other wetland",
  `wetland_restoration-coastal` = "Degraded mangrove, saltmarsh, or other wetland",
  agroforestry = "Agricultural or grazing land under traditional management",
  aforestation = "Agricultural field, pasture, grassland, shrubland, or wetland",
  beccs = "Traditionally managed cropland or plantation",
  wind = "Natural habitat or agricultural land without wind turbines",
  solar = "Natural habitat or agricultural land without solar panels"
)

table1 <-
  raw_data %>%
  filter(!is.na(lrr_richness)) %>%
  group_by(intervention, Intervention) %>%
  summarise(`Number of studies` = n_distinct(study_id),
            `Number of pair-wise comparisons` = n(), .groups = "drop") %>%
  mutate(`Reference habitat` = reference_habitats[intervention]) %>%
  select(Intervention, `Reference habitat`, `Number of studies`, `Number of pair-wise comparisons`) %>%
  arrange(match(Intervention, intervention_display_order))

table1 <- bind_rows(
  table1,
  data.frame(
    Intervention = "Total",
    `Reference habitat` = "",
    `Number of studies` = sum(table1$`Number of studies`),
    `Number of pair-wise comparisons` = sum(table1$`Number of pair-wise comparisons`),
    check.names = FALSE
  )
)

write_table(table1, "Table_1_interventions")

#################################################
## Study design: spatial vs temporal comparison
##
## From the extractors' free-text spatial_temporal
## field, for the Table 1 observations. Spatial: the
## intervention is compared with a separate reference
## site (including chronosequences and space-for-time
## substitutions), even when the same plots were
## resampled over time. Temporal: the same site
## before and after the intervention. Both: a
## BACI-type design, as recorded by the extractor.
## "Same as obs N" takes that observation's design.
## A study is "Both" if any observation is, or if it
## has spatial and temporal observations. Totals sum
## across interventions, as in Table 1. Every
## observation has a recorded design, so the "Not
## recorded" column is written only if one is missing.
#################################################

classify_spatial_temporal <- function(x){
  x <- tolower(trimws(x))
  design <- rep(NA_character_, length(x))
  design[grepl("^spatial|chronosequence|space-for-time|^paired|^same[- ]site|^same spatial", x)] <- "Spatial"
  design[grepl("^temporal|^before|^pre-treatment|^same as obs[0-9 ]*\\(temporal", x)] <- "Temporal"
  design[grepl("^both|spatial and temporal|spatial \\+ temporal|spatial-temporal|spatial/temporal|baci", x)] <- "Both"
  design
}

design_obs <- bind_rows(lapply(interventions, function(iv){
  d <- read_effects(iv)
  d <- d[!is.na(d$lrr_richness), ]
  d$design <- classify_spatial_temporal(d$spatial_temporal)
  same_as <- as.integer(sub(".*obs ?([0-9]+).*", "\\1",
                            ifelse(grepl("obs ?[0-9]+", tolower(d$spatial_temporal)),
                                   tolower(d$spatial_temporal), NA)))
  for(i in which(is.na(d$design) & !is.na(same_as))){
    j <- which(d$study_id == d$study_id[i] & d$observation_id == same_as[i])
    if(length(j) > 0) d$design[i] <- d$design[j[1]]
  }
  data.frame(Intervention = plot_labels[[iv]], study_id = d$study_id,
             design = ifelse(is.na(d$design), "Not recorded", d$design), stringsAsFactors = FALSE)
}))

study_design <- function(design){
  design <- design[design != "Not recorded"]
  if(length(design) == 0) return("Not recorded")
  if("Both" %in% design || all(c("Spatial", "Temporal") %in% design)) return("Both")
  design[1]
}

design_levels <- c("Spatial", "Temporal", "Both", "Not recorded")

design_studies <- design_obs %>%
  group_by(Intervention, study_id) %>%
  summarise(design = study_design(design), .groups = "drop") %>%
  count(Intervention, design) %>%
  pivot_wider(names_from = design, values_from = n, values_fill = 0)

design_comparisons <- design_obs %>%
  count(Intervention, design) %>%
  pivot_wider(names_from = design, values_from = n, values_fill = 0)

for(level in design_levels){
  if(!level %in% names(design_studies)) design_studies[[level]] <- 0L
  if(!level %in% names(design_comparisons)) design_comparisons[[level]] <- 0L
}

design_table <- design_studies %>%
  select(Intervention, all_of(design_levels)) %>%
  left_join(design_comparisons %>% select(Intervention, all_of(design_levels)),
            by = "Intervention", suffix = c("", "_comparisons")) %>%
  arrange(match(Intervention, intervention_display_order))

design_table <- bind_rows(design_table,
                          design_table %>% summarise(across(-Intervention, sum)) %>% mutate(Intervention = "Total"))

design_unrecorded <- sum(design_table$`Not recorded_comparisons`)

design_table <- design_table %>%
  transmute(
    Intervention,
    `Number of studies` = Spatial + Temporal + Both + `Not recorded`,
    `Spatial (control-impact)` = Spatial,
    `Temporal (before-after)` = Temporal,
    `Both (BACI-type or mixed)` = Both,
    `Not recorded` = `Not recorded`,
    `Comparisons: spatial / temporal / both` =
      paste(Spatial_comparisons, Temporal_comparisons, Both_comparisons, sep = " / "),
    `Comparisons: spatial / temporal / both / not recorded` =
      paste(Spatial_comparisons, Temporal_comparisons, Both_comparisons, `Not recorded_comparisons`, sep = " / ")
  )

## Same studies as Table 1
stopifnot(identical(design_table$`Number of studies`, table1$`Number of studies`))

if(design_unrecorded == 0){
  design_table <- design_table %>%
    select(-`Not recorded`, -`Comparisons: spatial / temporal / both / not recorded`)
} else {
  warning(design_unrecorded, " comparisons have no recorded spatial/temporal design; Table S2 keeps a Not recorded column")
  design_table <- design_table %>% select(-`Comparisons: spatial / temporal / both`)
}

write_table(design_table, "Table_S2_study_design")

#################################################
## Posterior quantile tables
##
## One for each main-text figure and for the taxa:
## main effects (Figure 2), biome (Figure 3A),
## forestation method (Figure 3B), climate zone
## (Figure 3C and Supplementary Figure 17), area and
## age slopes (Figures 4-5). The biome and zone tables
## list the nine interventions and the pooled
## forestation data.
#################################################

write_table(make_ci_table(main_effects, "Intervention"), "Table_S3_main_effects")

write_table(taxa_effects %>% rename(taxa = level) %>% make_ci_table(c("Intervention", "taxa")),
            "Table_S6_taxonomic_group")

write_table(zone_effects %>% rename(Zone = level) %>% make_ci_table(c("Intervention", "Zone")),
            "Table_S9_climate_zone")

write_table(biome_effects %>% rename(Biome = level) %>% make_ci_table(c("Intervention", "Biome")),
            "Table_S7_biome")

if(nrow(forestation_method_effects) > 0){
  write_table(forestation_method_effects %>% rename(Method = method, Biome = level) %>%
                make_ci_table(c("Method", "Biome")),
              "Table_S8_forestation_method")
}

write_table(all_model_coefficients %>%
              filter(model == "area_model", is_moderator_coefficient(parameter)) %>%
              make_ci_table("Intervention"),
            "Table_S10_area_slopes")

write_table(all_model_coefficients %>%
              filter(model == "time_model", is_moderator_coefficient(parameter)) %>%
              make_ci_table("Intervention"),
            "Table_S11_age_slopes")

#################################################
## Source of the species richness data (quoted in
## the Methods; written to output/summaries/)
##
## S_source, for the observations in the main
## models (S_aggregation mean or aggregated, and an
## effect size), so "All" matches Figure 2's N.
## n_studies counts studies with at least one
## observation from that source.
#################################################

data_source_labels <- c(
  "text or table numbers" = "Text or table",
  "figure extraction" = "Figure",
  "author supplied" = "Author provided",
  "calculated from raw data" = "Calculated from raw data"
)

data_source_obs <- bind_rows(lapply(interventions, function(iv){
  d <- read_effects(iv)
  d <- d[modelled_rows(d) & !is.na(d$lrr_richness), ]
  data.frame(
    Intervention = plot_labels[[iv]],
    study_id = d$study_id,
    source = unname(ifelse(is.na(d$S_source) | d$S_source == "", "Not recorded",
                           data_source_labels[d$S_source])),
    stringsAsFactors = FALSE
  )
}))

## Any other S_source value would be NA here
stopifnot(!anyNA(data_source_obs$source))

summarise_data_source <- function(df){
  df %>%
    group_by(source) %>%
    summarise(n_observations = n(), n_studies = n_distinct(study_id), .groups = "drop") %>%
    mutate(percent_observations = round(100 * n_observations / sum(n_observations), 1))
}

write_table(
  bind_rows(
    data_source_obs %>% group_by(Intervention) %>% group_modify(~ summarise_data_source(.x)) %>% ungroup(),
    summarise_data_source(data_source_obs) %>% mutate(Intervention = "All")
  ) %>%
    mutate(
      Intervention = factor(Intervention, levels = c(intervention_display_order, "All")),
      source = factor(source, levels = c(unname(data_source_labels), "Not recorded"))
    ) %>%
    arrange(Intervention, source) %>%
    select(Intervention, source, n_observations, percent_observations, n_studies),
  "data_source_summary", summary_dir
)

#################################################
## Prior schemes (main model and prior sensitivity)
##
## Read from the prior objects in R/models.R. SD
## priors are half-distributions (bounded at 0).
#################################################

describe_prior <- function(priors, cls){
  p <- priors$prior[priors$class == cls & priors$coef == "" & priors$group == ""]
  if(!length(p)) return("flat (brms default)")
  p <- sub("^cauchy", "Cauchy", sub("^normal", "Normal", p))
  if(cls == "sd") paste0("half-", p) else p
}

sigma_prior <- "half-Student-t(3, 0, 2.5) (brms default)"

write_table(
  data.frame(
    Parameter = c(
      "Intercept (overall log response ratio)",
      "Fixed effects (b; e.g. S_aggregation)",
      "Between-study SD (study_id random intercept)",
      "Residual SD (sigma)"
    ),
    Baseline = c(sapply(c("Intercept", "b", "sd"), describe_prior, priors = priors_moderator), sigma_prior),
    Tight = c(sapply(c("Intercept", "b", "sd"), describe_prior, priors = priors_tight), sigma_prior),
    Wide = c(sapply(c("Intercept", "b", "sd"), describe_prior, priors = priors_wide), sigma_prior),
    row.names = NULL
  ),
  "Table_S14_priors"
)

#################################################
## Forestation classifications
#################################################

write_table(
  data.frame(
    Category = c(
      "Author-defined intervention",
      "Biome type",
      "Forestation method",
      "Reference habitat structure"
    ),
    `Option 1` = c(
      paste(
        "Reforestation: active or passive restoration (natural regeneration)",
        "of forest on previously forested land, as described by the study",
        "authors. Excludes forest recovery after disturbances such as",
        "windthrow or wildfire."
      ),
      paste(
        "Forested: the site's coordinates fall in a biome historically",
        "dominated by forest in the RESOLVE Ecoregions 2017 map (Dinerstein",
        "et al. 2017): tropical and subtropical moist broadleaf, dry",
        "broadleaf, and coniferous forests; temperate broadleaf and mixed",
        "forests; temperate conifer forests; boreal forests/taiga; mangroves."
      ),
      paste(
        "Natural forest regrowth: forest established through passive",
        "regeneration or through planting of native species without",
        "subsequent harvest (coded as natural forest restoration)."
      ),
      paste(
        "Woody: the reference plot itself was dominated by trees (e.g.",
        "remnant or secondary forest, woodland, tree plantation), judged per",
        "observation from the study's site description."
      )
    ),
    `Option 2` = c(
      paste(
        "Afforestation: active planting of trees where the available evidence",
        "suggests trees did not grow in recent history (e.g. on grassland,",
        "agricultural fields or pasture, peatland or wetland, mine tailings),",
        "as described by the study authors."
      ),
      paste(
        "Open: the site's coordinates fall in any other biome: tropical and",
        "subtropical grasslands, savannas and shrublands; temperate grasslands,",
        "savannas and shrublands; flooded grasslands and savannas; montane",
        "grasslands and shrublands; tundra; Mediterranean forests, woodlands",
        "and scrub; deserts and xeric shrublands."
      ),
      paste(
        "Plantation: trees planted and managed for timber production",
        "(coded as timber plantation). Other methods (riparian corridors,",
        "windbreaks, fruit/nut plantations) are excluded from this comparison."
      ),
      paste(
        "Open: the reference plot was not tree-dominated (e.g. grassland,",
        "pasture, cropland, old field, shrubland, bare or mined land), judged",
        "per observation from the study's site description."
      )
    ),
    check.names = FALSE
  ),
  "Table_S5_habitat_categories"
)

#################################################
## Papers by language at each screening stage
## (quoted in the Methods; written to
## output/summaries/)
##
##   Title/abstract: every record screened, with
##     Web of Science's language field
##     (data/title_abstract_language.csv); a
##     record found by two searches counts once in
##     "All".
##   Full text: every study assessed at full text
##     (data/studies.csv), with the language read
##     from the paper itself
##     (data/paper_language.csv).
##   Included: studies in the main models.
## Re- and afforestation were screened together and
## are reported as Forestation. A study under several
## interventions counts once in "All".
##
## papers_by_language: every language, per stage.
## papers_english_vs_other: English vs all other
##   languages, per stage. Papers whose language is
##   not known (not recorded by Web of Science, no
##   file found, or not detectable) are counted
##   separately and left out of the percentages.
#################################################

if(file.exists(paper_language_file) && file.exists(title_abstract_language_file)){

  paper_language <-
    read.csv(paper_language_file, stringsAsFactors = FALSE) %>%
    mutate(key = tolower(study_id),
           language = ifelse(text_source == "no file found", "No file found", language)) %>%
    distinct(key, .keep_all = TRUE) %>%
    select(key, language)

  language_group_labels <- c(
    "[Re-af]forestation" = "Forestation",
    "Reforestation" = "Forestation",
    "Afforestation" = "Forestation",
    "BECCS (bioenergy)" = "Bioenergy cropping",
    "Wetland restoration - inland" = "Inland wetland restoration",
    "Wind energy" = "Wind"
  )

  relabel_group <- function(x){
    ifelse(x %in% names(language_group_labels), language_group_labels[x], x)
  }

  language_group_order <- c(
    "All", "Forestation", "Grassland restoration", "Inland wetland restoration",
    "Coastal wetland restoration", "Agroforestry", "Bioenergy cropping", "Wind", "Solar"
  )

  unknown_languages <- c("Not recorded", "Unknown", "No file found", "Not checked")

  title_abstract_records <-
    read.csv(title_abstract_language_file, stringsAsFactors = FALSE) %>%
    transmute(group = relabel_group(intervention), key = record_id, language)

  search_groups <- c(
    "Reforestation and afforestation" = "Forestation",
    "Grassland restoration" = "Grassland restoration",
    "Inland wetland restoration" = "Inland wetland restoration",
    "Coastal wetland restoration" = "Coastal wetland restoration",
    "Agroforestry" = "Agroforestry",
    "Bioenergy cropping" = "Bioenergy cropping",
    "Wind (2021 search)" = "Wind",
    "Wind (2026 search)" = "Wind",
    "Solar (2021 search)" = "Solar",
    "Solar (2026 search)" = "Solar"
  )

  screened_studies <-
    read.csv(studies_file, stringsAsFactors = FALSE) %>%
    transmute(group = unname(search_groups[search]), key = tolower(study_id)) %>%
    filter(!is.na(group)) %>%
    distinct()

  included_studies <-
    data_source_obs %>%
    transmute(group = relabel_group(Intervention), key = tolower(study_id)) %>%
    distinct()

  not_screened <- setdiff(included_studies$key, screened_studies$key)
  if(length(not_screened) > 0){
    log_message(length(not_screened), " modelled studies are not in studies.csv: ",
                paste(head(not_screened, 5), collapse = ", "))
  }

  add_language <- function(df){
    df %>%
      left_join(paper_language, by = "key") %>%
      mutate(language = ifelse(is.na(language), "Not checked", language))
  }

  count_language <- function(df, n_name){
    bind_rows(
      df %>% distinct(key, language) %>% mutate(group = "All"),
      df %>% distinct(group, key, language)
    ) %>%
      count(group, language, name = n_name)
  }

  language_counts <-
    count_language(title_abstract_records, "n_title_abstract") %>%
    full_join(count_language(add_language(screened_studies), "n_full_text"), by = c("group", "language")) %>%
    full_join(count_language(add_language(included_studies), "n_included"), by = c("group", "language")) %>%
    mutate(across(starts_with("n_"), ~ ifelse(is.na(.x), 0L, .x)),
           group = factor(group, levels = language_group_order))

  write_table(
    language_counts %>%
      group_by(group) %>%
      mutate(
        percent_title_abstract = round(100 * n_title_abstract / sum(n_title_abstract), 1),
        percent_full_text = round(100 * n_full_text / sum(n_full_text), 1),
        percent_included = round(100 * n_included / sum(n_included), 1),
        percent_of_full_text_included = ifelse(n_full_text > 0, round(100 * n_included / n_full_text, 1), NA)
      ) %>%
      ungroup() %>%
      arrange(group, desc(n_title_abstract), desc(n_full_text)) %>%
      select(Intervention = group, language, n_title_abstract, percent_title_abstract,
             n_full_text, percent_full_text, n_included, percent_included, percent_of_full_text_included),
    "papers_by_language", summary_dir
  )

  write_table(
    language_counts %>%
      pivot_longer(starts_with("n_"), names_to = "stage", values_to = "n") %>%
      mutate(stage = factor(stage, levels = c("n_title_abstract", "n_full_text", "n_included"),
                            labels = c("Title/abstract", "Full text", "Included")),
             class = ifelse(language == "English", "english",
                            ifelse(language %in% unknown_languages, "unknown", "other"))) %>%
      group_by(group, stage) %>%
      summarise(n_total = sum(n),
                n_english = sum(n[class == "english"]),
                n_other_language = sum(n[class == "other"]),
                n_language_unknown = sum(n[class == "unknown"]),
                other_languages = paste(language[class == "other" & n > 0][order(-n[class == "other" & n > 0])],
                                        collapse = ", "),
                .groups = "drop") %>%
      mutate(percent_english = ifelse(n_english + n_other_language > 0,
                                      round(100 * n_english / (n_english + n_other_language), 1), NA),
             percent_other_language = ifelse(n_english + n_other_language > 0,
                                             round(100 * n_other_language / (n_english + n_other_language), 1), NA)) %>%
      arrange(group, stage) %>%
      select(Intervention = group, stage, n_total, n_english, n_other_language, n_language_unknown,
             percent_english, percent_other_language, other_languages),
    "papers_english_vs_other", summary_dir
  )

}

#################################################
## Between-study heterogeneity of the main models
##
## Per intervention: tau (between-study SD, log
## scale) with its 95% CI and the share of
## study-level effects with the same sign as the
## main effect.
#################################################

study_heterogeneity <-
  read.csv(file.path(summary_dir, "all_study_heterogeneity.csv"), stringsAsFactors = FALSE)

tau_summary <- read.csv(file.path(summary_dir, "all_tau.csv"), stringsAsFactors = FALSE)

ci_text <- function(lo, hi, digits = 2) paste0(round(lo, digits), " - ", round(hi, digits))

same_sign <-
  study_heterogeneity %>%
  filter(row_type == "Study") %>%
  left_join(study_heterogeneity %>% filter(row_type == "Overall") %>% select(intervention, fixed_median = median),
            by = "intervention") %>%
  group_by(intervention) %>%
  summarise(pct_same_direction_as_fixed = round(100 * mean(sign(median) == sign(fixed_median)), 1),
            .groups = "drop")

write_table(
  tau_summary %>%
    left_join(same_sign, by = "intervention") %>%
    mutate(
      Intervention = plot_labels[intervention],
      tau = round(tau_median, 3),
      tau_ci_95 = ci_text(tau_q025, tau_q975)
    ) %>%
    select(Intervention, n_studies, n_observations, pct_same_direction_as_fixed, tau, tau_ci_95) %>%
    arrange(Intervention),
  "Table_S4_heterogeneity"
)

leave_one_out <- read.csv(file.path(summary_dir, "all_leave_one_out.csv"), stringsAsFactors = FALSE)

#################################################
## Study-level heterogeneity and leave-one-out
## figures
##
## One panel per intervention, x in log2(fold change)
## on a fixed symmetric axis. Heterogeneity: the
## distribution of the study-level effects.
## Leave-one-out: the distribution of the main-effect
## estimate across refits. Both show the main-model
## estimate (solid line) with its 95% CI (shaded) and
## no effect (dashed).
#################################################

make_combined_histogram_figure <- function(hist_df, overall_df, y_label, filename, x_breaks, binwidth_log2){

  hist_df$Intervention <- order_interventions(plot_labels[hist_df$intervention])
  overall_df$Intervention <- order_interventions(plot_labels[overall_df$intervention])

  p <-
    ggplot(hist_df, aes(x = log2_fold_change)) +
    geom_rect(data = overall_df,
              aes(xmin = q025_log2, xmax = q975_log2, ymin = -Inf, ymax = Inf, fill = Intervention),
              inherit.aes = FALSE, alpha = 0.12) +
    ## no outline, so empty bins draw no baseline
    geom_histogram(aes(y = after_stat(density), fill = Intervention),
                   binwidth = binwidth_log2, colour = NA, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(data = overall_df, aes(xintercept = median_log2, colour = Intervention), linewidth = 1) +
    facet_wrap(~ Intervention, ncol = 3, scales = "free_y", labeller = n_strip_labeller("main_model")) +
    scale_x_continuous(limits = log2(range(x_breaks)), breaks = log2(x_breaks),
                       labels = fold_fraction_labels(x_breaks), oob = scales::oob_keep) +
    scale_fill_manual(values = intervention_colors) +
    scale_colour_manual(values = intervention_colors) +
    theme_lbms() +
    theme(legend.position = "none") +
    labs(x = fold_change_label, y = y_label)

  save_lbms_plot(p, file.path(figure_dir, filename), width = 12, height = 10)

}

make_combined_histogram_figure(
  study_heterogeneity %>%
    filter(row_type == "Study") %>%
    mutate(log2_fold_change = median / log(2)),
  study_heterogeneity %>%
    filter(row_type == "Overall") %>%
    mutate(median_log2 = median / log(2), q025_log2 = q025 / log(2), q975_log2 = q975 / log(2)),
  "Density of studies",
  "Figure_S13_study_heterogeneity.png",
  x_breaks = c(1/4, 1/2, 1, 2, 4),
  binwidth_log2 = 0.2
)

make_combined_histogram_figure(
  leave_one_out %>% mutate(log2_fold_change = log2(fold_change)),
  main_effects %>%
    mutate(median_log2 = log(fold_change) / log(2), q025_log2 = q025 / log(2), q975_log2 = q975 / log(2)),
  "Density of leave-one-out refits",
  "Figure_S27_leave_one_study_out.png",
  x_breaks = c(1/2, 2/3, 1, 3/2, 2),
  binwidth_log2 = 0.03
)

log_message("Figures and tables complete.")
