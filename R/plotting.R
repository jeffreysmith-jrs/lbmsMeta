#################################################
## plotting.R
##
## Theme, colours, fold-change axes and figure
## saving for step 06.
#################################################

#################################################
## Theme and saving
#################################################

theme_lbms <- function(base_size = 14){

  theme_classic(base_size = base_size) +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = rel(0.95)),
      legend.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )

}

save_lbms_plot <- function(plot_object, filename, width = 7, height = 5, dpi = 600){

  ggsave(filename = filename, plot = plot_object, width = width, height = height,
         dpi = dpi, units = "in", bg = "white")

}

#################################################
## Colours
#################################################

intervention_colors <- c(
  "Reforestation" = "#117733",
  "Grassland restoration" = "#999933",
  "Inland wetland restoration" = "#88CCEE",
  "Coastal wetland restoration" = "#44AA99",
  "Agroforestry" = "#DDCC77",
  "Afforestation" = "#332288",
  "Bioenergy cropping" = "#882255",
  "Solar" = "#CC6677",
  "Wind" = "#AA4499"
)

scale_color_lbms <- function() scale_color_manual(values = intervention_colors)

scale_fill_lbms <- function() scale_fill_manual(values = intervention_colors)

taxa_colors <- c(
  "Plants" = "#2E8B57",
  "Invertebrates" = "#D95F02",
  "Vertebrates" = "#1F78B4",
  "Microbes" = "#6A3D9A"
)

biome_colors <- c("Forested" = "#2D6A4F", "Open" = "#C08552")

zone_colors <- c("Temperate" = "#F28E2B", "Tropical" = "#7B68EE")

structure_colors <- c("Woody" = "#2D6A4F", "Open" = "#C08552")

human_use_colors <- c(
  "Intact natural" = "#1B7837",
  "Degraded natural" = "#A6D96A",
  "Human use" = "#8C510A"
)

experiment_type_colors <- c("Natural" = "#2A6F97", "Manipulative" = "#BC4749")

## Sensitivity-analysis comparisons
aggregation_colors <- c("Aggregated" = "#009E73", "Mean" = "#D55E00", "Pooled" = "#0072B2")

weighting_colors <- c("Weighted" = "#D55E00", "Unweighted" = "#0072B2")

publication_bias_colors <- c(
  "Variance-weighted" = "#56B4E9",
  "Main model" = "#0072B2",
  "PEESE" = "#CC79A7"
)

high_n_colors <- c("Full sample" = "#0072B2", "Excl. top 10% studies" = "#D55E00")

one_per_study_colors <- c("Full sample" = "#0072B2", "One effect per study" = "#009E73")

prior_sensitivity_colors <- c(
  "Baseline prior" = "#0072B2",
  "Tight prior" = "#D55E00",
  "Wide prior" = "#009E73"
)

## Moderator-slope figures: strong support (95% CI
## excludes zero) vs not
slope_support_values <- c("Significant" = "dodgerblue4", "Not significant" = "grey70")

#################################################
## Fold-change axes
##
## Every fold-change axis is on a log2 scale,
## symmetric about 1 (1/L to L), with reciprocal
## ticks labelled as fractions (3/4 and 4/3, 1/2 and
## 2, ...). L is the smallest of fold_axis_steps that
## covers everything drawn on the axis, worked out
## per panel, and never more than 4. Pass L to fix
## the range instead. Values beyond the range are
## kept (histograms still bin them) and clipped at
## the panel edge.
#################################################

fold_axis_steps <- c(4/3, 3/2, 2, 3, 4)

fold_axis_breaks <- function(L){
  upper <- switch(
    as.character(match(TRUE, abs(fold_axis_steps - L) < 1e-9)),
    "1" = c(4/3),
    "2" = c(3/2),
    "3" = c(4/3, 5/3, 2),
    "4" = c(5/3, 2, 3),
    "5" = c(2, 4)
  )
  c(rev(1 / upper), 1, upper)
}

## 0.75 -> "3/4", 4/3 -> "4/3", 2 -> "2", 1/8 -> "1/8"
fold_fraction_labels <- function(x){
  vapply(x, function(v){
    if(is.na(v)) return(NA_character_)
    if(v < 1 && abs(1 / v - round(1 / v)) < 1e-6){
      return(paste0("1/", round(1 / v)))
    }
    for(q in 1:16){
      p <- v * q
      if(abs(p - round(p)) < 1e-6){
        return(if(q == 1) as.character(round(p)) else paste0(round(p), "/", q))
      }
    }
    format(signif(v, 2))
  }, character(1))
}

fold_axis_limit <- function(lims){
  lims <- lims[is.finite(lims) & lims > 0]
  if(!length(lims)) return(2)
  m <- max(abs(log2(lims)))
  hit <- fold_axis_steps[fold_axis_steps >= 2^m - 1e-9]
  if(length(hit)) hit[1] else max(fold_axis_steps)
}

## The breaks function sees the limits after ggplot's
## padding, so it rounds down to the step in use.
fold_axis_limit_down <- function(lims){
  m <- max(abs(log2(lims[is.finite(lims) & lims > 0])))
  hit <- fold_axis_steps[fold_axis_steps <= 2^m + 1e-9]
  if(length(hit)) max(hit) else min(fold_axis_steps)
}

scale_fold <- function(axis = "x", L = NULL, ...){
  scale_fun <- if(axis == "x") scale_x_continuous else scale_y_continuous
  scale_fun(
    trans = "log2",
    limits =
      if(is.null(L)) function(lims){ L <- fold_axis_limit(lims); c(1 / L, L) }
      else c(1 / L, L),
    breaks = function(lims) fold_axis_breaks(fold_axis_limit_down(lims)),
    labels = fold_fraction_labels,
    oob = scales::oob_keep,
    ...
  )
}
