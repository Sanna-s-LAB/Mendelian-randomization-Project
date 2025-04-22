library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)

# Upload useful data 
lifelines <- read_excel("~/Results tables/Table_LIFELINES_after_leaveoneout.xlsx", sheet = 1)
lifelines <- subset(lifelines, select=c("exposure","outcome","Consistency_GLGC_Lifelines"))
lifelines <- lifelines %>% distinct()
significant <- read_excel("~/Results tables/Table_of_significant_main_for_forest_plots.xlsx", sheet = 1)


significant <- merge(significant, lifelines, all.x = TRUE)
count_associations <- significant %>%
  dplyr::select("exposure", "outcome") %>%
  distinct()

significant$Cochran[is.na(significant$Cochran)] <- "no"
significant$Consistency_GLGC_Lifelines[is.na(significant$Consistency_GLGC_Lifelines)] <- "no"

# Trasnform dataframe in long format
significant_long <- significant %>%
  pivot_longer(
    cols = c(b.men, b.women, se.men, se.women, fdr.men, fdr.women),
    names_to = c(".value", "Sex"),
    names_pattern = "(.*)\\.(men|women)"
  ) %>%
  mutate(
    Sex = recode(Sex, men = "Men", women = "Women"),
    ci_lower = b - 1.96 * se,
    ci_upper = b + 1.96 * se,
    sig = case_when(
      (Sex == "Men" & (SignificantBysex_CISTRANS == 2 | SignificantBysex_CISTRANS == 3| SignificantBysex_CISTRANS == 4)) ~ "Significant",
      (Sex == "Women" & (SignificantBysex_CISTRANS == 1 | SignificantBysex_CISTRANS == 3| SignificantBysex_CISTRANS == 4)) ~ "Significant",
      is.na(fdr) ~ "NA",
      TRUE ~ "Not significant"
    ),
    ExposureOutcome = paste(Assay),
    highlight_GLGC = if_else(Consistency_GLGC_Lifelines == "yes", TRUE, FALSE),
    symbol = case_when(
      Sex == "Men" & Consistency_GLGC_Lifelines == "yes" & SignificantBysex_CISTRANS == 2 ~ "*",
      Sex == "Women" & Consistency_GLGC_Lifelines == "yes" & SignificantBysex_CISTRANS == 1 ~ "*",
      TRUE ~ ""
    ),
    highlight_Cochran = if_else(Cochran == "yes", TRUE, FALSE),
    symbol1 = case_when(
      Sex == "Men" & Cochran == "yes" ~ "†",
      Sex == "Women" & Cochran == "yes" ~ "†",
      TRUE ~ ""
    ),
    color_group = case_when(
      Sex == "Women" & (SignificantBysex_CISTRANS == 1 | SignificantBysex_CISTRANS == 3| SignificantBysex_CISTRANS == 4) ~ "Women_significant",
      Sex == "Men" & (SignificantBysex_CISTRANS == 2 | SignificantBysex_CISTRANS == 3| SignificantBysex_CISTRANS == 4) ~ "Men_significant",
      TRUE ~ "Not_significant"
    ),
    color_group = factor(color_group, levels = c("Women_significant", "Men_significant", "Not_significant"))
  ) %>%
  group_by(outcome) %>%
  mutate(ExposureOutcome = factor(ExposureOutcome, levels = unique(ExposureOutcome[order(b)]))) %>%
  ungroup()

# Create plots
plots <- significant_long %>%
  group_split(outcome) %>%
  lapply(function(significant_sub) {
    
    min_x <- min(significant_sub$ci_lower, na.rm = TRUE) - 0.05 * diff(range(significant_sub$ci_lower, na.rm = TRUE))
    
    ggplot(significant_sub, aes(x = b, y = ExposureOutcome, color = color_group)) +
      geom_text(aes(x = min_x, label = symbol), 
                color = "black", size = 7, vjust = 0.72) +
      geom_text(aes(x = min_x, label = symbol1), 
                color = "black", size = 5, vjust = 0.52, hjust=2) +
      geom_point(position = position_dodge(width = 0.6), size = 2.5) +
      geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                     position = position_dodge(width = 0.6),
                     height = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      scale_color_manual(
        values = c(
          "Women_significant" = "#e41a1c",  # Red
          "Men_significant" = "#377eb8",    # Blue
          "Not_significant" = "gray70"      # Grey
        ),
        labels = c(
          "Women_significant" = "Women (Significant)",
          "Men_significant" = "Men (Significant)",
          "Not_significant" = "Not significant"
        )
      ) +
      labs(
        title = unique(significant_sub$outcome),
        x = "Effect (IVW Beta ± 95% CI)",
        y = "Exposure",
        color = ""
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.y = element_text(size = 8),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_markdown(
          hjust = 0.5,
          size = 13,
          margin = margin(t = 10)
        )
      ) +
      guides(color = guide_legend(override.aes = list(size = 3)))
  })


plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
plots[[5]]

for (i in 1:5) {
  ggsave(
    filename = paste0("~/Results tables/plot_outcome_WMtogether", i, ".jpeg"),
    plot = plots[[i]],
    device = "jpeg",
    width = 7,
    height = 11,
    units = "in",
    dpi = 300
  )
}



