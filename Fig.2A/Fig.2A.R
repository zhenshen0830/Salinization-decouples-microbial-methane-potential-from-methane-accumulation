# =========================
# Two-way ANOVA
# mcrA
# =========================

rm(list = ls())


library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)

dat <- dat %>%
  mutate(
    Lake = factor(Lake,
                  levels = c("AFresh", "Brackish", "Saline"),
                  labels = c("Fresh", "Brackish", "Saline")),
    Size = factor(Size,
                  levels = c("A64", "B5", "C0.2"),
                  labels = c(">64", "5-64", "0.2-5"))
  )
pal_lake <- c(
  "Fresh"    = "#4cb8d0",
  "Brackish" = "#2b8cbe",
  "Saline"   = "#00558f"
)

pal_size <- c(
  ">64"   = "#83cbac",
  "5-64"  = "#756bb1",
  "0.2-5" = "#fc806f"
)

shape_lake <- c(
  "Fresh"    = 16,
  "Brackish" = 17,
  "Saline"   = 15
)

shape_size <- c(
  ">64"   = 16,
  "5-64"  = 17,
  "0.2-5" = 15
)

fit <- aov(mcrA ~ Lake * Size, data = dat)
anova_tab <- summary(fit)[[1]]
summary(fit)

sum_df <- dat %>%
  group_by(Lake, Size) %>%
  summarise(
    mean = mean(mcrA, na.rm = TRUE),
    sd   = sd(mcrA, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n),
    .groups = "drop"
  )

size_stats <- dat %>%
  group_by(Size) %>%
  summarise(
    med  = median(mcrA, na.rm = TRUE),
    mean = mean(mcrA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Size)

size_stats$x <- seq_len(nrow(size_stats))

lake_stats <- dat %>%
  group_by(Lake) %>%
  summarise(
    med  = median(mcrA, na.rm = TRUE),
    mean = mean(mcrA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Lake)

lake_stats$x <- seq_len(nrow(lake_stats))

y_min  <- min(dat$mcrA, na.rm = TRUE)
y_max  <- max(dat$mcrA, na.rm = TRUE)
y_rng  <- y_max - y_min

y_top  <- y_max + y_rng * 0.12
y_text <- y_max + y_rng * 0.055

sum_df_p1 <- sum_df %>% mutate(panel = "mcrA ~ Lake | Size")
sum_df_p4 <- sum_df %>% mutate(panel = "mcrA ~ Size | Lake")
dat_p2    <- dat %>% mutate(panel = "mcrA ~ Size | Size")
dat_p3    <- dat %>% mutate(panel = "mcrA ~ Lake | Lake")

theme_my <- theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    axis.title = element_text(face = "bold", color = "black", size = 11),
    axis.text = element_text(color = "black", size = 10),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    
    strip.background = element_rect(fill = "#BFBFBF", color = "#BFBFBF", linewidth = 0.5),
    strip.text = element_text(size = 10.5, face = "plain", color = "black",
                              margin = margin(3, 0, 3, 0)),
    legend.position = "top",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(3, "mm"),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5)
  )

p1 <- ggplot(sum_df_p1, aes(x = Lake, y = mean, group = Size, color = Size, shape = Size)) +
  geom_line(linewidth = 0.42) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.1, linewidth = 0.42) +
  scale_color_manual(
    values = pal_size,
    guide = guide_legend(order = 2, nrow = 1, byrow = TRUE)
  ) +
  scale_shape_manual(
    values = shape_size,
    guide = guide_legend(order = 2, nrow = 1, byrow = TRUE)
  ) +
  facet_wrap(~panel, nrow = 1) +
  labs(x = NULL, y = "mcrA") +
  coord_cartesian(ylim = c(y_min, y_top), clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.18))) +
  theme_my

p2 <- ggplot(dat_p2, aes(x = Size, y = mcrA, fill = Size)) +
  stat_boxplot(geom = "errorbar", width = 0.13, linewidth = 0.38, color = "black") +
  geom_boxplot(
    width = 0.54, alpha = 0.9,
    outlier.shape = NA, color = "black", linewidth = 0.45,
    fatten = 0
  ) +
  geom_jitter(aes(fill = Size),
              shape = 21, size = 0.95, width = 0.09,
              alpha = 0.42, stroke = 0.2, color = "gray35") +
  geom_segment(data = size_stats,
               aes(x = x - 0.27, xend = x + 0.27, y = med, yend = med),
               inherit.aes = FALSE, color = "white", linewidth = 0.72) +
  geom_point(data = size_stats,
             aes(x = x, y = mean),
             inherit.aes = FALSE,
             shape = 21, size = 1.85, fill = "white", color = "white", stroke = 0.2) +
  scale_fill_manual(values = pal_size, guide = "none") +
  facet_wrap(~panel, nrow = 1) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim = c(y_min, y_top), clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.18))) +
  theme_my

p3 <- ggplot(dat_p3, aes(x = Lake, y = mcrA, fill = Lake)) +
  stat_boxplot(geom = "errorbar", width = 0.13, linewidth = 0.38, color = "black") +
  geom_boxplot(
    width = 0.54, alpha = 0.9,
    outlier.shape = NA, color = "black", linewidth = 0.45,
    fatten = 0
  ) +
  geom_jitter(aes(fill = Lake),
              shape = 21, size = 0.95, width = 0.09,
              alpha = 0.42, stroke = 0.2, color = "gray35") +
  geom_segment(data = lake_stats,
               aes(x = x - 0.27, xend = x + 0.27, y = med, yend = med),
               inherit.aes = FALSE, color = "white", linewidth = 0.72) +
  geom_point(data = lake_stats,
             aes(x = x, y = mean),
             inherit.aes = FALSE,
             shape = 21, size = 1.85, fill = "white", color = "white", stroke = 0.2) +
  scale_fill_manual(values = pal_lake, guide = "none") +
  facet_wrap(~panel, nrow = 1) +
  
  labs(x = "Lake", y = "mcrA") +
  coord_cartesian(ylim = c(y_min, y_top), clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.18))) +
  theme_my

p4 <- ggplot(sum_df_p4, aes(x = Size, y = mean, group = Lake, color = Lake, shape = Lake)) +
  geom_line(linewidth = 0.42) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.1, linewidth = 0.42) +
  scale_color_manual(
    values = pal_lake,
    guide = guide_legend(order = 1, nrow = 1, byrow = TRUE)
  ) +
  scale_shape_manual(
    values = shape_lake,
    guide = guide_legend(order = 1, nrow = 1, byrow = TRUE)
  ) +
  facet_wrap(~panel, nrow = 1) +
  labs(x = "Size (μm)", y = NULL) +
  coord_cartesian(ylim = c(y_min, y_top), clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.18))) +
  theme_my

final_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(guides = "collect", widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(
    title = "mcrA: main effects and two-way interaction",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                                margin = margin(b = 3)),
      plot.margin = margin(1, 1, 1, 1)
    )
  ) &
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box = "horizontal"
  )


print(final_plot)

