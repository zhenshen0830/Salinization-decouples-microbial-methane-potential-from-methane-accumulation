# ============================================================
# SEM/path model for methane concentration
# ============================================================

library(lavaan)
library(semPlot)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ape)
library(patchwork)
library(cowplot)
library(magick)

rm(list = ls())


# ============================================================
# 1. Read data
# ============================================================

dat0 <- read.csv("Fig.6DEF.csv", fileEncoding = "GBK", check.names = FALSE)

dat <- dat0 %>%
  rename(
    salinity     = `salinity（‰）`,
    Ln_salinity  = `Ln salinity (‰)`,
    LnDOC        = `Ln DOC (mg/L)`,
    TN           = `TN (mg/L)`,
    TP           = `TP (mg/L）`,
    WT           = `WT (℃)`,
    C1           = `C1 (R.U.)`,
    LnmcrA       = LnmcrA,
    LnphnJ       = LnphnJ,
    LnCCH4       = `LnCCH4 (nmol/L)`,
    LnFCH4       = `LnFCH4 (umol/m2/d)`
  ) %>%
  select(
    site, Group,
    Ln_salinity, WT,
    TN, TP, NH4,
    LnDOC, HIX, C1,
    LnmcrA, LnphnJ,
    LnCCH4, LnFCH4
  ) %>%
  drop_na()

# ============================================================
# 2. Construct PC1 variables using PCA
# ============================================================

make_pc1 <- function(df, cols) {
  p <- prcomp(
    scale(df[, cols]),
    center = FALSE,
    scale. = FALSE
  )
  
  scores <- p$x[, 1]
  loadings <- p$rotation[, 1]
  
  # 
  if (sum(loadings) < 0) {
    scores <- -scores
    loadings <- -loadings
  }
  
  list(
    scores = as.numeric(scores),
    loadings = loadings,
    var = summary(p)$importance[2, 1]
  )
}

nut <- make_pc1(dat, c("TN", "TP", "NH4"))
org <- make_pc1(dat, c("LnDOC", "HIX", "C1"))
gen <- make_pc1(dat, c("LnmcrA", "LnphnJ"))

dat$Nutrient_PC1 <- nut$scores
dat$Organic_PC1  <- org$scores
dat$Gene_PC1     <- gen$scores

# 
pc1_loadings <- data.frame(
  Module = c(
    rep("Nutrient_PC1", length(nut$loadings)),
    rep("Organic_PC1", length(org$loadings)),
    rep("Gene_PC1", length(gen$loadings))
  ),
  Variable = c(
    names(nut$loadings),
    names(org$loadings),
    names(gen$loadings)
  ),
  Loading = c(
    as.numeric(nut$loadings),
    as.numeric(org$loadings),
    as.numeric(gen$loadings)
  )
)

pc1_variance <- data.frame(
  Module = c("Nutrient_PC1", "Organic_PC1", "Gene_PC1"),
  Explained_variance_PC1 = c(nut$var, org$var, gen$var)
)

write.csv(pc1_loadings, "PC1_loadings.csv", row.names = FALSE)
write.csv(pc1_variance, "PC1_explained_variance.csv", row.names = FALSE)


# ============================================================
# 3. Z score transformation
# ============================================================

zdat <- dat %>%
  mutate(across(
    c(
      Ln_salinity, WT,
      Nutrient_PC1, Organic_PC1, Gene_PC1,
      LnCCH4, LnFCH4
    ),
    ~ as.numeric(scale(.x)),
    .names = "z_{.col}"
  ))

write.csv(zdat, "SEM_zscore_plotdata.csv", row.names = FALSE)

# ============================================================
# 4. Updated lavaan SEM model
# ============================================================

model_sem <- '
  # Nutrient module
  z_Nutrient_PC1 ~ a1*z_Ln_salinity + a2*z_WT

  # Organic matter module
  z_Organic_PC1 ~ b1*z_WT + b2*z_Nutrient_PC1

  # Gene module
  z_Gene_PC1 ~ c1*z_WT + c2*z_Nutrient_PC1 + c3*z_Organic_PC1

  # Methane concentration
  z_LnCCH4 ~ d1*z_Ln_salinity + d2*z_WT + d3*z_Organic_PC1 + d4*z_Gene_PC1
'

Fit <- lavaan::sem(
  model_sem,
  data = zdat,
  meanstructure = TRUE,
  fixed.x = TRUE
)

summary(
  Fit,
  standardized = TRUE,
  fit.measures = TRUE,
  rsquare = TRUE
)

# ============================================================
# 5. Export SEM fit indices
# ============================================================

fit_out <- fitMeasures(
  Fit,
  c(
    "chisq", "df", "pvalue",
    "cfi", "tli", "nfi", "ifi",
    "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
    "srmr",
    "aic", "bic"
  )
)

fit_out
#write.csv(fit_out, "SEM_model_fit_indices.csv", row.names = TRUE)

# ============================================================
# 6. Export standardized path coefficients
# ============================================================

param_out <- parameterEstimates(
  Fit,
  standardized = TRUE
)

path_out <- param_out %>%
  filter(op == "~") %>%
  select(
    Response = lhs,
    Predictor = rhs,
    Estimate = est,
    SE = se,
    Z = z,
    P_value = pvalue,
    Std_all = std.all
  ) %>%
  mutate(
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      P_value < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  )

write.csv(
  path_out,
  "SEM_standardized_path_coefficients.csv",
  row.names = FALSE
)

# ============================================================
# 7. Export R-square
# ============================================================

r2_out <- inspect(Fit, "r2")

r2_out <- data.frame(
  Variable = names(r2_out),
  R2 = as.numeric(r2_out)
)

write.csv(r2_out, "SEM_Rsquare.csv", row.names = FALSE)

# ============================================================
# 8. Standardized direct effect on LnCH4
# ============================================================

direct_effects <- path_out %>%
  filter(Response == "z_LnCCH4") %>%
  mutate(
    Predictor_clean = recode(
      Predictor,
      z_Ln_salinity = "Ln salinity",
      z_WT = "WT",
      z_Organic_PC1 = "Organic PC1",
      z_Gene_PC1 = "Gene PC1"
    ),
    Predictor_clean = factor(
      Predictor_clean,
      levels = Predictor_clean[order(Std_all)]
    )
  )

write.csv(
  direct_effects,
  "Standardized_direct_effect_on_LnCH4.csv",
  row.names = FALSE
)

p_direct <- ggplot(
  direct_effects,
  aes(x = Predictor_clean, y = Std_all)
) +
  geom_col(width = 0.65, fill = "grey65", color = "black") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_text(
    aes(label = paste0(round(Std_all, 2), Significance)),
    hjust = ifelse(direct_effects$Std_all > 0, -0.15, 1.15),
    size = 4
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Standardized direct effect on LnCH4"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank()
  )

ggsave(
  "Standardized_direct_effect_on_LnCH4.png",
  p_direct,
  width = 6.5,
  height = 4.2,
  dpi = 600
)

ggsave(
  "Standardized_direct_effect_on_LnCH4.pdf",
  p_direct,
  width = 6.5,
  height = 4.2
)

# ============================================================
# 9. Standardized total effect on LnCH4
# ============================================================
# a1: salinity -> nutrient
# a2: WT -> nutrient
# b1: WT -> organic
# b2: nutrient -> organic
# c1: WT -> gene
# c2: nutrient -> gene
# c3: organic -> gene
# d1: salinity -> CH4
# d2: WT -> CH4
# d3: organic -> CH4
# d4: gene -> CH4

coef_std <- standardizedSolution(Fit) %>%
  filter(op == "~") %>%
  select(lhs, rhs, est.std)

get_beta <- function(lhs_name, rhs_name) {
  x <- coef_std %>%
    filter(lhs == lhs_name, rhs == rhs_name) %>%
    pull(est.std)
  
  if (length(x) == 0) {
    return(0)
  } else {
    return(x)
  }
}

a1 <- get_beta("z_Nutrient_PC1", "z_Ln_salinity")
a2 <- get_beta("z_Nutrient_PC1", "z_WT")

b1 <- get_beta("z_Organic_PC1", "z_WT")
b2 <- get_beta("z_Organic_PC1", "z_Nutrient_PC1")

c1 <- get_beta("z_Gene_PC1", "z_WT")
c2 <- get_beta("z_Gene_PC1", "z_Nutrient_PC1")
c3 <- get_beta("z_Gene_PC1", "z_Organic_PC1")

d1 <- get_beta("z_LnCCH4", "z_Ln_salinity")
d2 <- get_beta("z_LnCCH4", "z_WT")
d3 <- get_beta("z_LnCCH4", "z_Organic_PC1")
d4 <- get_beta("z_LnCCH4", "z_Gene_PC1")

# Total effect formulas
# salinity:
# direct: salinity -> CH4
# indirect:
# salinity -> nutrient -> organic -> CH4
# salinity -> nutrient -> gene -> CH4
# salinity -> nutrient -> organic -> gene -> CH4
total_salinity <- d1 +
  a1 * b2 * d3 +
  a1 * c2 * d4 +
  a1 * b2 * c3 * d4

# WT:
# direct: WT -> CH4
# indirect:
# WT -> organic -> CH4
# WT -> gene -> CH4
# WT -> nutrient -> organic -> CH4
# WT -> nutrient -> gene -> CH4
# WT -> organic -> gene -> CH4
# WT -> nutrient -> organic -> gene -> CH4
total_WT <- d2 +
  b1 * d3 +
  c1 * d4 +
  a2 * b2 * d3 +
  a2 * c2 * d4 +
  b1 * c3 * d4 +
  a2 * b2 * c3 * d4

# Nutrient PC1:
# nutrient -> organic -> CH4
# nutrient -> gene -> CH4
# nutrient -> organic -> gene -> CH4
total_nutrient <- b2 * d3 +
  c2 * d4 +
  b2 * c3 * d4

# Organic PC1:
# direct: organic -> CH4
# indirect: organic -> gene -> CH4
total_organic <- d3 + c3 * d4

# Gene PC1:
# direct only
total_gene <- d4

total_effects <- data.frame(
  Predictor = c(
    "Ln salinity",
    "WT",
    "Nutrient PC1",
    "Organic PC1",
    "Gene PC1"
  ),
  Total_effect = c(
    total_salinity,
    total_WT,
    total_nutrient,
    total_organic,
    total_gene
  )
) %>%
  mutate(
    Predictor = factor(
      Predictor,
      levels = Predictor[order(Total_effect)]
    )
  )

write.csv(
  total_effects,
  "Standardized_total_effect_on_LnCH4.csv",
  row.names = FALSE
)

p_total <- ggplot(
  total_effects,
  aes(x = Predictor, y = Total_effect)
) +
  geom_col(width = 0.65, fill = "grey65", color = "black") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_text(
    aes(label = round(Total_effect, 2)),
    hjust = ifelse(total_effects$Total_effect > 0, -0.15, 1.15),
    size = 4
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Standardized total effect on LnCH4"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank()
  )

ggsave(
  "Standardized_total_effect_on_LnCH4.png",
  p_total,
  width = 6.5,
  height = 4.2,
  dpi = 600
)

ggsave(
  "Standardized_total_effect_on_LnCH4.pdf",
  p_total,
  width = 6.5,
  height = 4.2
)

# ============================================================
# 10. SEM path diagram
# ============================================================

png(
  "SEM_lavaan_path_diagram_updated.png",
  width = 3200,
  height = 2200,
  res = 600
)

semPaths(
  Fit,
  what = "std",
  whatLabels = "std",
  style = "lisrel",
  layout = "tree2",
  residuals = FALSE,
  intercepts = FALSE,
  exoCov = FALSE,
  nCharNodes = 0,
  edge.label.cex = 1.1,
  sizeMan = 9,
  sizeMan2 = 5,
  curvePivot = TRUE,
  edge.color = "black",
  color = "white",
  border.width = 2,
  mar = c(6, 6, 6, 6)
)
semPaths
dev.off()

pdf(
  "SEM_lavaan_path_diagram_updated.pdf",
  width = 11,
  height = 7
)

semPaths(
  Fit,
  what = "std",
  whatLabels = "std",
  style = "lisrel",
  layout = "tree2",
  residuals = FALSE,
  intercepts = FALSE,
  exoCov = FALSE,
  nCharNodes = 0,
  edge.label.cex = 1.1,
  sizeMan = 9,
  sizeMan2 = 5,
  curvePivot = TRUE,
  edge.color = "black",
  color = "white",
  border.width = 2,
  mar = c(6, 6, 6, 6)
)

dev.off()

