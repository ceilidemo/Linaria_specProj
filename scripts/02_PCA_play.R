## 02_PCA : Tryna piece out PCA... what are we seeinggg
# Ceili DeMarais

library(dplyr)
library(ggplot2)
library(plotly)

# setup dat
datasets <- list(
  raw = list(dat = "data_work/spec_clean.csv", work = "data_work/PLSDA/raw"),
  processed = list(dat = "data_work/spec_clean_processed.csv", work = "data_work/PLSDA/processed")
)

for (dataset_name in names(datasets)) {
  
  dat  <- datasets[[dataset_name]]$dat
  work <- datasets[[dataset_name]]$work
  
  cat("\n\n", dataset_name, "\n")
  
  dati <- read.csv(dat) %>%
    mutate(
      treatment = factor(treatment),
      leafPos   = factor(leafPos),
      timePt    = factor(timePt)
    ) %>%
    dplyr::filter(leafPos == "bot", timePt == "data_1")
  
  wvl_cols    <- names(dati)[grepl("^X[0-9]", names(dati))]
  wvl_seq     <- wvl_cols[seq(1, length(wvl_cols), by = 10)]
  
  spec <- dati[, wvl_seq]
  
  # PCA
  pca_fit <- prcomp(spec, scale. = TRUE)
  
  # gt scores
  scores_df <- as.data.frame(pca_fit$x[, 1:3])
  scores_df$treatment <- dati$treatment
  
  # long format
  scores_long <- data.frame(
    PC1 = c(scores_df$PC1, scores_df$PC1, scores_df$PC1),
    PC2 = c(scores_df$PC2, scores_df$PC3, scores_df$PC3),
    PCx = c(rep("PC2", nrow(scores_df)), rep("PC3", nrow(scores_df)), rep("PC3", nrow(scores_df))),
    PCy = c(rep("PC1", nrow(scores_df)), rep("PC1", nrow(scores_df)), rep("PC2", nrow(scores_df))),
    facet = c(rep("PC1 vs PC2", nrow(scores_df)), rep("PC1 vs PC3", nrow(scores_df)), rep("PC2 vs PC3", nrow(scores_df))),
    treatment = rep(scores_df$treatment, 3)
  )
  
  cat("Variance explained PC1:", round(summary(pca_fit)$importance[2, 1] * 100, 1), "%\n")
  cat("Variance explained PC2:", round(summary(pca_fit)$importance[2, 2] * 100, 1), "%\n")
  cat("Variance explained PC3:", round(summary(pca_fit)$importance[2, 3] * 100, 1), "%\n")
  
  # Facet 2d plots
  pdf(file.path(work, "PCA_data1_faceted.pdf"), width = 14, height = 5)
  p1 <- ggplot(scores_long, aes(x = PC1, y = PC2, color = treatment, fill = treatment)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(aes(color = treatment), geom = "polygon", alpha = 0.2, type = "norm") +
    facet_wrap(~facet, scales = "free") +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 11, face = "bold")) +
    labs(title = paste("PCA data_1 -", dataset_name),
         x = "PC Axis", y = "PC Axis", color = "Treatment", fill = "Treatment")
  print(p1)
  dev.off()
  
  # 3d scatter interactive
  p2 <- plot_ly(scores_df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~treatment,
                type = "scatter3d", mode = "markers",
                marker = list(size = 5, opacity = 0.7)) %>%
    layout(title = paste("PCA 3D -", dataset_name),
           scene = list(xaxis = list(title = paste("PC1 (", round(summary(pca_fit)$importance[2, 1] * 100, 1), "%)", sep = "")),
                        yaxis = list(title = paste("PC2 (", round(summary(pca_fit)$importance[2, 2] * 100, 1), "%)", sep = "")),
                        zaxis = list(title = paste("PC3 (", round(summary(pca_fit)$importance[2, 3] * 100, 1), "%)", sep = ""))))
  
  htmlwidget_path <- file.path(work, "PCA_data1_3D.html")
  htmlwidgets::saveWidget(p2, htmlwidget_path)
  
  # Scree
  pdf(file.path(work, "PCA_data1_scree.pdf"), width = 7, height = 6)
  var_exp <- summary(pca_fit)$importance[2, 1:10]
  scree_df <- data.frame(PC = 1:10, var_exp = var_exp * 100)
  p3 <- ggplot(scree_df, aes(x = PC, y = var_exp)) +
    geom_point(size = 3) +
    geom_line() +
    theme_bw() +
    labs(title = paste("Scree Plot -", dataset_name),
         x = "Principal Component", y = "Variance Explained (%)")
  print(p3)
  dev.off()
  
  cat("Plots saved to:", work, "\n")
}

# 3D trajectory plot
for (dataset_name in names(datasets)) {
  
  dat  <- datasets[[dataset_name]]$dat
  work <- datasets[[dataset_name]]$work
  
  cat("\n\n", dataset_name, "\n")
  
  dati <- read.csv(dat) %>%
    mutate(
      treatment = factor(treatment),
      leafPos   = factor(leafPos),
      timePt    = factor(timePt)
    ) %>%
    dplyr::filter(leafPos == "bot")
  
  wvl_cols    <- names(dati)[grepl("^X[0-9]", names(dati))]
  wvl_seq     <- wvl_cols[seq(1, length(wvl_cols), by = 10)]
  
  # Perform PCA on all data combined
  spec_all <- dati[, wvl_seq]
  pca_fit <- prcomp(spec_all, scale. = TRUE)
  
  # Project all samples onto PCA space
  pca_scores <- as.data.frame(pca_fit$x[, 1:3])
  pca_scores$treatment <- dati$treatment
  pca_scores$timePt <- dati$timePt
  
  # Calculate mean trajectory for each treatment over time
  trajectory <- pca_scores %>%
    group_by(treatment, timePt) %>%
    summarise(PC1 = mean(PC1),
              PC2 = mean(PC2),
              PC3 = mean(PC3),
              .groups = "drop") %>%
    arrange(treatment, timePt)
  
  # Create numeric time ordering
  time_order <- data.frame(timePt = unique(as.character(trajectory$timePt)))
  time_order$time_num <- 1:nrow(time_order)
  
  trajectory <- trajectory %>%
    mutate(timePt = as.character(timePt)) %>%
    left_join(time_order, by = "timePt")
  
  # Get unique treatments
  treatments <- unique(trajectory$treatment)
  
  # Create color palette
  colors <- scales::hue_pal()(length(treatments))
  color_map <- setNames(colors, treatments)
  
  # Create 3D trajectory plot
  p_trajectory <- plot_ly()
  
  for (i in seq_along(treatments)) {
    trt <- treatments[i]
    traj_subset <- trajectory %>% filter(treatment == trt) %>% arrange(time_num)
    
    # Add line trajectory
    p_trajectory <- p_trajectory %>%
      add_trace(data = traj_subset,
                x = ~PC1, y = ~PC2, z = ~PC3,
                type = "scatter3d", mode = "lines+markers",
                line = list(color = color_map[trt], width = 4),
                marker = list(size = 8, opacity = 0.8),
                name = trt,
                showlegend = TRUE,
                hoverinfo = "text",
                text = ~paste("Treatment:", treatment, "<br>Time:", timePt,
                              "<br>PC1:", round(PC1, 2),
                              "<br>PC2:", round(PC2, 2),
                              "<br>PC3:", round(PC3, 2)))
  }
  
  p_trajectory <- p_trajectory %>%
    layout(title = paste("PCA 3D Trajectory Over Time -", dataset_name),
           scene = list(
             xaxis = list(title = paste("PC1 (", round(summary(pca_fit)$importance[2, 1] * 100, 1), "%)", sep = "")),
             yaxis = list(title = paste("PC2 (", round(summary(pca_fit)$importance[2, 2] * 100, 1), "%)", sep = "")),
             zaxis = list(title = paste("PC3 (", round(summary(pca_fit)$importance[2, 3] * 100, 1), "%)", sep = ""))),
           height = 700,
           width = 900)
  
  htmlwidget_path <- file.path(work, "PCA_3D_trajectory.html")
  htmlwidgets::saveWidget(p_trajectory, htmlwidget_path)
  
  cat("3D trajectory plot saved to:", htmlwidget_path, "\n")
}
