# =========================
# Libraries
# =========================
library(terra)
library(dplyr)
library(here)
library(sf)
library(stars)
library(viridis)
library(ade4)
library(factoextra)
library(patchwork)
library(tmap)
library(usdm)
library(NbClust)
library(clue)
library(maptree)
library(purrr)
library(readr)
library(RColorBrewer)
library(ggplot2)
library(tidyr)

# =========================
# Scenario
# =========================
scen <- "SrM_0-200_month"
processed_dir <- here("input/fPCA", scen, "stack/")

# =========================
# Raster settings
# =========================
res <- 0.05
ext <- ext(2, 11, 39.5, 44.7)

month_list <- c("jan","feb","mar","apr","may","jun",
                "jul","aug","sep","oct","nov","dec")

# =========================
# Loop over months
# =========================
for (month in month_list) {
  
  message("Processing: ", month)
  
  # -------------------------
  # 0. Output dir
  # -------------------------

  month_dir <- here("output", "congruence", scen, month)
  dataset_dir <- file.path(month_dir, "dataset")
  cluster_dir <- file.path(month_dir, "cluster")
  scree_dir <- file.path(month_dir, "opti_cluster")
  km_object_dir <- file.path(month_dir, "km_object")
  rast_dir <- file.path(month_dir, "rasters")
  map_dir <- file.path(month_dir, "maps")

  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(scree_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(km_object_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rast_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------
  # 1. Load + preprocess
  # -------------------------
  
  # Cropping
  clim_stack <- rast(file.path(processed_dir,
                               paste0("fPCA_3D_stack_", month, ".tif"))) %>%
    crop(ext)
  
  # Resampling
  ref_rast <- rast(ext(clim_stack),
                resolution = res,
                crs = crs(clim_stack))
  
  clim_resamp <- resample(clim_stack, ref_rast, method = "bilinear")
  
  # Transform in dataframe
  clim_df <- as.data.frame(clim_resamp, xy = TRUE) %>%
    na.omit()
  
  meta <- clim_df %>% select(x, y)
  clim <- clim_df %>% select(-x, -y)
  
  # -------------------------
  # 2. PCA variable selection
  # -------------------------
  pca <- dudi.pca(clim, center = TRUE, scale = TRUE, scannf = FALSE, nf = 5)
  
  contrib <- data.frame(
    variable = names(clim),
    pc1 = fviz_contrib(pca, choice = "var", axes = 1)@data$contrib,
    pc2 = fviz_contrib(pca, choice = "var", axes = 2)@data$contrib
  )
  
  # Select vars with the highest contribution to the two first axis
  selected_vars <- contrib %>%
    filter(pc1 >= quantile(pc1, 0.75, na.rm = TRUE) |
             pc2 >= quantile(pc2, 0.75, na.rm = TRUE)) %>%
    pull(variable)
  
  clim_scale <- scale(clim[, selected_vars, drop = FALSE]) %>%
    as.data.frame()
  
  # -------------------------
  # 3. VIF selection
  # -------------------------
  vif_selected_vars <- vifstep(clim_scale, th = 5)@results$Variables
  
  # To fix the selection of parameters : not useful with fpca outputs (always the same selected)
  
  # vif <- replicate(50, {
  #   res <- vifstep(clim_scale, th = 5)@results$Variables 
  #   paste(sort(res), collapse = ";") 
  #   }) 
  # 
  # all_vars <- lapply(vif, function(x) unlist(strsplit(x, ";"))) %>% unlist() # Calculate the frequency #
  # 
  # freq <- table(all_vars) %>% sort(decreasing = TRUE) # Find the dropping point across frequencies #
  # 
  # drops <- diff(freq) 
  # 
  # cut_index <- which.min(drops) 
  # 
  # vif_selected_vars <- names(freq)[1:cut_index]
  
  data_env <- clim_scale[, vif_selected_vars, drop = FALSE]
  
  # -------------------------
  # Save selected data
  # -------------------------
  
  saveRDS(data_env, file.path(dataset_dir, "data_env.rds"))
  saveRDS(meta, file.path(dataset_dir, "meta.rds"))
  saveRDS(ref_rast, file.path(dataset_dir, "ref_rast.rds"))
  
  # =========================
  # 4. Clustering analysis
  # =========================
  
  # ----------------------------------
  # Fix the maximum number of clusters
  # ----------------------------------
  
  compute_scree <- function(dist_method, clust_method, k_min, k_max, k_primary, k_secondary) {
    
    cat("\nComputing Hierarchical Clustering with", dist_method, "and", clust_method, "...\n")
    
    # =========================
    # Pre-computation
    # =========================
    S <- cov(data_env)
    dist_mat <- dist(data_env, method = dist_method)
    hc <- hclust(dist_mat, method = clust_method)
    
    cat("\nComputing Mahalanobis distances for k =", k_min, "to", k_max, "...\n")
    
    # =========================
    # Main loop over k
    # =========================
    resultats_scree <- map_dfr(k_min:k_max, function(k) {
      
      cluster <- cutree(hc, k)
      
      # Save clusters
      if (k <= k_secondary) {
        saveRDS(
          cluster,
          file.path(
            cluster_dir,
            paste0("hc_", dist_method, "_", clust_method, "_k", k, ".rds")
          )
        )
      }
      
      # Centroids
      centroides <- data_env %>%
        mutate(bioregion = cluster) %>%
        group_by(bioregion) %>%
        summarise(across(everything(), mean), .groups = "drop") %>%
        dplyr::select(-bioregion)
      
      paires <- combn(1:k, 2)
      
      # Mahalanobis distances
      distances <- apply(paires, 2, function(p) {
        sqrt(mahalanobis(
          as.numeric(centroides[p[1], ]),
          as.numeric(centroides[p[2], ]),
          S
        ))
      })
      
      list(
        scree = tibble(
          k = k,
          d_median = median(distances),
          d_mean = mean(distances)
        ),
        distances = tibble(
          k = k,
          cluster_i = paires[1, ],
          cluster_j = paires[2, ],
          distance = distances
        )
      )
      
    }, .progress = TRUE)
    
    # =========================
    # Save outputs
    # =========================
    write_csv(
      resultats_scree[["distances"]],
      file.path(scree_dir,
                paste0("distances_", dist_method, "_", clust_method, ".csv"))
    )
    
    write_csv(
      resultats_scree[["scree"]],
      file.path(scree_dir,
                paste0("screeplot_mahalanobis_", dist_method, "_", clust_method, ".csv"))
    )
    
    # =========================
    # Diagnostics
    # =========================
    delta_d_median <- diff(resultats_scree[["scree"]]$d_median)
    
    cat(sprintf(
      "Selected partitions: k = %d (primary) and k = %d (secondary)\n",
      k_primary, k_secondary
    ))
    
    cat(k_primary, ": primary elbow of d_mean — large-scale provinces\n")
    cat(k_secondary, ": secondary peak of d_mean/d_median — sub-regional structures\n")
    
    # =========================
    # Screeplot
    # =========================
    p_scree <- resultats_scree[["scree"]] %>%
      tidyr::pivot_longer(-k, names_to = "metric", values_to = "distance") %>%
      mutate(metric = factor(metric,
                             levels = c("d_median", "d_mean"),
                             labels = c("Median", "Mean"))) %>%
      ggplot(aes(k, distance, color = metric, linetype = metric)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.5) +
      
      geom_vline(xintercept = k_primary, linetype = "dashed",
                 color = "#E63946", linewidth = 0.7) +
      annotate("text",
               x = k_primary + 0.4,
               y = max(resultats_scree[["scree"]]$d_mean) * 0.98,
               label = paste0("k = ", k_primary, "\n(large provinces)"),
               hjust = 0, size = 3.2, color = "#E63946") +
      
      geom_vline(xintercept = k_secondary, linetype = "dashed",
                 color = "#2D6A4F", linewidth = 0.7) +
      annotate("text",
               x = k_secondary + 0.4,
               y = max(resultats_scree[["scree"]]$d_mean) * 0.98,
               label = paste0("k = ", k_secondary, "\n(sub-regions)"),
               hjust = 0, size = 3.2, color = "#2D6A4F") +
      
      scale_x_continuous(breaks = seq(k_min, k_max, by = 2)) +
      scale_color_manual(values = c("#E63946", "#457B9D")) +
      scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
      
      labs(
        title = "Screeplot of Mahalanobis distances between bioregion centroids",
        subtitle = paste0(
          "Two retained partitions: k = ", k_primary,
          " (primary) and k = ", k_secondary
        ),
        x = "Number of bioregions (k)",
        y = "Mahalanobis distance between centroids",
        color = "Metric",
        linetype = "Metric",
        caption = paste0(
          "Variables: ", ncol(data_env),
          " | n = ", format(nrow(data_env), big.mark = ",")
        )
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "grey90"),
        plot.caption = element_text(color = "grey50", size = 8)
      )
    
    # =========================
    # Delta plot
    # =========================
    p_delta <- tibble(
      k = resultats_scree[["scree"]]$k[-1],
      delta = delta_d_median
    ) %>%
      ggplot(aes(k, delta)) +
      geom_col(fill = "#457B9D", alpha = 0.7) +
      geom_hline(yintercept = 0, color = "grey40") +
      scale_x_continuous(breaks = seq(k_min + 1, k_max, by = 2)) +
      labs(
        title = "Rate of change of median Mahalanobis distance",
        subtitle = expression(paste(Delta, "d_median")),
        x = "k",
        y = expression(Delta * " d_median")
      ) +
      theme_classic(base_size = 12)
    
    # =========================
    # Save figure
    # =========================
    ggsave(
      file.path(
        scree_dir,
        paste0("screeplot_mahalanobis_", dist_method, "_", clust_method, "_", month, ".png")
      ),
      plot = p_scree / p_delta + patchwork::plot_layout(heights = c(2, 1)),
      width = 14,
      height = 10,
      dpi = 300,
      bg = "white"
    )
  }
  
  # Apply compute scree function
  result_scree_eucli_ward <- compute_scree(dist_method = "euclidean", clust_method = "ward.D2", k_min = 2, k_max = 50, k_primary = 2, k_secondary = 15)
  result_scree_manhattan_ward <- compute_scree(dist_method = "manhattan", clust_method = "ward.D2", k_min = 2, k_max = 50, k_primary = 2, k_secondary = 15)
  
  # =========================
  # CAUTION : to change in function of the screeplot
  # =========================
  
  kmax_euc_ward <- 15
  kmax_man_ward <- 15 

  # =========================
  # K-means clustering
  # =========================
  
  k_max <- max(kmax_euc_ward, kmax_man_ward)
  
  for (k in 2:k_max) {
    
    message("k-means k = ", k)
    
    km <- kmeans(data_env, centers = k, nstart = 50, iter.max = 100)
    
    # Save cluster labels
    saveRDS(km$cluster,
            file.path(cluster_dir, paste0("km_", k, ".rds")))
    
    # Save full model
    saveRDS(km,
            file.path(km_object_dir, paste0("km_object_", k, ".rds")))
  }
  
  # =========================
  # Rasterization
  # =========================

  
  cluster_files <- list.files(cluster_dir, pattern = "\\.rds$", full.names = TRUE)
  
  for (f in cluster_files) {
    
    cluster <- readRDS(f)
    
    cluster_df <- cbind(meta, cluster = cluster) %>%
      mutate(clust = as.factor(cluster))
    
    cluster_v <- vect(cluster_df,
                      geom = c("x", "y"),
                      crs = "EPSG:4326")
    
    cluster_rast <- rasterize(cluster_v, ref_rast, field = "clust")
    
    out_name <- paste0(tools::file_path_sans_ext(basename(f)), ".tif")
    
    writeRaster(cluster_rast,
                file.path(rast_dir, out_name),
                overwrite = TRUE)
  }
  
  # =========================
  # Mapping
  # =========================
  
  rast_files <- list.files(rast_dir, pattern = "\\.tif$", full.names = TRUE)
  
  partition_list <- lapply(rast_files, rast)
  names(partition_list) <- tools::file_path_sans_ext(basename(rast_files))
  
  coastline_med <- st_read(here("input/Med_contours/coastline_med.shp"))
  
  palette_cols <- c(
    "#A8E6CF", "#D1B3E0", "#FFE082", "#F4A3A3", "lightblue1",
    "#E0E0E0", "darkseagreen3", "#8C97D6", "#FFD3A5", "#FAD1E0",
    "slategray2", "#B8A9E3", "#80CBC4", "#FFF59D", "#E6CCB2",
    "#6FAED6", "#F48FB1", "pink3", "#B0BEC5", "darkseagreen2",
    "#FF8A65", "#D7CCC8", "#A3C4F3", "thistle1", "lightsalmon"
  )
  
  barplot(rep(1, length(palette_cols)),
          col = palette_cols,
          border = NA,
          space = 0,
          main = "Palette (25 couleurs)",
          axes = FALSE)
  
  # =========================
  # Map function
  # =========================
  
  map_partition <- function(r, title) {
    
    tm_shape(r, is.main = TRUE) +
      tm_raster(
        col.scale = tm_scale_categorical(values = palette_cols),
        col.legend = tm_legend(title = "Regions")
      ) +
      
      tm_shape(coastline_med) +
      tm_polygons(fill = "grey94", col = "black", lwd = 0.2) +
      
      tm_graticules(
        lines = TRUE,
        col = "grey60",
        lwd = 0.5,
        labels.size = 2
      ) +
      
      tm_layout(
        legend.outside = TRUE,
        legend.outside.position = c("right", "top"),
        legend.stack = "horizontal",
        legend.title.size = 1.6,
        legend.text.size = 1.6,
        legend.frame = FALSE,
        inner.margins = 0,
        main.title = title,
        main.title.size = 2.2,
        main.title.position = "center"
      )
  }
  
  # =========================
  # Generate maps
  # =========================
  
  map_list <- lapply(names(partition_list), function(nm) {
    map_partition(partition_list[[nm]], nm)
  })
  
  names(map_list) <- names(partition_list)
  
  # =========================
  # Save maps
  # =========================
  
  for (nm in names(map_list)) {
    
    tmap_save(
      map_list[[nm]],
      filename = file.path(map_dir, paste0("map_", nm, "_", month, ".png")),
      width = 27,
      height = 15
    )
  }
  
}