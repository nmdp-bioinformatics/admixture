library(sf)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(rnaturalearth)
library(scatterpie)
library(stringr)

# -------- SETTINGS --------
K_values <- 2:12
out_prefix <- "nmf_out_k"

# -------- Load map data ONCE --------
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(4326)

bbox <- st_bbox(c(xmin = -12, ymin = 20, xmax = 65, ymax = 48), crs = st_crs(4326))
world_crop <- st_crop(world, bbox)

centroids <- world %>%
  mutate(centroid = st_centroid(geometry),
         lon = st_coordinates(centroid)[,1],
         lat = st_coordinates(centroid)[,2]) %>%
  st_drop_geometry() %>%
  select(iso_a3, name_long, lon, lat)

iso_lookup <- tibble::tribble(
  ~key,          ~iso3,
  "CY",          "CYP",
  "HR",          "HRV",
  "ES",          "ESP",
  "FR",          "FRA",
  "GR",          "GRC",
  "IT",          "ITA",
  "TR",          "TUR",
  "IL",          "ISR",
  "LIBYA",       "LBY",
  "MOROCCO",     "MAR",
  "TUNISIA",     "TUN",
  "IRAQ",        "IRQ",
  "IRAN",        "IRN",
  "YEMEN",       "YEM",
  "ISRAEL",      "ISR"
)

# -------- helper: cosine-order populations --------
cosine_order <- function(W_df) {
  M <- as.matrix(W_df %>% select(starts_with("K")))
  rownames(M) <- W_df$population
  M <- M / sqrt(rowSums(M^2))
  S <- M %*% t(M)
  hc <- hclust(as.dist(1 - S), method = "average")
  hc$labels[hc$order]
}

# -------- LOOP --------
for (K in K_values) {
  
  cat("Processing K =", K, "\n")
  
  W_path <- paste0(out_prefix, K, "/W.populations_by_K.csv")
  if (!file.exists(W_path)) {
    cat("  missing:", W_path, " — skipping\n")
    next
  }
  
  W <- read_csv(W_path, show_col_types = FALSE)
  
  k_cols <- grep("^K\\d+$", names(W), value = TRUE)
  
  # ======================
  # stacked bar chart
  # ======================
  
  pop_order <- cosine_order(W)
  
  W_long <- W %>%
    pivot_longer(all_of(k_cols), names_to = "component", values_to = "prop") %>%
    mutate(
      population = factor(population, levels = pop_order),
      component  = factor(component, levels = k_cols)
    )
  
  p_bar <- ggplot(W_long, aes(population, prop, fill = component)) +
    geom_col(width = 0.9) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = paste0("NMF admixture proportions (K=", K, ")"),
      x = NULL, y = "Mixture proportion"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank())
  
  ggsave(
    paste0("nmf_stackedbars_K", K, ".png"),
    p_bar, width = 10, height = 8, dpi = 300
  )
  
  # ======================
  # map pie chart prep
  # ======================
  
  W_loc <- W %>%
    mutate(
      suffix = ifelse(str_detect(population, "^IL-"),
                      str_replace(population, "^IL-", ""),
                      NA_character_),
      plot_key = case_when(
        str_detect(population, "^IL-") & suffix %in% iso_lookup$key ~ suffix,
        TRUE ~ population
      ),
      label = population
    ) %>%
    left_join(iso_lookup, by = c("plot_key" = "key")) %>%
    mutate(
      iso3 = ifelse(is.na(iso3) & str_detect(population, "^IL-"), "ISR", iso3)
    ) %>%
    left_join(centroids, by = c("iso3" = "iso_a3"))
  
  # US-NAFR anchor
  W_loc <- W_loc %>%
    mutate(
      lon = ifelse(population == "US-NAFR", 2.0, lon),
      lat = ifelse(population == "US-NAFR", 35.0, lat)
    )
  
  # deterministic offsets for overlapping Israel subpops
  W_loc <- W_loc %>%
    mutate(
      lon = case_when(
        population == "IL-ARAB"  ~ lon - 0.9,
        population == "IL-DRUZE" ~ lon + 0.9,
        TRUE ~ lon
      ),
      lat = case_when(
        population == "IL-ARAB"  ~ lat - 0.9,
        population == "IL-DRUZE" ~ lat + 0.9,
        TRUE ~ lat
      )
    )
  
  W_loc <- W_loc %>% filter(!is.na(lon), !is.na(lat)) %>%
    mutate(r = 0.8)
  
  # ======================
  # map pie chart
  # ======================
  
  p_map <- ggplot() +
    geom_sf(data = world_crop, linewidth = 0.3, fill = "white") +
    scatterpie::geom_scatterpie(
      data = W_loc,
      aes(x = lon, y = lat, r = r),
      cols = k_cols,
      alpha = 0.9
    ) +
    geom_text(
      data = W_loc,
      aes(x = lon, y = lat, label = label),
      size = 4,
      vjust = -2
    ) +
    coord_sf(xlim = c(-12, 65), ylim = c(20, 48), expand = FALSE) +
    labs(
      title = paste0("NMF admixture map (K=", K, ")"),
      subtitle = "Pie = component mixture"
    ) +
    theme_void()
  
  ggsave(
    paste0("nmf_piemap_K", K, ".png"),
    p_map, width = 12, height = 7, dpi = 300
  )
  
}

# ============================================================
# Combined PDF: stacked bars (left, small) + pie map (right, large)
# ============================================================
library(png)
library(grid)
library(gridExtra)

pdf("nmf_combined_K2-12.pdf", width = 16, height = 9)

for (K in K_values) {
  bar_file <- paste0("nmf_stackedbars_K", K, ".png")
  map_file <- paste0("nmf_piemap_K", K, ".png")

  if (!file.exists(bar_file) || !file.exists(map_file)) {
    cat("  Skipping K =", K, "— missing PNG(s)\n")
    next
  }

  bar_img <- rasterGrob(readPNG(bar_file), interpolate = TRUE)
  map_img <- rasterGrob(readPNG(map_file), interpolate = TRUE)

  # Layout: stacked bars 1/3 width (left), pie map 2/3 width (right)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 3)))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(bar_img)
  popViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2:3))
  grid.draw(map_img)
  popViewport()

  popViewport()
}

dev.off()
cat("Written: nmf_combined_K2-12.pdf\n")
