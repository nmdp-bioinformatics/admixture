library(sf)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(scales)

# ---------- settings ----------
# Usage: Rscript chloropleth.R [freq_dir] [target_haplo]
# Example: Rscript chloropleth.R freqs "A*33:01~C*08:02~B*14:02~DRB1*01:02"
args <- commandArgs(trailingOnly = TRUE)
freq_dir <- if (length(args) >= 1) args[1] else "freqs"
target_haplo <- if (length(args) >= 2) args[2] else stop("Usage: Rscript chloropleth.R <freq_dir> <target_haplo>")

# ISO3 lookup for your codes + IL-suffix countries
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

# ---------- helper: parse key from filename ----------
parse_key <- function(fn) {
  # freqs.CY.csv -> CY
  # freqs.IL-LIBYA.csv -> IL-LIBYA
  # freqs.US-NAFR.cav -> US-NAFR
  key <- str_replace(basename(fn), "^freqs\\.", "")
  key <- str_replace(key, "\\.[^.]+$", "")
  key
}

# ---------- read all files, extract frequency for target haplotype ----------
files <- list.files(freq_dir, pattern = "^freqs\\.", full.names = TRUE)

hap_freqs <- lapply(files, function(f) {
  key <- parse_key(f)
  
  # Skip unspecified anchor-like populations for choropleth (optional)
  if (key == "US-NAFR") return(NULL)
  
  df <- suppressMessages(read_csv(f, show_col_types = FALSE))
  
  # Expected columns: Haplo, Count, Freq (extras ok)
  if (!all(c("Haplo", "Freq") %in% names(df))) {
    warning("Skipping (missing Haplo/Freq): ", f)
    return(NULL)
  }
  
  freq <- df %>%
    filter(Haplo == target_haplo) %>%
    summarise(freq = sum(as.numeric(Freq), na.rm = TRUE)) %>%
    pull(freq)
  
  if (length(freq) == 0) freq <- 0
  
  tibble::tibble(population = key, freq = freq)
}) %>% bind_rows()

# ---------- map population keys to ISO3 ----------
hap_freqs2 <- hap_freqs %>%
  mutate(
    suffix = ifelse(str_detect(population, "^IL-"),
                    str_replace(population, "^IL-", ""),
                    NA_character_),
    plot_key = case_when(
      str_detect(population, "^IL-") & suffix %in% iso_lookup$key ~ suffix,
      # IL-ARAB / IL-DRUZE are not countries -> Israel
      str_detect(population, "^IL-") ~ "ISRAEL",
      TRUE ~ population
    )
  ) %>%
  left_join(iso_lookup, by = c("plot_key" = "key")) %>%
  filter(!is.na(iso3))

# If multiple populations map to the same country (e.g., Israel + IL-ARAB),
# take the max (or sum). Max is usually better for “presence/intensity” maps.
country_freq <- hap_freqs2 %>%
  group_by(iso3) %>%
  summarise(freq = max(freq, na.rm = TRUE), .groups = "drop")

# ---------- base map ----------
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(4326)

bbox <- st_bbox(c(xmin = -12, ymin = 20, xmax = 65, ymax = 48), crs = st_crs(4326))
world_crop <- st_crop(world, bbox)

map_df <- world_crop %>%
  left_join(country_freq, by = c("iso_a3" = "iso3")) 
# %>%  mutate(freq = ifelse(is.na(freq), 0, freq))

# ---------- plot ----------
p <- ggplot(map_df) +
  geom_sf(aes(fill = freq), color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-12, 65), ylim = c(20, 48), expand = FALSE) +
  scale_fill_viridis_c(
    trans = "sqrt",  # nicer for rare haplotypes; try "log10" if you prefer
    labels = label_number(accuracy = 0.0001),
    name = "Freq"
  ) +
  labs(
    title = "Choropleth of haplotype frequency",
    subtitle = target_haplo
  ) +
  theme_void() +
  theme(legend.position = "right")

print(p)
out_file <- paste0("haplo_choropleth_", gsub("[*:~]", "_", target_haplo), ".png")
ggsave(out_file, p, width = 12, height = 7, dpi = 300)
