# Redraws manuscript Figure 3a: the 158 analysed streams on the map of
# Switzerland (streams without a complete winter-summer DOC pair are dropped,
# so map and analysis show the same set).
suppressMessages({library(tidyverse); library(here)})

# the 158 analysed streams (paired winter + summer sampling), as in the plan
m <- read.csv(here("data", "processed", "m12.csv"))
w <- m |> filter(season == "winter") |> select(site, beaver_territory)
s <- m |> filter(season == "summer") |> select(site, beaver_territory)
m12 <- read.csv(here("data", "processed", "m12.csv"))
w2 <- m12 |> filter(season == "winter") |> select(site, beaver_territory, up_w = DOC_input, dn_w = DOC_out)
s2 <- m12 |> filter(season == "summer") |> select(site, beaver_territory, up_s = DOC_input, dn_s = DOC_out)
d <- inner_join(w2, s2, by = c("site", "beaver_territory")) |>
  filter(complete.cases(up_w, dn_w, up_s, dn_s))

# coordinates and the producer-abundance subset from m6
m6 <- read.csv(here("data", "processed", "m6.csv"))
co <- m6 |> distinct(site, .keep_all = TRUE) |> select(site, lon = x_upstream, lat = y_upstream)
# NOTE: the 16 producer-abundance streams are a separate campaign; their
# coordinates are not in the repository (x/y are NA in m6.csv), so they cannot
# be drawn yet -- add their coordinates to data/processed to place the black dots.
pts <- d |> distinct(site) |> left_join(co, by = "site")
cat("DOC streams on map:", nrow(pts), "\n")

# Switzerland outline from the maps package
ch <- map_data("world", region = "Switzerland")
p <- ggplot() +
  geom_polygon(data = ch, aes(long, lat, group = group), fill = "grey96", colour = "grey40", linewidth = 0.4) +
  geom_point(data = pts, aes(lon, lat), colour = "red", size = 1.7, alpha = 0.85) +
  coord_quickmap() +
  annotate("text", x = 5.99, y = 47.75, label = "(a)", size = 6.5, hjust = 0) +
  theme_void()
ggsave(here("results", "figures", "plan", "fig_map_158.png"), p, width = 8, height = 5.2, dpi = 300, bg = "white")
