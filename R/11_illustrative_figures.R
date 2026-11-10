# =============================================================================
# 11_illustrative_figures.R
# -----------------------------------------------------------------------------
# Illustrations only (no classification, no results):
#   Figure 1    one example network from each generative model, n = 100,
#               with the parameters given in the Figure 1 caption;
#   Fig. S12    SBM panels, n = 100, k = 2-6 communities, lower and higher
#               between-community probability (supplement, network comparison
#               figure), same design as sbm-sim-different-communities.R.
# Each panel is one random draw, so a rerun gives a different-looking network
# with the same parameters.
#
# Outputs: figures/{er,sbm,sf,sw,sp}_plot22.png
#          figures/n-100-k-{2..6}-{lower,higher}.png
# =============================================================================

source(here::here("R", "00_setup.R"))
set.seed(22)

save_net <- function(g, file, colour = "#56B4E9") {
  grDevices::png(here::here("figures", file), width = 900, height = 900, res = 150)
  graphics::par(mar = rep(0.2, 4))
  plot(g, vertex.size = 4, vertex.label = NA, vertex.color = colour,
       edge.color = grDevices::adjustcolor("grey40", 0.4), layout = igraph::layout_with_fr(g))
  grDevices::dev.off()
  message(sprintf("%-28s nodes %d, edges %d", file, igraph::vcount(g), igraph::ecount(g)))
}

## ---- Figure 1 (parameters from the caption) ----------------------------------
n <- 100
save_net(igraph::sample_gnp(n, p = 0.24), "er_plot22.png")
## SBM: two blocks of 0.4n and 0.6n. The caption does not give the connection
## probabilities; 0.4 within and 0.08 between give about 1,210 edges.
b <- c(0.4 * n, 0.6 * n)
save_net(igraph::sample_sbm(n, pref.matrix = matrix(c(0.4, 0.08, 0.08, 0.4), 2), block.sizes = b),
         "sbm_plot22.png")
save_net(igraph::sample_pa(n, power = 2, m = 13, directed = FALSE), "sf_plot22.png")
save_net(igraph::sample_smallworld(dim = 1, size = n, nei = 3, p = 0.01), "sw_plot22.png")
save_net(igraph::sample_grg(n, radius = 0.32), "sp_plot22.png")

## ---- SBM panels: n = 100, k = 2-6 ---------------------------------------------
sbm_draw <- function(k, q_range) {
  pm <- matrix(0, k, k)
  for (i in 1:(k - 1)) for (j in (i + 1):k) pm[i, j] <- pm[j, i] <- stats::runif(1, q_range[1], q_range[2])
  diag(pm) <- stats::runif(k, 0.5, 0.8)
  sizes <- rep(floor(n / k), k); sizes[1] <- sizes[1] + n - sum(sizes)
  g <- igraph::sample_sbm(n, pref.matrix = pm, block.sizes = sizes)
  list(g = g, community = rep(seq_len(k), times = sizes))
}
pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")
for (sc in list(c("lower", 0.01, 0.04), c("higher", 0.1, 0.4))) {
  for (k in 2:6) {
    d <- sbm_draw(k, as.numeric(sc[2:3]))
    save_net(d$g, sprintf("n-100-k-%d-%s.png", k, sc[1]), colour = pal[d$community])
  }
}

message("11_illustrative_figures.R complete.")
