# Step 9b. NMDS (Bray-Curtis, k = 2, metaMDS, 20 random starts, no transformation) of subfamily composition per family,
# samples with reads on the family; used for Supplementary Figure 7. usage: Rscript s06b_subfamily_nmds.R <s06 dir>
suppressMessages(library(vegan)); dir <- commandArgs(TRUE)[1]
for (f in c("RDase", "HyDase", "GSTDase", "OxDase", "MetDase", "FlaHase", "VanHase", "NHFeHase")) {
  d <- read.csv(file.path(dir, paste0("main_", f, "_vegan_input.csv")), check.names = FALSE)
  M <- as.matrix(d[, -(1:3)]); k <- rowSums(M) > 0; d <- d[k, ]; M <- M[k, , drop = FALSE]
  set.seed(42); m <- metaMDS(M, distance = "bray", k = 2, trymax = 20, autotransform = FALSE, trace = 0, parallel = 2)
  out <- data.frame(sample_id = d$sample_id, host_system = d$host_system, NMDS1 = m$points[, 1], NMDS2 = m$points[, 2], stress = m$stress)
  write.csv(out, file.path(dir, paste0("main_", f, "_nmds.csv")), row.names = FALSE); cat(f, nrow(out), round(m$stress, 3), "\n")
}
