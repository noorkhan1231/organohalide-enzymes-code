# Step 13c. Geography: does compositional dissimilarity of the ten families increase with geographic distance within each host system?
# usage: Rscript s10_geography.R <family_rpkg_main.csv> <sample_coordinates.tsv> <out_prefix>
# Samples with site coordinates (ENA run record or BioSample attribute) and at least one target read.
# (1) sample level: Mantel test (Spearman, 999 permutations) between Bray-Curtis dissimilarity and great-circle distance;
# (2) study level (one point per study: mean composition, mean coordinates), which removes within-study pseudo-replication;
# (3) absolute latitude vs family RPKG at study level (Spearman, BH-adjusted).
suppressMessages(library(vegan))
a <- commandArgs(TRUE)
d <- read.csv(a[1], check.names = FALSE); co <- read.delim(a[2]); out <- a[3]
FAM <- c("RDase", "HyDase", "GSTDase", "OxDase", "DehyDase", "MetDase", "FlaHase", "VanHase", "NHFeHase", "SAMHase")
co <- co[co$coord_source %in% c("ENA run record", "BioSample attribute"), c("sample_id", "lat", "lon")]
d <- merge(d, co, by = "sample_id"); d <- d[rowSums(d[, FAM]) > 0, ]
hav <- function(la, lo) { r <- pi / 180; la <- la * r; lo <- lo * r; n <- length(la); m <- matrix(0, n, n)
  for (i in seq_len(n)) { dl <- la - la[i]; dg <- lo - lo[i]
    h <- sin(dl / 2)^2 + cos(la[i]) * cos(la) * sin(dg / 2)^2; m[i, ] <- 2 * 6371 * asin(pmin(1, sqrt(h))) }
  as.dist(m) }
res <- list(); lat <- list()
for (h in c("Chicken", "Coral", "Rhizosphere", "Sponge")) {
  x <- d[d$host_system == h, ]
  if (nrow(x) >= 10) {
    set.seed(42); m <- mantel(vegdist(x[, FAM], "bray"), hav(x$lat, x$lon), method = "spearman", permutations = 999)
    res[[length(res) + 1]] <- data.frame(host = h, level = "sample", n = nrow(x), studies = length(unique(x$study)), r = m$statistic, P = m$signif,
                                         max_km = max(hav(x$lat, x$lon)))
  }
  s <- aggregate(x[, c(FAM, "lat", "lon")], by = list(study = x$study), FUN = mean)
  if (nrow(s) >= 8) {
    set.seed(42); m2 <- mantel(vegdist(s[, FAM], "bray"), hav(s$lat, s$lon), method = "spearman", permutations = 999)
    res[[length(res) + 1]] <- data.frame(host = h, level = "study", n = nrow(s), studies = nrow(s), r = m2$statistic, P = m2$signif, max_km = max(hav(s$lat, s$lon)))
    for (f in FAM) { ct <- suppressWarnings(cor.test(abs(s$lat), s[[f]], method = "spearman", exact = FALSE))
      lat[[length(lat) + 1]] <- data.frame(host = h, family = f, n_studies = nrow(s), rho = unname(ct$estimate), P = ct$p.value) }
  }
}
R <- do.call(rbind, res); R$FDR <- p.adjust(R$P, "BH"); write.csv(R, paste0(out, "_mantel.csv"), row.names = FALSE)
L <- do.call(rbind, lat); L$FDR <- p.adjust(L$P, "BH"); write.csv(L, paste0(out, "_latitude.csv"), row.names = FALSE)
print(R); print(L[L$FDR < 0.05, ])
