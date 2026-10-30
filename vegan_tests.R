# vegan tests on one abundance table.  usage: Rscript vegan_tests.R <table.csv> <id_col> <group_col> <block_col|-> <descriptor_col|-> <out_prefix> <nmds_k|0>
# table: one row per sample, abundance columns = every column not named in id/group/block/descriptor
suppressMessages(library(vegan))
a <- commandArgs(TRUE); d <- read.csv(a[1], check.names = FALSE, stringsAsFactors = FALSE)
idc <- a[2]; gc <- a[3]; bc <- a[4]; dc <- a[5]; out <- a[6]; k <- as.integer(a[7])
meta <- c(idc, gc, bc, dc); X <- d[, setdiff(names(d), meta), drop = FALSE]
keep <- rowSums(X) > 0; d <- d[keep, ]; X <- X[keep, , drop = FALSE]
g <- factor(d[[gc]]); res <- list()
bray <- vegdist(X, "bray")
set.seed(42); p <- adonis2(bray ~ g, permutations = 999, parallel = 2)
res[[length(res) + 1]] <- data.frame(test = "host_sample", n = nrow(X), df = p$Df[1], R2 = p$R2[1], F = p$F[1], P = p$`Pr(>F)`[1])
set.seed(42); bd <- betadisper(bray, g); pt <- permutest(bd, permutations = 999, parallel = 2)
res[[length(res) + 1]] <- data.frame(test = "permdisp_sample", n = nrow(X), df = pt$tab$Df[1], R2 = NA, F = pt$tab$F[1], P = pt$tab$`Pr(>F)`[1])
an <- anova(bd); res[[length(res) + 1]] <- data.frame(test = "permdisp_sample_parametric", n = nrow(X), df = an$Df[1], R2 = NA, F = an$`F value`[1], P = an$`Pr(>F)`[1])
write.csv(data.frame(group = names(tapply(bd$distances, g, mean)), mean_dist = as.numeric(tapply(bd$distances, g, mean))), paste0(out, "_dispersion.csv"), row.names = FALSE)
mode <- Sys.getenv("VEGAN_MODE", "full")
# note: a host test with permutations restricted within study is not informative here (each study belongs to one host), so it is not run
if (bc != "-") {
  b <- factor(d[[bc]])
  if (mode == "full") {
    set.seed(42); p2 <- adonis2(bray ~ g / b, permutations = 999, by = "terms", parallel = 2)
    res[[length(res) + 1]] <- data.frame(test = "block_nested", n = nrow(X), df = p2$Df[2], R2 = p2$R2[2], F = p2$F[2], P = p2$`Pr(>F)`[2])
  }
  st <- aggregate(X, by = list(g = d[[gc]], b = d[[bc]]), FUN = mean); SX <- st[, -(1:2)]
  b2 <- vegdist(SX, "bray"); set.seed(42); p4 <- adonis2(b2 ~ g, data = st, permutations = 9999, parallel = 2)
  res[[length(res) + 1]] <- data.frame(test = "host_study", n = nrow(SX), df = p4$Df[1], R2 = p4$R2[1], F = p4$F[1], P = p4$`Pr(>F)`[1])
  set.seed(42); pt2 <- permutest(betadisper(b2, factor(st$g)), permutations = 999, parallel = 2)
  res[[length(res) + 1]] <- data.frame(test = "permdisp_study", n = nrow(SX), df = pt2$tab$Df[1], R2 = NA, F = pt2$tab$F[1], P = pt2$tab$`Pr(>F)`[1])
}
if (dc != "-") {
  ds <- factor(paste(d[[gc]], d[[dc]]))
  set.seed(42); p5 <- adonis2(bray ~ g / ds, permutations = 999, by = "terms", parallel = 2)
  res[[length(res) + 1]] <- data.frame(test = "descriptor_nested", n = nrow(X), df = p5$Df[2], R2 = p5$R2[2], F = p5$F[2], P = p5$`Pr(>F)`[2])
}
write.csv(do.call(rbind, res), paste0(out, "_tests.csv"), row.names = FALSE)
if (k > 0) {
  set.seed(42); nm <- metaMDS(X, distance = "bray", k = k, trymax = 20, trace = 0, parallel = 2)
  sc <- as.data.frame(scores(nm, display = "sites")); sc$id <- d[[idc]]; sc$group <- d[[gc]]; sc$stress <- nm$stress
  write.csv(sc, paste0(out, "_nmds.csv"), row.names = FALSE)
}
