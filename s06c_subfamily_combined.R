# Combined subfamily composition across all families: Bray-Curtis on RPKG of 55 labelled subfamilies, PCoA, PERMANOVA, PERMDISP.
# usage: Rscript s06c_subfamily_combined.R <FINAL dir>
suppressMessages(library(vegan))
b <- commandArgs(TRUE)[1]; set.seed(1)
L <- read.csv(file.path(b, 'results/s04_build/subfamily_rpkg_long_main.csv')); L <- L[L$subfamily != 'unassigned', ]
F <- read.csv(file.path(b, 'results/s04_build/family_rpkg_main.csv'))
L$k <- paste(L$family, L$subfamily, sep = '::')
X <- xtabs(RPKG ~ sample_id + k, data = L); X <- unclass(X); X <- X[rowSums(X) > 0, , drop = FALSE]
host <- F$host_system[match(rownames(X), F$sample_id)]; study <- F$study[match(rownames(X), F$sample_id)]
d <- vegdist(X, 'bray'); a <- adonis2(d ~ host, permutations = 999); bd <- betadisper(d, host); pb <- permutest(bd, permutations = 999)
pc <- cmdscale(d, k = 3, eig = TRUE); ev <- pc$eig[pc$eig > 0]; ev <- ev / sum(ev)
out <- data.frame(sample_id = rownames(X), host_system = host, PCo1 = pc$points[, 1], PCo2 = pc$points[, 2], PCo3 = pc$points[, 3])
write.csv(out, file.path(b, 'results/s06_subfamily/combined_pcoa_coordinates.csv'), row.names = FALSE)
# study level: study means
Xs <- rowsum(X, study) / as.vector(table(study)[rownames(rowsum(X, study))]); hs <- host[match(rownames(Xs), study)]
ds <- vegdist(Xs, 'bray'); as <- adonis2(ds ~ hs, permutations = 9999); bs <- permutest(betadisper(ds, hs), permutations = 999)
res <- data.frame(level = c('sample', 'study'), n = c(nrow(X), nrow(Xs)), n_subfamilies = ncol(X), R2 = c(a$R2[1], as$R2[1]), F = c(a$F[1], as$F[1]),
  P = c(a$`Pr(>F)`[1], as$`Pr(>F)`[1]), permdisp_F = c(pb$tab$F[1], bs$tab$F[1]), permdisp_P = c(pb$tab$`Pr(>F)`[1], bs$tab$`Pr(>F)`[1]),
  PCo1_var = ev[1], PCo2_var = ev[2])
write.csv(res, file.path(b, 'results/s06_subfamily/combined_permanova.csv'), row.names = FALSE); print(res)
