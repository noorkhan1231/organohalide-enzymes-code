suppressMessages(library(vegan)); a <- commandArgs(TRUE); dir <- a[1]
g <- read.csv(file.path(dir, "groups.csv")); res <- list()
for (k in seq_len(nrow(g))) {
  grp <- g$group[k]; D <- as.matrix(read.csv(file.path(dir, paste0(grp, "_dist.csv")), row.names = 1, check.names = FALSE))
  m <- read.csv(file.path(dir, paste0(grp, "_meta.csv"))); tab <- table(m$host)
  if (length(tab) < 2) next
  d <- as.dist(D); h <- factor(m$host)
  if (nrow(m) <= length(tab)) { res[[length(res) + 1]] <- data.frame(group = grp, n = nrow(m), n_host = length(tab), min_per_host = min(tab), R2 = 1, P = NA, permdisp_P = NA, R2_given_length = NA, P_given_length = NA, R2_length = NA); next }
  set.seed(42); p <- adonis2(d ~ h, permutations = 999)
  set.seed(42); bd <- tryCatch(permutest(betadisper(d, h), permutations = 999), error = function(e) list(tab = data.frame(`Pr(>F)` = NA, check.names = FALSE)))
  r <- data.frame(group = grp, n = nrow(m), n_host = length(tab), min_per_host = min(tab), R2 = p$R2[1], P = p$`Pr(>F)`[1], permdisp_P = bd$tab$`Pr(>F)`[1], R2_given_length = NA, P_given_length = NA, R2_length = NA)
  if (all(!is.na(m$length))) { ll <- log10(m$length); set.seed(42); p2 <- tryCatch(adonis2(d ~ ll + h, permutations = 999, by = "margin"), error = function(e) NULL)
    if (!is.null(p2)) { r$R2_given_length <- p2$R2[2]; r$P_given_length <- p2$`Pr(>F)`[2]; r$R2_length <- p2$R2[1] } }
  res[[length(res) + 1]] <- r }
R <- do.call(rbind, res); write.csv(R, file.path(dir, "permanova.csv"), row.names = FALSE); print(R, digits = 3)
