# Regenerate the two demonstration figures from CSV evidence using base R only.
bundle_root <- normalizePath(Sys.getenv("INAPEST_BUNDLE_ROOT", unset = "."), mustWork = TRUE)
results_dir <- Sys.getenv("INAPEST_TEST_OUT", unset = file.path(bundle_root, "results", "rerun"))
fig_dir <- file.path(bundle_root, "figures", "rerun")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

x <- read.csv(file.path(results_dir, "mortality_by_fecundity_results.csv"))
png(file.path(fig_dir, "mortality_fecundity_effect.png"), width=1400, height=900, res=160)
par(mar=c(5.2,5.4,2.5,1.5), las=1)
morts <- sort(unique(x$mortality)); frs <- sort(unique(x$fecundity_reduction))
ymax <- max(x$q95) * 1.06
plot(range(frs), c(0,ymax), type="n", xlab="Fecundity reduction under management", ylab="Final total population", xaxt="n")
axis(1, at=frs, labels=frs)
for (i in seq_along(morts)) {
  z <- x[x$mortality==morts[i],]
  lines(z$fecundity_reduction, z$final_population_mean, type="b", pch=15+i, lty=i, lwd=2)
  arrows(z$fecundity_reduction, z$q05, z$fecundity_reduction, z$q95, angle=90, code=3, length=.04)
}
legend("topright", legend=paste("Mortality =", morts), lty=seq_along(morts), pch=15+seq_along(morts), lwd=2, bty="n")
dev.off()

lu <- read.csv(file.path(results_dir, "landuse_mortality_fecundity_results.csv"))
z <- lu[lu$mortality==0.3 & lu$pattern=="landuse_gradient",]
vals <- c(z$node1_LU1_mean, z$node2_LU2_mean, z$node3_LU3_mean)
fr <- c(z$FR_LU1, z$FR_LU2, z$FR_LU3)
png(file.path(fig_dir, "landuse_fecundity_gradient_effect.png"), width=1200, height=800, res=160)
par(mar=c(5.2,5.4,2.5,1.5), las=1)
bp <- barplot(vals, names.arg=paste0("Land use ",1:3,"\nFR = ",fr), ylab="Mean final population", ylim=c(0,max(vals)*1.2))
text(bp, vals, labels=format(round(vals,1), nsmall=1), pos=3)
dev.off()
cat("Figures written to", fig_dir, "\n")
