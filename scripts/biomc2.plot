#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: biomc2.plot [breakpoint-profile.tsv] [recomb_freq.pdf]")
input <- if (length(args)) args[1] else "breakpoint-profile.tsv"
output <- if (length(args) > 1L) args[2] else "recomb_freq.pdf"
d <- read.delim(input)
pdf(output, width = 8, height = 6)
layout(matrix(seq_len(3)), heights = c(0.16, 1, 1))
par(mar = rep(0, 4))
plot.new()
legend("center", c(paste("Mean", d$metric[1], "distance"), "Breakpoint probability", "Pooled intervals"),
       col = c("grey55", "#0072B2", "#B34700"), lty = c(1, NA, 1),
       pch = c(NA, 16, 18), bty = "n", horiz = TRUE, cex = 0.8)
par(mar = c(4, 4, 2, 1))
for (group in c("all", "modal")) {
  p <- d[d$group == group, ]
  if (!nrow(p)) {
    plot.new()
    title("No candidate breakpoints")
    next
  }
  p <- p[order(p$position), ]
  prob <- p$break_weight / p$weight
  dist <- p$distance_sum / p$weight
  plot(p$position, dist, type = "h", col = "grey55", ylim = c(0, max(1, dist)),
       xlim = c(1, p$alignment_length[1]), xlab = "Alignment position",
       ylab = "Posterior mean", main = if (group == "all") "All samples" else
         sprintf("Samples with %d breakpoints", p$modal_count[1]))
  points(p$position, prob, pch = 16, col = "#0072B2", cex = 0.55)
  # Partition the cumulative breakpoint mass into one interval per modal breakpoint.
  mode <- p$modal_count[1]
  if (mode > 0L && sum(prob) > 0) {
    mass <- cumsum(prob) / sum(prob)
    at <- function(q) p$position[which(mass >= q)[1L]]
    for (j in seq_len(mode) - 1L) {
      x <- vapply((j + c(0.025, 0.5, 0.975)) / mode, at, numeric(1))
      segments(x[1], 0.96 * max(1, dist), x[3], 0.96 * max(1, dist), col = "#B34700")
      points(x[2], 0.96 * max(1, dist), pch = 18, col = "#B34700")
    }
  }
  if (any(prob == 0)) rug(p$position[prob == 0], col = "grey70")
}
invisible(dev.off())
