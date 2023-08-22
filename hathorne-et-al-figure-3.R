# Code to reproduce Figure 3 in the commentary "Seasonal and inter-annual
# marine climate proxy data" Ed Hathorne, Andrew M. Dolman, and Thomas Laepple

# Required R packages -------

# available on CRAN
library(dplyr)
library(tidyr)
library(RSEIS)
library(ggplot2)

# Only required to create the composite figure
library(patchwork)

# R package PaleoSpec is available from GitHub
# https://github.com/EarthSystemDiagnostics/paleospec

# install.packages("remotes")
# remotes::install_github("EarthSystemDiagnostics/paleospec")
library(PaleoSpec)

## R package pfields is available from GitHub
# https://github.com/EarthSystemDiagnostics/pfields
# remotes::install_github("EarthSystemDiagnostics/pfields")
library(pfields)


# Custom functions ------

##' Peak searching developed for the isotope cycle project
##'
##' @title Analyse the peaks of a timeseries
##' @param depth depth in m
##' @param observable parameter to be analysed, e.g. O18
##' @param res Resolution in m on which the analysis will take place
##' @param spans Spans = minimum width of the peak, real width = res*spans
##' @param bPlot provide debug plots
##' @return List with cycle length (depth in m, length in m), and yMax and yMin, each one (depth in m, max/min value in the original units)
##' @author Thomas Laepple
AnalysePeaks <- function(depth, observable, res = 1, spans = 3, bPlot = FALSE, ...) {
  observable <- scale(observable, scale = FALSE)

  # add numerical noise to avoid missing plateaus caused by rounding errors
  observable <- observable + runif(length(observable), -1e-5, 1e-5)

  # Approximate to resolution
  temp <- approx(depth, observable, seq(from = min(depth), to = max(depth), by = res))


  maxV <- temp$y[which(peaks(temp$y, span = spans))]
  minV <- temp$y[which(peaks(-1 * temp$y, span = spans))]

  pMax <- temp$x[which(peaks(temp$y, span = spans))]
  pMin <- temp$x[which(peaks(-1 * temp$y, span = spans))]

  if (bPlot) {
    plot(depth, observable, type = "l", ...)
    abline(v = pMax, col = "red")

    abline(h = 0, col = "grey", lty = 2)
  }

  maxPos <- (pMax[-1] + pMax[-length(pMax)]) / 2
  minPos <- (pMin[-1] + pMin[-length(pMin)]) / 2


  cyclelength <- cbind(c(minPos, maxPos), c(diff(pMin), diff(pMax)))
  yMax <- cbind(c(pMax), c(maxV))
  yMin <- cbind(c(pMin), c(minV))

  ts.filtered <- data.frame(depth = depth, observable = observable)

  return(list(
    cyclelength = cyclelength, yMax = yMax, yMin = yMin,
    ts.filtered = ts.filtered, peaks = data.frame(peaks = pMax)
  ))
}


# Function that does the filtering and peak counting --------

#' Filter timeseries, count peaks and tune age model to peaks
#'
#' @param vec timeseries
#' @param filt a filter (numeric vector of weights)
#' @param ntYear number of timepoints sampled per year
#'
#' @return a list
FilterAndCountPeaks <- function(vec, filt, ntYear = 10) {

  # Filter the timeseries and keep middle 2/3 to minimize artefacts on the edges
  # ~100 years = 1000 points in this example

  n <- length(vec)

  ts.bp <- window(stats::filter(ts(vec), filt), n * 1/6, n * 5/6)

  # also cut the raw time-series to the same time period
  ts.raw <- window(ts(vec), n * 1/6, n * 5/6)

  # Get the maxima
  result <- AnalysePeaks(seq(ts.bp), c(ts.bp), bPlot = FALSE)

  n_peaks <- nrow(result$yMax)

  print(paste0(n_peaks, " peaks detected"))

  ## Indices where the peaks have been found
  index.max <- result$yMax[, 1]

  ## Interpolation to new time-scale
  newAge <- approx(index.max, (1:length(index.max)) * ntYear - (ntYear / 2), seq(ts.raw), rule = 2)$y
  newY <- approx(newAge, c(ts.raw), seq(ts.bp), rule = 2)$y

  ts.tuned <- pTs(newY)

  AnalysePeaks(seq(ts.bp), c(ts.bp),
               bPlot = FALSE,
               main = "bandpass filtered with peak detection", xlab = "depth units"
  )

  # power spectra of timeseries
  spec.raw <- SpecMTM(ts(ts.raw, deltat = 1/ntYear), plot = FALSE)
  spec.tuned <- SpecMTM(ts(ts.tuned, deltat = 1/ntYear), plot = FALSE)

  timeseries <- data.frame(
    time = 1:length(ts.tuned),
    ts.tuned = as.numeric(ts.tuned),
    ts.raw = as.numeric(ts.raw)
  )

  out <- list(
    ts.tuned = ts.tuned,
    ts.raw = ts.raw,
    timeseries = timeseries,
    spec.tuned = spec.tuned,
    spec.raw = spec.raw,
    result = result
  )
  return(out)
}



# Example with a simulated power-law noise timeseries --------

## Sets the random number generator seed so that the results are the same each time
## Change this to get a different (random) result
set.seed(20230129)

# Assumptions and parameters
# 150 yr of simulated data with a powerlaw of beta = 0.7 as estimated from their figure
# 10 samples per year
# 1.5 mm per year sediment accumulation rate

N <- 150
ntYear <- 10
acc_rate <- 1.5

# Simulate a timeseries
raw <- SimPowerlaw(beta = 0.7, N = N * ntYear)

# Make a symmetric bandpass filter around 1yr = 1/10 depth units in a similar
# width as in Napier et al.

bp <- Bandpass(omega.upper = 0.7 / (ntYear), omega.lower = 1.3 / (ntYear), n = 51)

dat <- FilterAndCountPeaks(vec = raw, filt = bp, ntYear = ntYear)


## plot the results ----

### The raw simulated timeseries
fig.raw <- dat$timeseries %>%
  ggplot(aes(x = time / ntYear * acc_rate, y = ts.raw)) +
  geom_line(alpha = 1, colour = "#1b9e77") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  expand_limits(x = 100) +
  labs(x = "Depth [mm]", y = "Raw values")

### Peaks detected after Bandpass filtering
fig_peaks <- dat$result$ts.filtered %>%
  ggplot(aes(x = depth / ntYear * acc_rate, y = observable)) +
  geom_line() +
  geom_vline(
    data = dat$result$peaks / ntYear * acc_rate, aes(xintercept = peaks),
    colour = "red", lty = 1, alpha = 0.5
  ) +
  labs(colour = "", x = "Depth [mm]", y = "Bandpass filtered") +
  theme_bw() +
  expand_limits(x = 100) +
  theme(panel.grid = element_blank(), legend.position = "none")

### The simulated timeseries after tuning to the detected peaks
fig.tuned <- dat$timeseries %>%
  ggplot(aes(x = time / (ntYear), y = ts.tuned)) +
  geom_line(alpha = 1, colour = "#d95f02") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Tuned age model [years]", y = "Interpolated values")

fig.spec <- gg_spec(dat[4:5]) +
  geom_vline(xintercept = 1, lty = 2, alpha = 0.5) +
  coord_cartesian(x = c(0.2, 2)) +
  scale_colour_brewer("", labels = c("Raw", "Tuned"), palette = "Dark2") +
  theme(legend.position = "right")


fig.ts.comp <- patchwork::wrap_plots(A = fig.raw,
                                     B = fig_peaks,
                                     C = fig.tuned,
                                     D = fig.spec,
                                     ncol = 1) +
  patchwork::plot_annotation(tag_levels = "a")

fig.ts.comp

ggsave("fig.ts.comp.png", fig.ts.comp, height = 9, width = 7, dpi = 300)
ggsave("fig.ts.comp.pdf", fig.ts.comp, height = 9, width = 7, dpi = 300)

# for svg the package svglite is required
# ggsave("fig.ts.comp.svg", fig.ts.comp, height = 9, width = 7, dpi = 300)
