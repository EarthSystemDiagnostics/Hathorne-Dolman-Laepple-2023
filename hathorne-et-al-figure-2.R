# Code to reproduce Figure 2 in the commentary "Seasonal and inter-annual
# marine climate proxy data" Ed Hathorne, Andrew M. Dolman, and Thomas Laepple


# Required R packages -----

## All packages from CRAN ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(sedproxy)

# Monthly ERSST for 1900 – 2000 as used in Napier et al. (2022)
# ERSST v5, 24 N, 66 E
ersst_Napier <- read.csv("ersst_Napier.csv")


# Sedproxy on monthly resolution  -----

## Illustration of effect of mixing in the water column before settling 
## assuming ~3 months mixing

# sediment accumulation rate of 1.5 mm / year in units of cm/kyr
s_cm_kyr <- 1.5 * 1000 * 0.1

# months of mixing in water column
n_months_mix <- 3


## input timeseries for forward modelling alkenone proxy -------

# reverse time series so youngest is "on top" like a sediment core 
clim.in <- ts(cbind(rev(ersst_Napier$SST_ERA5)),
              deltat = 1, start = min(ersst_Napier$Age.months)
)

# time axis
tpts <- ersst_Napier$Age.months

# measurement error for Alkenone in units of temperature
sigma.alk <- 0.23


## Sets the random number generator seed so that the results are the same each time
## Change this to get a different (random) result
set.seed(20230822)
sed1 <- ClimToProxyClim(clim.in, tpts,
                        bio.depth = 3 * s_cm_kyr / 1000, 
                        sed.acc.rate = s_cm_kyr, 
                        layer.width = 0,
                        plot.sig.res = 1,
                        sigma.meas = sigma.alk)



# Create dataframe of results and plot --------
ersst2 <- ersst_Napier %>% 
  filter(Year > 1900, Year < 2000) %>% 
  mutate(timepoints = Age.months) %>% 
  left_join(., sed1$simulated.proxy) %>% 
  select(Year, Month, Decimal.year, timepoints, SST_ERA5, reconstructed.climate) %>% 
  rename(ERSST = SST_ERA5, 
         `Forward modelled SST` = reconstructed.climate) %>% 
  pivot_longer(cols = c(ERSST, `Forward modelled SST`))


fig_ersst_sedproxy <- ersst2 %>% 
  ggplot(aes(x = `Month`, y = value, group = Year, colour = Year)) +
  geom_line() +
  scale_colour_viridis_c(limits = c(1900, 2000), option = "viridis") +
  labs(y = "SST [°C]") +
  scale_x_continuous(breaks = seq(2, 12, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(~name)

fig_ersst_sedproxy

ggsave(filename = "fig.ersst.sedproxy.png", fig_ersst_sedproxy, width = 6, height = 4, dpi = 300)
ggsave(filename = "fig.ersst.sedproxy.pdf", fig_ersst_sedproxy, width = 6, height = 4)

# for svg the package svglite is required
# ggsave(filename = "fig.ersst.sedproxy.svg", fig_ersst_sedproxy, width = 6, height = 4)
