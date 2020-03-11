# INstall Stan
# See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started 


# # #

# #
# # # This step is optional, but it can result in compiled Stan programs that execute much faster than they otherwise would. Simply paste the following into R once
# dotR <- file.path(Sys.getenv("HOME"), ".R")
# if (!file.exists(dotR)) dir.create(dotR)
# M <- file.path(dotR, "Makevars.win")
# if (!file.exists(M)) file.create(M)
# cat("\nCXX14FLAGS=-O3 -march=corei7 -mtune=corei7",
#     "CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y",
#     "CXX11FLAGS=-O3 -march=corei7 -mtune=corei7",
#     file = M, sep = "\n", append = TRUE)
# # 
# # # check if rstan exists and delete it.
# remove.packages("rstan")
# if (file.exists(".RData")) file.remove(".RData")
# 
# .rs.restartR() # restart R
# install.packages("rstan", type = "source")
# # # # # install
# # install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
# # # #If this line ultimately returns TRUE, then your C++ toolchain is properly installed and you can jump to the next section.
# pkgbuild::has_build_tools(debug = TRUE)

# list.of.packages <- c("pkgbuild", "parallel", "dplyr", "rstanarm", "reshape2", "sjPlot", "extrafont", 
#                       "bayestestR", "cowplot", "ggrepel", 'htmlTable', "tidyverse", "brms", "tidybayes")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)


library(tidybayes)
library(htmlTable)
library(pkgbuild)
library(parallel)
library(dplyr)
library(rstanarm)
library(reshape2)
library(sjPlot)
library(patchwork)
library(extrafont)
library(bayestestR)
require(cowplot)
library(ggrepel)
library(tidyverse)


#  ------ load rstan ----------------
library(rstan)
options(mc.cores = parallel::detectCores())      # enabling parallel processing
rstan_options(auto_write = TRUE)                 # allows you to automatically save a bare version of a compiled Stan program to the hard disk so that it does not need to be recompiled (unless you change it).
# Sys.setenv(LOCAL_CPPFLAGS = '-march=native')     # Not needed if the comment out section above is run...

# if errors occur running models try:
# file.edit("~/.R/Makevars")
# And enter the following: CXX14 = g++ -std=c++1y -Wno-unused-variable -Wno-unused-function -fPIC

# ---- ggplot theme setting ------
# font_add("Arial", "Arial")

# showtext_auto() 
# windows()

theme_set(theme_bw() + 
            theme(legend.title=element_blank(),
                  legend.position = 'bottom',
                  legend.key.size = unit(0.75, "cm"),
                  text = element_text(size=12)))

cbbPalette <- c("grey50", "grey10")
darkPalette <- c( "#458DB6", "#C48802")

# ----- Load Data ----------

zones_df = read_csv('./data/CWC_can_change_df.csv')
tsNames = unique(zones_df$time_step)
ts.list <- zones_df %>%
  group_by(time_step) %>%
  group_split()

# ------ Set up functions ------------
bin_bayes <- dget('Bin_Bayes_function.R')

add_ts_function <- function(stanout, idx){
  df <- as.data.frame(stanout) %>%
    select (-lp__) %>%
    mutate(time_step = tsNames[idx])
}


# - run model for gains -----
StanOutGains <- lapply(ts.list, function(x) bin_bayes(x, 'gain'))

gg_list <- list()
for (i in 1:length(StanOutGains )){
  gg_list[[i]] <- traceplot(StanOutGains [[i]]) +
    labs(title = paste(tsNames[i], 'gains'))
}

(gg_list[[1]]+gg_list[[2]])/
 (gg_list[[3]]+gg_list[[4]])/
 (gg_list[[5]]+gg_list[[6]])/
 (gg_list[[7]])


OutGains <- mapply(add_ts_function, StanOutGains, seq_along(StanOutGains), SIMPLIFY = FALSE) %>%
  bind_rows()%>%
  melt()%>%
  rename(Zone = variable) %>%
  mutate(CategoryB = ifelse(Zone=='Yes_Beaver_Rate', 'Foraging Present', 'Foraging Absent' )) %>%
  mutate(time_step = fct_rev(time_step))


Gains_plot <- ggplot(data = OutGains, aes(y= time_step ,x=value, fill=CategoryB, color=CategoryB)) + 
  stat_halfeyeh(alpha=0.6, size=1.5) +
  scale_fill_manual("Zone", values=darkPalette, labels=c('Foraging Absent', 'Foraging Present')) +
  scale_color_manual("Zone",values=cbbPalette, labels=c('Foraging Absent', 'Foraging Present') )+
  # coord_cartesian(xlim = c(0, 0.7)) +
  labs(x = 'Probability of canopy gain')+
  theme(axis.title.y=element_blank())
Gains_plot
ggsave('R_plots/gain_prob_All_TS.jpg', plot = Gains_plot, dpi=300, width = 174, height = 150, units="mm")

# - run model for loss --------

StanOutLoss <- lapply(ts.list, function(x) bin_bayes(x, 'loss'))

gg_list2 <- list()
for (i in 1:length(StanOutLoss)){
  gg_list2[[i]] <- traceplot(StanOutLoss [[i]]) +
    labs(title = paste(tsNames[i], 'loss'))
}

(gg_list2[[1]]+gg_list2[[2]])/
  (gg_list2[[3]]+gg_list2[[4]])/
  (gg_list2[[5]]+gg_list2[[6]])/
  (gg_list2[[7]])


OutLoss<- mapply(add_ts_function, StanOutLoss, seq_along(StanOutLoss), SIMPLIFY = FALSE) %>%
  bind_rows()%>%
  melt()%>%
  rename(Zone = variable) %>%
  mutate(CategoryB = ifelse(Zone=='Yes_Beaver_Rate', 'Foraging Present', 'Foraging Absent' )) %>%
  mutate(time_step = fct_rev(time_step))


Loss_plot <- ggplot(data = OutLoss, aes(y= time_step ,x=value, fill=CategoryB, color=CategoryB)) + 
  stat_halfeyeh(alpha=0.6, size=1.5) +
  scale_fill_manual("Zone", values=darkPalette, labels=c('Foraging Absent', 'Foraging Present')) +
  scale_color_manual("Zone",values=cbbPalette, labels=c('Foraging Absent', 'Foraging Present') )+
  # coord_cartesian(xlim = c(0, 0.7)) +
  labs(x = 'Probability of canopy loss')+
  theme(axis.title.y=element_blank())
Loss_plot
ggsave('R_plots/loss_prob_All_TS.jpg', plot = Loss_plot, dpi=300, width = 174, height = 150, units="mm")
