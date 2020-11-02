# -------------------- load libraries ---------------
library(tidyverse)
library(broom)
library(emmeans)
library(ggfortify)
library(performance)
library(gt)
library(grid)
library(gridExtra)
library(png)
#  ---------- set plotting themes -----------

theme_set(theme_bw() + 
            theme(legend.title=element_blank(),
                  legend.position = 'bottom',
                  legend.key.size = unit(0.75, "cm"),
                  text = element_text(size=12),
                  strip.background =element_rect(fill="grey99")))


# ----------- Load Data -------------------------

zones_df <- read_csv('./data/CWC_can_change_df.csv') %>%
  filter(time_step != "Dec16 - Mar18") %>%
  filter(time_step != "Dec16 - Mar18") %>%
  filter(time_step !="Feb17 - Mar18") %>%
  filter(time_step !="Jan18 - Mar18") %>%
  filter(time_step !="Dec16 - Feb17") %>%
  filter(time_step !="Feb17 - Jan18")


# --------------- Box Plots ----------------------------
# create df for regions where any change is observed
All_df <- zones_df %>%
  mutate(loss_gain = 'ALL') %>%
  filter(canopy_change != 0)

# create df with positive and negative change uniquely labeled and bind to All_df
Box_p_df <- zones_df %>%
  mutate(loss_gain = ifelse(canopy_change <0, "LOSS", ifelse(canopy_change > 0, "GAIN", "NO_CHANGE"))) %>%
  filter(loss_gain != "NO_CHANGE") %>%
  mutate(canopy_change_abs = abs(canopy_change)) %>%
  bind_rows(All_df) %>%
  mutate(loss_gain = fct_relevel(loss_gain, 'LOSS', 'GAIN'))%>%
  mutate(signs_YNf = fct_relevel(signs_YNf, 'No Foraging'))


# function to plot boxplot without outliers and maintain scale. From: https://stackoverflow.com/questions/25124895/no-outliers-in-ggplot-boxplot-with-facet-wrap
calc_boxplot_stat <- function(x) { 
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}


# this adds colours to the right hand facets
colour_right_facets <- function(plot){
  g <- ggplot_gtable(ggplot_build(plot))
  stripr <- which(grepl('strip-r', g$layout$name))
  fills <- c("grey90","grey70","grey50")
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  grid.arrange(g)
}

# Veg change box plots
change_plot <- function(.data){
p <- ggplot(.data, aes(x=signs_YNf, y= canopy_change, colour=signs_YNf, fill=signs_YNf)) +
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", alpha = 0.8, colour='#000000') +
  facet_grid(loss_gain~ time_step,  scales = "free") +
  labs(y=bquote('Canopy Elevation change  '~(m/m^2)))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual("Zone", values=c('#377eb8', '#ff7f00')) 

colour_right_facets(p)
}


change_raste_plot <- change_plot(Box_p_df)

ggsave('NewPlots/boxplot_by_change.jpg', change_raste_plot, dpi=600, width=6, height = 6)



# ---------- Regression Analysis ------------


reg_tab <- function(.data){ 
  
  tab_col <- 'grey80'
  
  lg <- pull(.data, var = loss_gain)[1]
  if (lg == 'LOSS'){
    tit <- 'Loss Regions Regression (Gamma GLM):'
    tab_col <- 'grey90'
  } else if (lg == 'GAIN'){
    tit <- 'Gain Regions Regression (Gamma GLM):'
    tab_col <- 'grey70'
  } else {
    tit <- 'Net Change Regression (LM):'
    tab_col <- 'grey50'
  }
  
  gt_tab <- .data %>%
    group_by(time_step)%>%
    select(!loss_gain) %>%
    gt() %>%
    tab_style(style = cell_fill(tab_col), locations = cells_title()) %>%
    tab_style(style = cell_fill('grey95'), locations = cells_row_groups()) %>%
    tab_header(title = md(sprintf('**<div style="text-align: left"> %s </div>**', tit))) %>%
    gtsave(.,tempfile('tab', fileext = '.png')) %>%
    readPNG(.) %>%
    rasterGrob( interpolate=TRUE, width = unit(10,"cm"))
  return(gt_tab)
}

join.hori <- function(.data) {
  x11()
  g <- do.call(grid.arrange, c(.data, ncol=3))
  dev.off()
  return(g)
}

# regression tibble for Net change
group_regr_1 <- Box_p_df %>% 
  filter(loss_gain=='ALL')%>%
  group_by(loss_gain,time_step) %>%
  do(fitForage = tidy(lm(formula = canopy_change ~ signs_YNf, data = .))) %>% 
  unnest(fitForage) %>%
  mutate_if(is.numeric, round, 5) 

# regression for positive and negative change then join to net change and create gt table
group_regr_all <- Box_p_df %>% 
  filter(loss_gain!='ALL')%>%
  group_by(loss_gain, time_step) %>%
  do(fitForage = tidy(glm(formula = canopy_change_abs ~ signs_YNf, data = ., family=Gamma(link = 'identity')))) %>% 
  unnest(fitForage) %>%
  mutate_if(is.numeric, round, 3) %>%
  bind_rows(group_regr_1) %>%
  mutate(estimate = ifelse(loss_gain == 'LOSS', estimate *-1, estimate)) %>%
  mutate(statistic = ifelse(loss_gain == 'LOSS', statistic *-1, estimate))%>%
  mutate(term = ifelse(term=='signs_YNfForaging Observed', 'Foraging Observed', term)) %>%
  mutate(p.value = ifelse(p.value < 0.001, '< 0.001 **', 
                          ifelse(p.value < 0.05, paste(formatC(p.value,format = "f", 3), '*', sep = " "),
                                 ifelse(p.value < 0.1, paste(formatC(p.value,format = "f", 3), '.', sep = " "),
                                        formatC(p.value,format = "f", 3))))) %>%
  rename(T.statistic = statistic) %>%
  group_by(loss_gain) %>%
  group_split() %>%
  purrr::map(., ~reg_tab(.)) %>%
  join.hori(.) %>%
  grid.arrange() %>%
  ggsave(filename = 'NewPlots/regression_summary.jpg', dpi = 600, height = 10, width = 30, units = 'cm')


# --------------- Model Diagnostics --------------------------------

diag.plots <- function(.model){
  par(mfrow = c(2, 2))
  resids <- resid(.model)
  hist(resids)
  plot(.model, which = 2)
  plot(.model, which = 1)
  plot(.model, which = 5)
  par(mfrow = c(1, 1))
}

diag.plots2 <- function(.model){
  autoplot(.model, which = 1:6, ncol = 3, label.size = 3) 
}


check_gamma_reg <- function(.data) {
  reg <- glm(formula = canopy_change ~ signs_YNf, data = .data, family=Gamma(link = 'identity'))
  lg <- pull(.data, var = loss_gain)[1]
  ts <- pull(.data, var = time_step)[1]
  message(sprintf("%s zone regression for %s", lg, ts))
  print(tidy(reg))
  diag.plots2(reg)
}
check_gauss_reg <- function(.data) {
  reg <- lm(formula = canopy_change ~ signs_YNf, data = .data)
  lg <- pull(.data, var = loss_gain)[1]
  ts <- pull(.data, var = time_step)[1]
  message(sprintf("%s zone regression for %s", lg, ts))
  print(tidy(reg))
  diag.plots2(reg)
}

#check gain and loss regressions

Box_p_df %>% 
  filter(loss_gain!='ALL')%>%
  mutate(canopy_change = abs(canopy_change)) %>%
  group_by(loss_gain, time_step) %>%
  group_split() %>%
  purrr::map(~check_gamma_reg(.))

# check net change regressions
Box_p_df %>% 
  filter(loss_gain=='ALL')%>%
  group_by(time_step) %>%
  group_split() %>%
  purrr::map(~check_gauss_reg(.))
