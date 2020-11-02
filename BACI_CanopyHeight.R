# BACI analysis using GLM on canopy height daata

library(tidyverse)
library(broom)
library(emmeans)
library(patchwork)
library(gt)
library(grid)
library(gridExtra)
library(png)

# ---- ggplot theme setting ------

theme_set(theme_bw() + 
            theme(legend.title=element_blank(),
                  legend.position = 'bottom',
                  legend.key.size = unit(0.75, "cm"),
                  text = element_text(size=12),
                  strip.background =element_rect(fill="grey99")))


# ------ Import and Clean Data ----------

zones_df <- read_csv('./data/CWC_can_heights_df.csv') %>%
  filter(canopy_height > 0)

BACI_df <- zones_df %>%
  filter(time_step == "Dec16" | time_step == "Jan18" | time_step == "Sep17" | time_step == "Sep18") %>%
  mutate(signs_YNf = fct_relevel(signs_YNf, 'No Foraging')) %>%
  mutate(season = ifelse(time_step == "Dec16" | time_step =="Jan18", 'Dec16 - Jan18',
                         ifelse(time_step == "Sep17" | time_step =="Sep18", 'Sep17 - Sep18', NA)))

# Summer_BACI_df <- zones_df %>%
#   filter(time_step == "Sep17" | time_step == "Sep18") %>%
#   mutate(time = ifelse(time_step=="Sep17",0,1)) %>%
#   mutate(signs_YNf = fct_relevel(signs_YNf, 'No Foraging')) %>%
#   mutate(season = "Summer")

# combined_BACI_df <- Winter_BACI_df %>%
#   bind_rows()

# ---------- Plot Data --------------



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

compare_canopies_plot <- function(.data){
  p <- .data %>%
    ggplot(., aes(x= canopy_height, group=time_step, fill=time_step)) +
    geom_density(colour = "black", alpha = 0.8)+
    facet_grid(season~signs_YNf) +
    labs(x = "Canopy Height  (m)", y='Density') +
    # scale_fill_brewer(type="qual", palette = "Dark2", direction=-1)
    scale_fill_manual(values=c('#d95f02', '#1b9e77', '#e6ab02', "#7570b3"))
  
  colour_right_facets(p)
  
}
# "#56B4E9", "#E69F00"
grid_plot <- compare_canopies_plot(BACI_df)

ggsave('NewPlots/CanopyHeightCompare.jpg', grid_plot, dpi=600, width= 6, height = 6)

# ---------- Regression ----------------

fit.model <- function(.data){
  glm(canopy_height ~ signs_YNf * time_step, family=Gamma(link='identity'),data = .data)
}


diag.plots <- function(.model){
  par(mfrow = c(2, 2))
  # plot(.model)
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


reg_sum_tab <- function(.model, tit, tab.col){
  
  tidy(.model) %>%
    mutate_at(vars(estimate, std.error, statistic), round,3) %>%
    mutate(p.value = ifelse(p.value < 0.001, '< 0.001 **', 
                            ifelse(p.value < 0.05, paste(formatC(p.value,format = "f", 3), '*', sep = " "),
                                   ifelse(p.value < 0.1, paste(formatC(p.value,format = "f", 3), '.', sep = " "),
                                          formatC(p.value,format = "f", 3))))) %>%
    rename(T.statistic = statistic) %>%
    mutate(term = ifelse(term=='signs_YNfForaging Observed', 'Foraging Observed', term)) %>%
    mutate(term = ifelse(term=='time_stepSep18', 'Sep 2018', term)) %>%
    mutate(term = ifelse(term=='time_stepJan18 ', 'Jan 2018', term)) %>%
    mutate(term = ifelse(term=='signs_YNfForaging Observed:time_stepSep18', 'Foraging Observed:Sep 2018', term)) %>%
    mutate(term = ifelse(term=='signs_YNfForaging Observed:time_stepJan18', 'Foraging Observed:Jan 2018', term)) %>%
    gt()%>%
    tab_header(title = md(sprintf('**<div style="text-align: left"> %s </div>**', tit))) %>%
    tab_style(style = cell_fill(tab.col), locations = cells_title())%>%
    gtsave(.,tempfile('tab', fileext = '.png')) %>%
    readPNG(.) %>%
    rasterGrob( interpolate=TRUE, width = unit(5,"cm"))
  
}

emmeans_tab <- function(.model, tit, tab.col){
  tidy(emmeans(.model, ~signs_YNf * time_step)) %>%
    mutate_at(vars(estimate, std.error, z.ratio), round,3) %>%
    select(-c(df, p.value)) %>%
    rename(Date = time_step) %>%
    rename(Zone = signs_YNf) %>%
    gt()%>%
    tab_header(title = md(sprintf('**<div style="text-align: left"> %s </div>**', tit))) %>%
    tab_style(style = cell_fill(tab.col), locations = cells_title()) %>%
    gtsave(.,tempfile('tab', fileext = '.png')) %>%
    readPNG(.) %>%
    rasterGrob( interpolate=TRUE, width = unit(5,"cm"))
}

join_tabs_vert <- function(tab_list){
    x11()
    g <- do.call(grid.arrange, tab_list)
    dev.off()
    return(grid.arrange(g))
}

reg <- BACI_df %>%
  group_by(season) %>%
  group_split()%>%
  purrr::map(., ~fit.model(.))

winter_reg <- reg[[1]]
summer_reg <- reg[[2]]

wint_reg_tab <- reg_sum_tab(reg[[1]], "Regression Summary: Dec 2016 - Jan 2018", 'grey90')
sum_reg_tab <- reg_sum_tab(reg[[2]], "Regression Summary: Sep 2017 - Sep 2018", 'grey70')


BACI_comb_tab <- join_tabs_vert(list(wint_reg_tab, sum_reg_tab)) %>%
  ggsave(., filename = 'NewPlots/BACI_regSumm.jpg', dpi = 600, height = 7, width = 6, units = 'cm')

wint_marmean_tab <- emmeans_tab(reg[[1]], "Marginal Means: Dec 2016 - Jan 2018", 'grey90')
sum_marmean_tab <- emmeans_tab(reg[[2]], "Marginal Means: Sep 2017 - Sep 2018", 'grey70')


BACI_margmeans <- join_tabs_vert(list(wint_marmean_tab, sum_marmean_tab)) %>%
  ggsave(., filename = 'NewPlots/BACI_MarMeans.jpg', dpi = 600, height = 7, width = 6, units = 'cm')

# ---------- Regression Diagnostics
summary(winter_reg)

purrr::map(reg, diag.plots)

diag.plots(winter_reg)
tidy(winter_reg)
tidy(emmeans(winter_reg, ~signs_YNf * time_step))


hist(Summer_BACI_df$canopy_height)
