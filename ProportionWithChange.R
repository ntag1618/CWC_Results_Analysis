# # Proportion of Change in area
#
# # -------------------- load libraries ---------------
library(tidyverse)
library(gt)
library(grid)
library(gridExtra)
library(png)
library(broom)
library(ggfortify)
#
# #  ---------- set plotting themes -----------

theme_set(theme_bw() +
            theme(legend.title=element_blank(),
                  legend.position = 'bottom',
                  legend.key.size = unit(0.75, "cm"),
                  text = element_text(size=12),
                  strip.background =element_rect(fill="grey99")))


#
# # ----------- Load Data -------------------------
#
zones_df <- read_csv('./data/CWC_can_change_df.csv') %>%
  filter(time_step != "Dec16 - Mar18") %>%
  filter(time_step != "Dec16 - Mar18") %>%
  filter(time_step !="Feb17 - Mar18") %>%
  filter(time_step !="Jan18 - Mar18") %>%
  filter(time_step !="Dec16 - Feb17") %>%
  filter(time_step !="Feb17 - Jan18")



# -------------- Summarise Data and plot--------------------

Prop_change <- zones_df %>%
  mutate(loss_gain = ifelse(canopy_change <0, "LOSS", ifelse(canopy_change > 0, "GAIN", "NO_CHANGE"))) %>%
  group_by(time_step, signs_YNf, loss_gain) %>%
  summarise(Area_m2 = n()*0.25) %>%
  ungroup()%>%
  group_by(time_step, signs_YNf) %>%
  mutate(Zone_area_m2 = sum(Area_m2)) %>%
  ungroup() %>%
  mutate(Perc_Area = Area_m2/Zone_area_m2*100) %>%
  mutate(loss_gain = fct_relevel(loss_gain, 'LOSS', 'GAIN'))%>%
  mutate(signs_YNf = fct_relevel(signs_YNf, 'No Foraging')) 

Prop_change %>%
  filter(loss_gain != "NO_CHANGE") %>%
  ggplot(., aes(x=loss_gain, y = Perc_Area, fill = signs_YNf)) +
  geom_bar(width=0.8, stat="identity", position=position_dodge(width=0.9), colour='black', alpha=0.8)+
  geom_text(aes(label=paste(round(Perc_Area, 1), "%")),position=position_dodge(width=0.95), vjust=1.5, size=3) +
  facet_wrap(~time_step) +
  labs(y='Area with Elevation Change (%)', x= 'Direction of Elevation Change') +
  scale_fill_manual("Zone", values=c('#AC4EC8', '#60C84E')) 

ggsave('NewPlots/Amount_of_change.jpg', dpi=600, width=6, height = 6)


# --------- Create summary table of areas with gain loss and no change for each group. ------------------- 

make_table <- function(.data){

  .data %>%
    mutate_if(is.numeric, round, 1) %>%
    gt(rowname_col = "loss_gain") %>%
    summary_rows(
      groups = TRUE,
      columns = vars(Area_m2),
      fns = list(total = "sum")) %>%
    tab_style(style = cell_fill('grey80'), locations = cells_title()) %>%
    tab_style(style = cell_fill('grey95'), locations = cells_row_groups()) %>%
    cols_label(
      Area_m2 = md('Area (m<sup>2</sup>)'),
      Perc_Area = md('% Area')) # %>%
    # tab_header(title = md(sprintf('**<div style="text-align: left"> %s </div>**',))) 
}


Prop_change%>%
  select(!Zone_area_m2)%>%
  group_by(time_step, signs_YNf) %>%
  make_table() %>%
  gtsave(., filename = normalizePath('NewPlots/Areas_Summary.png', mustWork=FALSE))
  # gtsave(., filename = normalizePath('NewPlots/Areas_Summary.html', mustWork=FALSE))

# ------------ Logistic regression ------------------


log_reg <- function(.data){
  ts <- pull(.data, var = time_step)[1]
  
  .mod1 <- tidy(glm(formula = lossTRUE ~ signs_YNf, data = .data, family = binomial(link = "logit")), conf.int=TRUE, exponentiate=TRUE) %>%
    mutate(time_step = ts) %>%
    mutate(loss_gain = 'LOSS')
  
  # message(sprintf("GAIN: Logistic Regression Table for %s", ts ))
  .mod2 <- tidy(glm(formula = gainTRUE ~ signs_YNf, data = .data, family = binomial(link = "logit")), conf.int=TRUE, exponentiate=TRUE)%>%
    mutate(time_step = ts) %>%
    mutate(loss_gain = 'GAIN')
  
  comb.tab <- bind_rows(.mod1, .mod2)
  
  return(comb.tab)
  
}


LR_zones_df <- zones_df %>%
  mutate(signs_YNf = fct_relevel(signs_YNf, 'No Foraging')) %>%
  mutate(lossTRUE = ifelse(canopy_change < 0, 1, 0)) %>%
  mutate(gainTRUE = ifelse(canopy_change > 0, 1, 0)) %>%
  group_by(time_step)%>%
  group_split()%>%
  purrr::map(., ~log_reg(.)) %>%
  bind_rows() %>%
  mutate(term = ifelse(term == 'signs_YNfForaging Observed', 'Foraging Observed', term)) 
  
LR_zones_df %>%
  mutate(p.value = ifelse(p.value < 0.001, '< 0.001 **', 
                          ifelse(p.value < 0.05, paste(formatC(p.value,format = "f", 3), '*', sep = " "),
                                 ifelse(p.value < 0.1, paste(formatC(p.value,format = "f", 3), '.', sep = " "),
                                        formatC(p.value,format = "f", 3))))) %>%
  mutate_if(is.numeric, round, 2)%>%
  select(time_step, loss_gain, term, estimate,conf.low, conf.high, std.error, statistic, p.value) %>% 
  group_by(time_step, loss_gain) %>%
  gt() %>%
  tab_style(style = cell_fill('grey95'), locations = cells_row_groups()) %>%
  gtsave(., filename = normalizePath('NewPlots/log_reg_tabs.html', mustWork=FALSE))
  
summary(lr1)
