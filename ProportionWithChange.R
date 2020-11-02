# # Proportion of Change in area
#
# # -------------------- load libraries ---------------
library(tidyverse)
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



# -------------- Summarise Data --------------------

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
  mutate(signs_YNf = fct_relevel(signs_YNf, 'No Foraging')) %>%
  filter(loss_gain != "NO_CHANGE")

Prop_change

ggplot(Prop_change, aes(x=loss_gain, y = Perc_Area, fill = signs_YNf)) +
  geom_bar(width=0.8, stat="identity", position=position_dodge(width=0.9), colour='black', alpha=0.8)+
  geom_text(aes(label=paste(round(Perc_Area, 1), "%")),position=position_dodge(width=0.95), vjust=1.5, size=3) +
  facet_wrap(~time_step) +
  labs(y='Area with Elevation Change (%)', x= 'Direction of Elevation Change') +
  scale_fill_manual("Zone", values=c('#377eb8', '#ff7f00')) 

ggsave('NewPlots/Amount_of_change.jpg', dpi=600, width=6, height = 6)
