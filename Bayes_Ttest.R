# ---- Load Packages -----
library(BEST)
library(tidyverse)
library(reshape2)
library(bayestestR)
library(ggrepel)
library(rstanarm)
library(tictoc)
library(foreach)
library(doParallel)
library(tidybayes)

# ---- ggplot theme setting ------



theme_set(theme_bw() + 
            theme(legend.title=element_blank(),
                  legend.position = 'bottom',
                  legend.key.size = unit(0.75, "cm"),
                  text = element_text(size=12)))

cbbPalette <- c("grey50", "grey10")
darkPalette <- c( "#56B4E9", "#E69F00")

# ------ Import and Clean Data ----------

zones_df = read_csv('./data/CWC_can_change_df.csv')

zones_df= zones_df[complete.cases(zones_df), ] # remove all rows with NAs

ts_names <- (unique(zones_df$time_step)) # available time steps
print(ts_names)

# ----- Import Bayes T test funtion ----
bayes.ttest = dget('Bayes_TT_Function.R')

# --- funtion to return bayes summary plots etc. ---
bestcheck <- function(best.obj){
  summary(best.obj)
  plot(best.obj)
  plot(best.obj, "sd")
  plotPostPred(best.obj)
  plotAll(best.obj, credMass=0.95, ROPEm=c(-0.1,0.1),
          ROPEeff=c(-0.1,0.1), compValm=0.5, showCurve=TRUE)
}

ts.list <- zones_df %>%
  group_by(time_step) %>%
  group_split()


# ----- setup parallel backend to use many processors -------
# cores=detectCores()
# cores = cores[1]-1
cores=7
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)


# ---- for loop to create loss, gain, and total change bayes t tests and plots ----


tic('TOTAL')
All_The_Bayes <- foreach(ts_all = ts.list, .combine=list, .multicombine=TRUE,
                         .packages = c('tidyverse', 'BEST', 'reshape2'))  %dopar% {
                           # for (ts_all in ts.list){ # only needed for debugging
                             set.seed(1234)
                           ts_name = as.character(ts_all[1,4])                         
                           
                           # --- Subset data ----
                           #FOR REAL
                           # ts_gain = ts_all [ts_all $canopy_change > 0 ,]     #  Subset gains
                           # ts_loss = ts_all [ts_all $canopy_change < 0 ,]    #  Subset losses
                           
                           #FOR TESTING!!!!!
                           ts_gain = ts_all [ts_all $canopy_change > 0 ,] [sample(nrow( ts_all [ts_all $canopy_change > 0 ,]), 1000, replace = FALSE, prob = NULL),]     #  Subset gains
                           ts_loss = ts_all [ts_all $canopy_change < 0 ,] [sample(nrow(ts_all [ts_all $canopy_change < 0 ,]), 1000, replace = FALSE, prob = NULL),]   #  Subset losses
                           ts_all = ts_all [sample(nrow(ts_all), 1000, replace = FALSE, prob = NULL),]
                           
                           # ---gains section -------------------------------------------------------------------------
                           
                           
                           no_beaver = ts_gain$canopy_change[ts_gain$signs_YN==0] # Subset Gains by No beavers
                           beaver = ts_gain$canopy_change[ts_gain$signs_YN==1]    # Subset Gains by Yes beavers
                           
                           # ----- Run Bayes T Test ------
                           bayes_out = bayes.ttest(no_beav = no_beaver, beav = beaver)
                           
                           # ----- Look at the results -----
                           best.post1 = bayes_out$BEST_results
                           
                           # bestcheck(best.post1) # check out the summary/validation plots...
                           
                           gain_plot <- bayes_out$ggplot_obj
                           
                           p1 <- gain_plot +
                             labs(x = bquote('Canopy Volume '~( m^3 )~ 'gain / '~m^2), y = "Point Density", main = ts_name)
                           
                           ggsave(sprintf('R_plots/%s_gain_post.png', ts_name), plot = p1)
                           
                           # save posteriors to disk:
                           write.csv(best.post1, file=sprintf('save_posteriors/%s_Gain_Post.csv', ts_name), row.names = FALSE)
                           
                           
                           # --- loss section --------------------------------------------------------------------------------
                           
                           no_beaver = abs(ts_loss$canopy_change[ts_loss$signs_YN==0]) # Subset losses by No beavers
                           beaver = abs(ts_loss$canopy_change[ts_loss$signs_YN==1])    # Subset losses by Yes beavers
                           
                           # ----- Run Bayes T Test ------
                           bayes_out = bayes.ttest(no_beav = no_beaver, beav = beaver)
                           
                           # ----- Look at the results -----
                           best.post2 = bayes_out$BEST_results
                           
                           # bestcheck(best.post2) # check out the summary/validation plots...
                           
                           loss_plot <- bayes_out$ggplot_obj
                           
                           p2 <- loss_plot +
                             labs(x = bquote('Canopy Volume '~( m^3 )~ 'loss / '~m^2), y = "Point Density", main = ts_name)
                           
                           ggsave(sprintf('R_plots/%s_loss_post.png', ts_name), plot = p2)
                           
                           # save posteriors to disk:
                           write.csv(best.post2, file=sprintf('save_posteriors/%s_Loss_Post.csv', ts_name), row.names = FALSE)
                           
                          
                           # --- All Section --------------------------------------------------------------------------
                           
                           no_beaver = ts_all$canopy_change[ts_all$signs_YN==0] # Subset losses by No beavers
                           beaver = ts_all$canopy_change[ts_all$signs_YN==1]    # Subset losses by Yes beavers
                           
                           # ----- Run Bayes T Test ------
                           bayes_out = bayes.ttest(no_beav = no_beaver, beav = beaver)
                           
                           # ----- Look at the results -----
                           best.post3 = bayes_out$BEST_results
                           
                           # bestcheck(best.post3) # check out the summary/validation plots...
                           
                           all_plot <- bayes_out$ggplot_obj
                           
                           p3 <- all_plot +
                             labs(x = bquote('Canopy Volume '~( m^3 )~ 'Change / '~m^2), y = "Point Density", main = ts_name) 
                           
                           ggsave(sprintf('R_plots/%s_Allchange_post.png', ts_name), plot = p3)
                           
                           # save posteriors to disk:
                           write.csv(best.post3, file=sprintf('save_posteriors/%s_AllChange_Post.csv', ts_name), row.names = FALSE)
                           
                           
                           return(list('plot1' = p1,'plot2' =  p2,'plot3' =  p3, ts_name, 'BEST_r1' = best.post1, 'BEST_r2' = best.post2, 'BEST_r3' = best.post3))
                         }
stopCluster(cl)

toc()

# function to return dataframe with timestep column

return_dir_df <- function(x, dir){
  timestep <- x[[4]]
  df <- as_tibble(x[[dir]])
  df$time_step <- timestep
  return(df)
}

# --------- gains reshape and plot -------------
reshapedGains <- All_The_Bayes %>%
  lapply(function(x) return_dir_df(x, 5)) %>%
  bind_rows() %>%
  select(-c(nu, sigma1, sigma2))%>%
  melt()%>%
  rename(Zone = variable) %>%
  mutate(CategoryB = ifelse(Zone=='mu2', 'Foraging Present', 'Foraging Absent' )) %>%
  mutate(time_step = fct_rev(time_step))


Gains_plot <- ggplot(data = reshapedGains, aes(y= time_step ,x=value, fill=CategoryB, color=CategoryB)) + 
  stat_halfeyeh(alpha=0.6, size=1.5) +
  scale_fill_manual("Zone", values=darkPalette, labels=c('Foraging Absent', 'Foraging Present')) +
  scale_color_manual("Zone",values=cbbPalette, labels=c('Foraging Absent', 'Foraging Present') )+
  coord_cartesian(xlim = c(0, 0.175)) +
  labs(x = bquote('Canopy Volume ('~m^3~ ') Gain / '~m^2))+
  theme(axis.title.y=element_blank())
Gains_plot
ggsave('R_plots/Gains_All_TS.jpg', plot = Gains_plot, dpi=300, width = 174, height = 150, units="mm")

# --------- loss reshape and plot -------------

reshapedLoss <- All_The_Bayes %>%
  lapply(function(x) return_dir_df(x, 6)) %>%
  bind_rows() %>%
  select(-c(nu, sigma1, sigma2))%>%
  melt()%>%
  rename(Zone = variable) %>%
  mutate(CategoryB = ifelse(Zone=='mu2', 'Foraging Present', 'Foraging Absent' )) %>%
  mutate(time_step = fct_rev(time_step))


Loss_plot <- ggplot(data = reshapedLoss, aes(y= time_step ,x=value, fill=CategoryB, color = CategoryB)) + 
  stat_halfeyeh(alpha=0.6, size=1.5) +
  scale_fill_manual("Zone", values=darkPalette, labels=c('Foraging Absent', 'Foraging Present')) +
  scale_color_manual("Zone",values=cbbPalette, labels=c('Foraging Absent', 'Foraging Present') )+
  # coord_cartesian(xlim = c(0, 0.7)) +
  labs(x = bquote('Canopy Volume ('~m^3~ ') Loss / '~m^2))+
  theme(axis.title.y=element_blank())
Loss_plot
ggsave('R_plots/Loss_All_TS.jpg', plot = Loss_plot, dpi=300, width = 174, height = 150, units="mm")

# --------- overall reshape and plot -------------
reshapedTotal <- All_The_Bayes %>%
  lapply(function(x) return_dir_df(x, 7)) %>%
  bind_rows() %>%
  select(-c(nu, sigma1, sigma2))%>%
  melt()%>%
  rename(Zone = variable) %>%
  mutate(CategoryB = ifelse(Zone=='mu2', 'Foraging Present', 'Foraging Absent' )) %>%
  mutate(time_step = fct_rev(time_step))


Total_plot <- ggplot(data = reshapedTotal, aes(y= time_step ,x=value, fill=CategoryB, colour = CategoryB)) + 
  stat_halfeyeh(alpha=0.6, size=1.5) +
  scale_fill_manual("Zone", values=darkPalette, labels=c('Foraging Absent', 'Foraging Present')) +
  scale_color_manual("Zone",values=cbbPalette, labels=c('Foraging Absent', 'Foraging Present') )+
  # coord_cartesian(xlim = c(0, 0.7)) +
  labs(x = bquote('Net Canopy Volume ('~m^3~ ') Change / '~m^2))+
  theme(axis.title.y=element_blank())
Total_plot
ggsave('R_plots/Total_All_TS.jpg', plot = Total_plot, dpi=300, width = 174, height = 150, units="mm")
