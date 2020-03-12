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

zones_df = read_csv('./data/CWC_can_heights_df.csv')

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
                           
                            
                           #FOR TESTING!!!!!

                           ts_all = ts_all [sample(nrow(ts_all), 5000, replace = FALSE, prob = NULL),]


                           # --- All Section --------------------------------------------------------------------------
                           
                           no_beaver = ts_all$canopy_height[ts_all$signs_YN==0] # Subset losses by No beavers
                           beaver = ts_all$canopy_height[ts_all$signs_YN==1]    # Subset losses by Yes beavers
                           
                           # ----- Run Bayes T Test ------
                           bayes_out = bayes.ttest(no_beav = no_beaver, beav = beaver)
                           
                           # ----- Look at the results -----
                           best.post3 = bayes_out$BEST_results
                           
                           # bestcheck(best.post3) # check out the summary/validation plots...
                           
                           all_plot <- bayes_out$ggplot_obj
                           
                           p1 <- all_plot +
                             labs(x = 'Canopy Height (m)', y = "Point Density", main = ts_name) 
                           
                           ggsave(sprintf('sd_plots/%s_AllHeight_post.png', ts_name), plot = p1)
                           
                           # save posteriors to disk:
                           write.csv(best.post3, file=sprintf('save_posteriors/%s_AllHeight_Post.csv', ts_name), row.names = FALSE)
                           
                           
                           return(list('plot1' = p1,'ts_name'= ts_name, 'BEST_r1' = best.post3))
                         }
stopCluster(cl)

toc()

# function to return dataframe with timestep column

return_dir_df <- function(x, dir){
  timestep <- x[[2]]
  df <- as_tibble(x[[dir]])
  df$time_step <- timestep
  return(df)
}

# --------- height sd plot -------------
reshapedTotal <- All_The_Bayes %>%
  lapply(function(x) return_dir_df(x, 3)) %>%
  bind_rows() %>%
  select(-c(nu, mu1, mu2))%>%
  melt()%>%
  rename(Zone = variable) %>%
  mutate(CategoryB = ifelse(Zone=='sigma2', 'Foraging Present', 'Foraging Absent' )) %>%
  mutate(time_step = fct_relevel(time_step, "Sep18", "Mar18", "Jan18", "Sep17","Feb17","Dec16"))


height_sd_plot <- ggplot(data = reshapedTotal, aes(y= time_step ,x=value, fill=CategoryB, colour = CategoryB)) + 
  stat_halfeyeh(alpha=0.6, size=1.5) +
  scale_fill_manual("Zone", values=darkPalette, labels=c('Foraging Absent', 'Foraging Present')) +
  scale_color_manual("Zone",values=cbbPalette, labels=c('Foraging Absent', 'Foraging Present') )+
  # coord_cartesian(xlim = c(0, 0.7)) +
  labs(x = 'Canopy Height standard deviation (m)')+
  theme(axis.title.y=element_blank())
height_sd_plot
ggsave('sd_plots/height_sd_all.jpg', plot = Total_plot, dpi=300, width = 174, height = 150, units="mm")


# --------- mean reshape and plot -------------
reshapedTotal <- All_The_Bayes %>%
  lapply(function(x) return_dir_df(x, 3)) %>%
  bind_rows() %>%
  select(-c(nu, sigma1, sigma2))%>%
  melt()%>%
  rename(Zone = variable) %>%
  mutate(CategoryB = ifelse(Zone=='mu2', 'Foraging Present', 'Foraging Absent' )) %>%
  mutate(time_step = fct_relevel(time_step, "Sep18", "Mar18", "Jan18", "Sep17","Feb17","Dec16"))


height_mean_plot <- ggplot(data = reshapedTotal, aes(y= time_step ,x=value, fill=CategoryB, colour = CategoryB)) + 
  stat_halfeyeh(alpha=0.6, size=1.5) +
  scale_fill_manual("Zone", values=darkPalette, labels=c('Foraging Absent', 'Foraging Present')) +
  scale_color_manual("Zone",values=cbbPalette, labels=c('Foraging Absent', 'Foraging Present') )+
  # coord_cartesian(xlim = c(0, 0.7)) +
  labs(x = 'Mean Canopy Height (m)')+
  theme(axis.title.y=element_blank())
height_mean_plot
ggsave('sd_plots/height_mean_all.jpg', plot = Total_plot, dpi=300, width = 174, height = 150, units="mm")


