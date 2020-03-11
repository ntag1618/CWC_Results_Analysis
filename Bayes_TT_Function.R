function(no_beav,beav){
  # ---- ggplot theme setting ------
  custPalette <- c( "#56B4E9", "#E69F00")
  theme_set(theme_bw() + 
              theme(legend.title=element_blank(),
                    legend.position = 'bottom',
                    legend.key.size = unit(1.2, "cm")))
  
  BESTout <- BESTmcmc(no_beav, beav, parallel=TRUE)
  
  p.sum.df <- as.data.frame(summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                                    ROPEeff=c(-0.2,0.2))) %>%
    rownames_to_column()
  
  # get posterior draws only 
  All_posterior <- BESTout[c(1,2)]
  
  
  ############### Plotting The posterior distribution - Density plots ##################################
  
  # ----- Reshape Data for plotting ---------------
  reshapedPost <- melt(All_posterior)
  names(reshapedPost)[names(reshapedPost) == "variable"] <- "Category"
  reshapedPost$CategoryB = 'Foraging Absent'
  reshapedPost$CategoryB[reshapedPost$Category=='mu2'] = 'Foraging Present'
  
  
  # ---- Generate ggplot object ----------------
  post_plot <- ggplot(data = reshapedPost, aes(x=value, fill=CategoryB)) + 
    geom_density(alpha=0.4,colour = "grey50", size = 0.8) +
    scale_fill_manual("Zone", values=custPalette) 
  
  # post_plot
  # ggsave('R_plots/Sep_Sep_Post.jpg')
  
  #---------- Return List -------------
  return_list = list('BEST_results' = BESTout, 'summary_tab' = p.sum.df, 'ggplot_obj' = post_plot)
  
  return(return_list)
  
  
}