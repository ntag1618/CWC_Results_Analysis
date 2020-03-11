function(df, direction){
  
  No_Beaver_df <- df %>% filter(signs_YN == 0)
  
  
  Yes_Beaver_df <- df %>% filter(signs_YN == 1)

  
  
  if (direction == 'gain'){
    
    All_zone_list <- list(nA = nrow(No_Beaver_df), 
                          sA = nrow(No_Beaver_df %>% filter(canopy_change > 0)),
                          
                          nB = nrow(Yes_Beaver_df),
                          sB = nrow(Yes_Beaver_df %>% filter(canopy_change > 0))
                          
                          ) 
    
  } else if (direction == 'loss') {
    
    All_zone_list <- list(nA = nrow(No_Beaver_df), 
                          sA = nrow(No_Beaver_df %>% filter(canopy_change < 0)),
                          
                          nB = nrow(Yes_Beaver_df),
                          sB = nrow(Yes_Beaver_df %>% filter(canopy_change < 0))
                          
                          ) 
    
  } else {
    stop(' direction argument must beeither gain or loss' )
  }
  
  
  ######### The Stan model as a string.##########
  model_string <- "
// Here we define the data we are going to pass into the model
data {
  // Number of trials
  int nA;
  int nB;


  // Number of successes
  int sA;
  int sB;

  
  //int n; // Number of trials     # single model version
  // int s;  // Number of successes # single model version
}
// Here we define what 'unknowns' aka parameters we have.
parameters {
  real<lower=0, upper=1> No_Beaver_Rate;
  real<lower=0, upper=1> Yes_Beaver_Rate;


}
// The generative model
model {
  No_Beaver_Rate ~ uniform(0, 1);
  Yes_Beaver_Rate ~ uniform(0, 1);

  
  sA ~ binomial(nA, No_Beaver_Rate);
  sB ~ binomial(nB, Yes_Beaver_Rate);

}
"


######## Run the Stan Model ################

Zone_samples <- stan(model_code = model_string, data = All_zone_list, iter = 4000)

return(Zone_samples)
}