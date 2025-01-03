
#mod of contour_xyz to remove dependence on the order of x and y 
contour2_xyz <- function( x , y , z , ... ) {
  df <- data.frame(x=x,y=y,z=z)
  df <- df[order(y,x),]
  
  ux <- unique(df[order(df$x),]$x) 
  uy <- unique(df[order(df$y),]$y) 
  
  m <- matrix( df$z , nrow=length(ux) , ncol=length(uy) )
  
  contour( ux , uy , m , ... )
}

#Function to use rethinking::link with a dataframe and output a dataframe
link_df <- function( fit , data  , ... ) {
  require(dplyr)
  
  #remove this later
  #fit <- quap_model_test
  #data <- d2
  
  
  #How many formulas are we working with?  
  num_formulas = length(fit@links)
  
  #tibble of formula variable names
  varnames<- tibble(names=rep("null",num_formulas))
  for(i in 1:num_formulas){
    varnames$names[[i]] <- fit@links[[i]][[1]]
  }
  
  
  
  #Explicitly index the data so we can use the index to collect results later
  data_indexed <<- as_tibble(data) %>%
    rowid_to_column("index") 
  
  
  #do the linking
  link_output <- as_tibble(link(fit,data=data,...))
  
  if(num_formulas==1){
    #create a conversion table
    convert_tib <- tibble(key=colnames(as_tibble(link_output)),index=data_indexed$index)
    
    #reshape the link output
    my_data <- list()
    for(i in varnames$names){
      tmp=as_tibble( link_output  ) %>%
        gather(1:ncol(link_output ),key="key",value ="linked_value") %>%
        inner_join(convert_tib,"key") %>%
        inner_join(data_indexed,"index") %>%    
        rename(!!i :=  linked_value)
      
      my_data[[i]] <- tmp
    }
    
    
    for(i in varnames$names[[1]]){
      joined_dat=as_tibble( my_data[[i]] ) %>%
        select(-i)}
    
    for(i in varnames$names){
      joined_dat <- bind_cols(joined_dat,select(my_data[[i]],i))
    }
    
    
  }
  
  if(num_formulas>1){ 
    #create a conversion table
    convert_tib <- tibble(key=colnames(as_tibble(link_output[[1]])),index=data_indexed$index)
    
    #reshape the link output
    my_data <- list()
    for(i in varnames$names){
      tmp=as_tibble( link_output[[i]] ) %>%
        gather(1:ncol(link_output[[1]]),key="key",value ="linked_value") %>%
        inner_join(convert_tib,"key") %>%
        left_join(data_indexed,"index") %>%    
        rename(!!i :=  linked_value)
      
      my_data[[i]] <- tmp
    }
    
    
    for(i in varnames$names[[1]]){
      joined_dat=as_tibble( my_data[[i]] ) %>%
        select(-i)}
    
    for(i in varnames$names){
      joined_dat <- bind_cols(joined_dat,select(my_data[[i]],i))
    }
    
  }
  #return the tibble
  select(joined_dat ,-key)
  
  
  
}



#Function to use rethinking::sim with a dataframe and output a dataframe
sim_df <- function( fit , data  , ... ) {
  require(dplyr)
  
  
  #extract the name of the data column
  data_name <- fit@formula[[1]][[2]]
  
  #Explicitly index the data so we can use the index to collect results later
  data_indexed <<- as_tibble(data) %>%
    rowid_to_column("index")
  
  #Remove the column with the observations if present
  if(paste0(data_name) %in% colnames(data_indexed)){
    data_indexed <-  select(data_indexed, -paste0(data_name))
  }
  
  
  #do the simulation
  sim_output <- as_tibble(sim(fit,data=data,...))
  
  #create a conversion table
  convert_tib <- tibble(key=colnames(sim_output),index=data_indexed$index)
  
  
  #reshape the simulation output
  sim_output_tibble <- as_tibble(sim_output)%>%
    gather(1:ncol(sim_output),key="key",value ="sampled_value") %>%
    inner_join(convert_tib,"key") %>%
    #    left_join(data_indexed,"index") %>%    
    inner_join(data_indexed,"index") %>%    
    
    rename(!!data_name :=  sampled_value)
   
  
  #return the tibble
  select(sim_output_tibble ,-key)
  
  
  
}



#Versions that require tidybayes

link_df_tb = function(fit, data,  n = NULL, seed = NULL, ... ) {
  if(  !("tidybayes" %in% installed.packages()) ){
    stop("Requires tidybayes package")
  } 
  
  
  require(tidybayes)
  
  if (!(inherits(fit, "map") || inherits(fit, "quap"))){
    stop("Can only run on fitted models from quap or map out of the rethinking package")
  }
    
  
  set.seed(seed)
  value = ".value"

   
  
  
  if (is.null(n)) { n = 1000  }
  
  # get the draws from the link-level predictors
  num_formulas = length(fit@links)
  
  #tibble of formula variable names
  varnames<- tibble(names=rep("null",num_formulas))
  for(i in 1:num_formulas){
    varnames$names[[i]] <- fit@links[[i]][[1]]
  }
  
  
  draws_list = rethinking::link(fit, data, n = n, flatten = FALSE , ... )
# use tidybayes to make long versions of the rethinking::link output 
# loop through all the variables that have values and give them their actual
# names. 
  for(i in 1:(length(draws_list))){
    draws <- add_draws(data, draws_list[[i]], value = value) %>%
      rename(index=".row",draw.number=".draw",!!varnames$names[i] := ".value")
    if(i==1){draws_tot <- draws}
    else{draws_tot <-left_join(draws_tot,draws)}
  }
  
  ungroup(draws_tot)
}

sim_df_tb = function(fit, data,  n = NULL, seed = NULL, ... ) {

  
  if(  !("tidybayes" %in% installed.packages()) ){
    stop("Requires tidybayes package")
    } 
  
  require(tidybayes)
  
  if (!(inherits(fit, "map") || inherits(fit, "quap"))){
    stop("Can only run on fitted fits from quap or map out of the rethinking package")
  }
  
  
  set.seed(seed)
  value = ".value"
  
  if (is.null(n)) { n = 1000  }
  
  #extract the name of the data column
  data_name <- fit@formula[[1]][[2]]
  
  #from tidybayes
  sims = rethinking::sim(fit, data, n = n)
  draws = add_draws(data, sims, value = paste0(data_name,".","sim"))%>%
    rename(index=".row",draw=".draw")
  
ungroup(draws)  
}
