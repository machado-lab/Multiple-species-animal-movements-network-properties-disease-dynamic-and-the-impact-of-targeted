# Multiple species animal movements: network properties, disease dynamic and theimpact of targeted control actions 

Supplementary material from the manuscript **Multiple species animal movements: network properties, disease dynamic and the impact of targeted control actions**.

> the following code is in fact in R language and is available to everyone, if you want to use please cite the manuscript preprint avaliable here [preprint avaliable here](https://arxiv.org/abs/2107.10108?context=q-bio).

## load the necessary packages
```
if(!require(igraph)){install.packages("igraph")};library(igraph)
if(!require(tidyverse)){install.packages("tidyverse")};library(tidyverse)
if(!require(SimInf)){install.packages("SimInf")};library(SimInf)
```
## Creates a function to run the model and the desfragmetation of the model 
the following code uses the tidyverse [tidyverse](https://github.com/tidyverse/tidyverse) package and uses the [SimInf](https://github.com/stewid/SimInf) framkework you can visit their repositories for more information about it.

```r
nico_sier_granjas_random_fun <- function(sim.number = NA, 
                                         rm_par, # Social network parameter to be tested 
                                         beta, #  tramission coefficient value
                                         epsilon, # incubation period if requiried 
                                         nrm, # Number of nodes to remove from the network 
                                         events, # programmed events accordign with the SimInf package
                                         infectados.ini, # Farm premises ID of the initial infected farms  
                                         granjas, # population of all farms each line represent one farm 
                                         Specie.to.select, # specie to initiate the infection
                                         granjas.especie,# dataframe where each line represent one farm and the also indicating the specie 
                                         tspan # period of the simulation 
) {
  
  # assing the object to the function 
  sim.number <- sim.number
  rm_par <-  rm_par
  nrm <-  nrm
  events <- events
  Specie.to.select 
  granjas <- granjas
  granjas.especie <- granjas.especie
  tspan <- tspan 
  
  #recovery the id fot the farms 
  swine.id2 <- granjas.especie %>% filter(Specie == "Swine") %>% pull(node)
  bovine.id2 <- granjas.especie %>% filter(Specie == "Bovine") %>% pull(node)
  sheep.id2 <-  granjas.especie %>% filter(Specie == "Sheep/Goat") %>% pull(node)
  
  
  swine.id <- swine.id2[!swine.id2 %in% unique(c(bovine.id2, sheep.id2))]
  bovine.id <- bovine.id2[!bovine.id2 %in% unique(c(swine.id2, sheep.id2))]
  sheep.id <- sheep.id2[!sheep.id2 %in% unique(c(swine.id2, bovine.id2))]
  
  
  #Select the farms accoridign  with the specie
  if (Specie.to.select == "All species") {
    infectados.ini.random <- granjas %>%
      sample_n(., 1000 ) %>% pull(network.id)
    
  } else if (Specie.to.select == "Bovine") {
    
    infectados.ini.random <- granjas %>%
      filter(network.id %in% bovine.id) %>%
      sample_n(., 1000 ) %>% pull(network.id)
    
  } else if (Specie.to.select == "Swine") {
    
    infectados.ini.random <- granjas %>%
      filter(network.id %in% swine.id) %>%
      sample_n(., 1000 ) %>% pull(network.id)
    
  } else if (Specie.to.select == "Sheep") {
    
    infectados.ini.random <- granjas %>%
      filter(network.id %in% sheep.id) %>%
      sample_n(., 1000 ) %>% pull(network.id)
    
  } else if (Specie.to.select == "multiespecies") {
    
    infectados.ini.random <- granjas %>%
      filter(!network.id %in% unique(c(bovine.id, swine.id, sheep.id)))%>%
      sample_n(., 1000 ) %>% pull(network.id)
    
  } 
  
  
  # creates the U0 object accoriding with the SimIinf package 
  
  u0 <- data.frame(S = granjas$pop_total,
                   I = rep(0, nrow(granjas)),
                   R = rep(0, nrow(granjas)))
  
  u0$I[infectados.ini.random] <- round(u0$S[infectados.ini.random] * p_anim_inf) # internal prevalence
  u0$I[infectados.ini.random][u0$I[infectados.ini.random] == 0] <- 1 # (no existen ej 0.5 animales tonces redondeamos a 1  )
  u0 <- u0 %>% mutate(S = S-I ) 
  

  # creates a backup
  u01<- u0
  beta <- beta
  
  # select the farms according to the SNA parameter
  
  if(rm_par == "pg"){
    step_rm <- page_rank[1:nrm]
  } else if(rm_par == "prdx"){
    step_rm <- paradox1[1:nrm]
  }  else if(rm_par == "no_control"){
    step_rm <- "xx"
  }else if(rm_par == "dg"){
    step_rm <- all_degree[1:nrm]
  } else if(rm_par == "rm"){
    step_rm <- random.nodes[1:nrm]
  } else if(rm_par == "btw"){
    step_rm <- btw[1:nrm]
  } else if(rm_par == "cls"){
    step_rm <- cls[1:nrm] } 
  
  # Intervened properties do not send or receive animals
  node.to.remove <- step_rm # nodes  to remove 
  events_remov <-  events %>% mutate( n= ifelse(node %in% node.to.remove, 0, n ),
                                      n = ifelse(dest %in% node.to.remove, 0, n))
  
  #ahora un modelo SI ----
  model1 <- SIR(u0 = u01,                  # initial status of farms 
                tspan = tspan,             # simulation periodf 
                beta = beta, # The transmission rate from susceptible to Infected
                gamma = 0,                 # The recovery rate from infected to recovered. 
                events = events_remov)     # scheduled movements whitin farms 
  
  # run the model 
  result <- run(model = model1, threads = 4)  # select the number of cores in your computer
  # counting infected farms per day
  tr <- trajectory(model = result)
  tr <- setDT(tr)
  tr <- tr[I > 0 ]
  tr[, Swine := ifelse(node %in% swine.id2, 1, 0)]
  tr[, bovine := ifelse(node %in% bovine.id2, 1, 0)]
  tr[, sheep := ifelse(node %in% sheep.id2, 1, 0)]
  tr <- tr[, inf_leng_nodes:=length(unique(node)), by=list(time, Swine, bovine, sheep)]
  myVector <- c("time","Swine", "bovine", "sheep", "inf_leng_nodes")
  tr <- tr[, ..myVector]
  tr <- unique(tr)
  tr <- as.data.frame(tr)   
  
  #add  model parameters
  tr$Parameter <- rm_par
  tr$removed_nodes <- nrm
  tr$beta <- beta
  tr$sim_number <- sim.number
  return(tr)  }


```
