# Multiple species animal movements: network properties, disease dynamic and the impact of targeted control actions 

Supplementary material from the manuscript **Multiple species animal movements: network properties, disease dynamic and the impact of targeted control actions**.

> the following code is in fact in R language and is available to everyone, if you want to use please cite the [manuscript](https://veterinaryresearch.biomedcentral.com/articles/10.1186/s13567-022-01031-2#citeas) avaliable here and the [SimInf](https://github.com/stewid/SimInf) package .

## Load the necessary packages
```
if(!require(igraph)){install.packages("igraph")};library(igraph)
if(!require(tidyverse)){install.packages("tidyverse")};library(tidyverse)
if(!require(SimInf)){install.packages("SimInf")};library(SimInf)
```
## Creates a function to run the model and the implement control action by node removal from the network
The following code uses the [tidyverse](https://github.com/tidyverse/tidyverse) package and uses the Siminf framkework you can visit their repositories for more information about it.

```r
network_control_actions <- function(sim.number = NA, 
                                         rm_par,                #SNA parameter
                                         beta,                  #  tramission coefficient value
                                         nrm,                   # Number of nodes to remove from the network 
                                         events,                # programmed events accordign with the SimInf package
                                         infectados.ini,        # Farm premises ID of the initial infected farms  
                                         granjas,               # population of all farms each line represent one farm 
                                         Specie.to.select,      # specie to initiate the infection
                                         granjas.especie,       # list of farms with the specie description c("Bovine", "Swine" etc..) 
                                         tspan                  # period of the simulation 
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
  
  
  # Creates the U0 object accoriding with the SimIinf package 
  
  u0 <- data.frame(S = granjas$pop_total,
                   I = rep(0, nrow(granjas)),
                   R = rep(0, nrow(granjas)))
  
  u0$I[infectados.ini.random] <- round(u0$S[infectados.ini.random] * p_anim_inf) # internal prevalence
  u0$I[infectados.ini.random][u0$I[infectados.ini.random] == 0] <- 1 # (no existen ej 0.5 animales tonces redondeamos a 1  )
  u0 <- u0 %>% mutate(S = S-I ) 
  

  # creates a backup
  u01<- u0
  beta <- beta
  
  # Select the farms according to the SNA parameter
  
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
                beta = beta,               # The transmission rate from susceptible to Infected
                gamma = 0,                 # The recovery rate from infected to recovered. 
                events = events_remov)     # scheduled movements whitin farms 
  
  # Run the model 
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

## Set the parameters of simulation
 
 ```r
tspan <-as.numeric(min(events$time)):as.numeric(max(events$time))          # 0 to 1000 days of simulation  
beta <- 0.7                                                                # tramission coefficient 
p_anim_inf <- 0.1                                                          # within farm prevalence 
num.simulaciones <- 100                                                    # number of stochastic simulations
list_measures <- c( "no_control")                                          # parametes to be tested "pg" = pagerank, "dg" = degreee, 
                                                                           # "rm" = random nodes,"btw"= betweennes, 
                                                                           # "cls" = ClusterCoefficient, "no_control" for any control action 
```
## Run the model 


```r
    result <- network_control_actions(rm_par = "no_control",               # Social network parameter to be tested 
                                     beta= beta,                           # tramission coefficient
                                     nrm= 1000,                            # nodes to remove
                                     events =events,                       # echeduled events
                                     infectados.ini = infectados.ini,      # cinitial onfected farms 
                                     sim.number = 1,                       # id of the simulation 
                                     granjas = population,                 # dataframe with the farm population list 
                                     Specie.to.select = index.especie,     # specie to initiate the infection
                                     granjas.especie= granjas.especie,     # dataframe indicating the specie c( "Swine", "Sheep", "Bovine", "multiespecies")  
                                     tspan = tspan)                        # period of time of the simulation

```



## Authors

Nicolas Cespedes Cardenas [![ORCIDiD](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-7884-2353),
Abagael Sykes [![ORCIDiD](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-3751-1798),
Francisco PN Lopes,
and Gustavo Machado [![ORCIDiD](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-7552-6144)

## Citation
Cardenas, N.C., Sykes, A.L., Lopes, F.P.N. et al. Multiple species animal movements: network properties, disease dynamics and the impact of targeted control actions. Vet Res 53, 14 (2022). https://doi.org/10.1186/s13567-022-01031-2


