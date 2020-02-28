###Spatially explicit stage based population model used to evaluate the effects of ecological traps on metapopulations###
###Model is based on Wood Frog biology, but can be adapted to any species###
###Model has the ability to incorporate environmental stochasicity at the larval stage by manipulating hydroperiod for each population
###In V11 the dispersal matrix can be a distance matrix instead of a probability matrix

library(popbio)
library(ggplot2)
library(gridExtra)
library(msm)
library(data.table)
library(igraph)
library(pse)
library(dplyr)

setwd("N:/Population Model Results") #lab computer
setwd("Z:/cedge/Population Model Results") #laptop computer

#Population model
Amphib.hydro.1L.pop <- function(Pej1, Pj1j2, Pj2a, Paa, Fj2, Fa, k=5000, npop=4, area=NULL, suit=NULL, h1=NULL, h2=NULL, nYears, PerDisp = 0.05, DM=NULL, write = "F", filename){
  #Starting population matrix used to calculate starting populations
  s.matrix.1L <- matrix(c(0, 0, Fj2[1], Fa[1],
                          Pej1[1], 0, 0, 0,
                          0, Pj1j2[1], 0, 0,
                          0, 0, Pj2a[1], Paa[1]), nrow=4,ncol=4,byrow=T)
  
  lambda.1L <- lambda(s.matrix.1L) #Calculate Lambda for 1L
  StabStage.1L <- stable.stage(s.matrix.1L) #Calculate Stable Stage for 1L
  Elast.1L <- elasticity(s.matrix.1L) #Calculate Elasticity for 1L
  
  #calculating starting populations based on stable stage matrix
  A.1L <- 500/StabStage.1L[4] #based on a population of 500 adults
  InitAbund.1L <- c(A.1L * StabStage.1L[1], A.1L * StabStage.1L[2], A.1L * StabStage.1L[3], 500) #set intitial abundances to match stable stage
  
  
  #define area, suitability, and hydroperiod values, if they are not provided then use defualts
  if (is.null(area)){
    area.vals <- rep(10,npop)
  }else{
    area.vals <- area
  }
  
  if (is.null(suit)){
    suit.vals <- rep(1,npop)
  }else{
    suit.vals <- suit
  }
  
  if (is.null(h1)){
    h1.vals <- rep(0,npop)
  }else{
    h1.vals <- h1
  }
  
  if (is.null(h2)){
    h2.vals <- rep(0,npop)
  }else{
    h2.vals <- h2
  }
  
  #Creating the disperal matrix. If it is not specified then dispersal is uniform
  if (is.null(DM)){
    DispDist <- matrix(1, nrow=npop, ncol=npop)
    for (i in 1:npop){
      DispDist[i,i] <- 0
    }
  }else{
    DispDist <- DM
  }
  
  D1 <- apply(DispDist, 1:2, function(i) 2*exp(-i)) #Turn distances in to dispersal probability
  DispProb <- replace(D1, D1 == 2, 0) #Individuals don't disperse from patch i to patch i, and no dispersal between unconnected nodes
  DispPerc <- t(apply(DispProb, 1, function(i) (i/sum(i)))) #Turn dispersal probability into percent of individuals
  DispPerc <- replace(DispPerc, is.na(DispPerc), 0) #This is necessary if there is no dispersal to a population as NaN are returned in line above
  
  
  #Create the storage arrays for each population and do one time calculations
  k.vals <- rep(1,npop)
  allyears <- list()
  for (i in 1:npop){
    k.vals[i] <- k * area.vals[i] #carrying capacity
    
    allyears[[i]] <- matrix(0,nrow=nrow(s.matrix.1L)+1,ncol=nYears+1) #storage array for each population
    allyears[[i]][1:4,1] <- InitAbund.1L #set intitial abundance of each population to stable stage
    allyears[[i]][5,1] <- 2
  }
  
  #Creating a storage array for when wetlands dry, 0=before 1/2 year, 1=before end of year, 2=does not dry  
  Perm <- matrix(2, nrow=npop, ncol=nYears+1) 
  
  #Calculating populations in each year
  for (t in 2:(nYears+1)){
    #determine which wetland dry
    hydro <- rtnorm(1, 0.5, 0.20, 0, 0.9) #the hydroperiod this year
    
    #Calculate density dependance correction and the hydroperiod for each population
    DenDep <- rep(1,npop)
    pop.matrix <- list()
    emigrate <- rep(0,npop)
    dispersal <- matrix(0, nrow=npop, ncol=npop)
    for (i in 1:npop){
      DenDep[i] <- 1 - exp(-k.vals[i]/(allyears[[i]][3,t-1] + allyears[[i]][4,t-1]))
      Perm[i,t] <- ifelse(hydro < h1.vals[i], 0, ifelse(hydro < h2.vals[i], 1, 2))
      
      #Create a new matrix for each population      
      pop.matrix[[i]] <- matrix(
        c(
          0, 0, ifelse(allyears[[i]][3,t-1] < 5, 0, round(rtnorm(n=1, mean=Fj2[1], sd=Fj2[2], lower=0, upper=Inf))), ifelse(allyears[[i]][4,t-1] < 5, 0, round(rtnorm(n=1, mean=Fa[1], sd=Fa[2], lower=0, upper=Inf))),
          ifelse(Perm[i,t-1]>0, rtnorm(n = 1, mean = Pej1[1], sd = Pej1[2], lower = 0, upper = 1) * suit.vals[i], 0), 0, 0, 0,
          0, rtnorm(n = 1, mean = Pj1j2[1], sd = Pj1j2[2], lower = 0, upper = 1) * DenDep[i], 0, 0,
          0, 0, rtnorm(n = 1, mean = Pj2a[1], sd = Pj2a[2], lower = 0, upper = 1) * DenDep[i], rtnorm(n = 1, mean = Paa[1], sd = Paa[2], lower = 0, upper = 1) * DenDep[i]
        )
        ,nrow=4,ncol=4,byrow=T)
      
      allyears[[i]][1:4,t] <-  round(pop.matrix[[i]] %*% allyears[[i]][1:4,t-1]) #calculate the new population size, rounded to whole number
      
      #Number that immigrate to other populations stored as a matrix row disperses to column
      emigrate[i] <- round(allyears[[i]][3,t] * PerDisp) #total number of emmigrate from the population, rounded to whole number
      for (j in 1:npop){
        dispersal[i,j] <- round(emigrate[i]*DispPerc[i,j])
      }
      allyears[[i]][3,t] <- allyears[[i]][3,t] - emigrate[i]
      
      allyears[[i]][5,t] <- Perm[i,t] #add the hydroperiod for the year
      
    }#end of calculations for a population in a year
    
    #adds juvenile dispersers into population
    for (i in 1:npop){
      allyears[[i]][3,t] <- allyears[[i]][3,t] + colSums(dispersal)[i]
    }
  } # this is the end of the annual population calculations
  
  
  
  #transpose all the population data, renames the lists, and name the columns
  allyearstrans <- lapply(allyears, t)
  allyearstrans <- lapply(allyearstrans, as.data.frame)
  pops <- rep(paste("Pop",1:npop, sep=""))
  names(allyearstrans) <- pops
  allyearstrans <- Map(cbind, allyearstrans, new_clumn = names(allyearstrans))
  colna <- c("Eggs", "Juv1", "Juv2", "Adult", "Hydroperiod", "Pop")
  allyearstrans <- lapply(allyearstrans, setNames, colna)
  
  alldata <- as.data.frame(rbindlist(allyearstrans)) #combine all the lists for ploting and summary
  alldata$Year <- rep(1:(nYears+1), npop)
  alldata$Hydroperiod <- ifelse(alldata$Hydroperiod == 0, "Ephemeral", ifelse(alldata$Hydroperiod == 1,"Temporary", "Permanent"))
  alldata$Extinct <- ifelse(alldata$Eggs + alldata$Juv1 + alldata$Juv2 + alldata$Adult == 0, 1, 0) #if the population is extinct then value is 1
  
  #summarize all the data
  popsummary <- alldata %>% 
    group_by(Year) %>% 
    summarize(TotalEgg = sum(Eggs), TotalJuv1 = sum(Juv1), TotalJuv2 = sum(Juv2), TotalAdult = sum(Adult), PerExtinct = sum(Extinct)/npop)
  
  #Calculate stochastic growth rate, and determine when extinction occurs
  popsummary <- popsummary %>%
    mutate(r = ifelse(Year == 1 | TotalAdult == 0 | lag(TotalAdult) == 0, NA, log(TotalAdult/lag(TotalAdult))), Extinct = ifelse(TotalEgg + TotalJuv1 + TotalJuv2 + TotalAdult == 0, 1, 0))
  logr <- sum(popsummary$r, na.rm=T)*(1/length(which(!is.na(popsummary$r))))
  
  popsummary <- data.frame(popsummary)
  extinct <- ifelse(sum(popsummary$Extinct == 1) == 0, nYears + 1, popsummary[with(popsummary, order(-Extinct, Year)), ][1,1])

  #plot the data
  popplot <- ggplot(popsummary, aes(x=Year, y=TotalAdult))
  popplot <- (popplot +geom_point() + scale_y_log10() + geom_smooth())
  
  #Write .csv if required
  if (write == T){
    write.csv(alldata, file = filename) #Write a .csv with the population data over time
  }
  
  AP_Plot <- ggplot(alldata, aes(x=Year, y=(Adult+1), colour=factor(Pop))) #Plot adults over time
  #EP_Plot <- ggplot(Allpopulations, aes(x=Time, y=Egg, colour=factor(Population)))  #Plot eggs over time
  #TP_Plot <- ggplot(Allpopulations, aes(x=Time, y=Tadpole2nd, colour=factor(Population)))  #Plot Tadpole1 over time
  #JP1_Plot <- ggplot(Allpopulations, aes(x=Time, y=Juvenile1, colour=factor(Population)))
  #JP2_Plot <- ggplot(Allpopulations, aes(x=Time, y=Juvenile2, colour=factor(Population)))
  APlot <- (AP_Plot + geom_point()) # + theme(legend.position = "none"))
  #EPlot <- (EP_Plot + geom_point() + theme(legend.position = "none"))
  #TPlot <- (TP_Plot + geom_point() + theme(legend.position = "none"))
  #J1Plot <- (JP1_Plot + geom_point() + theme(legend.position = "none"))
  #J2Plot <- (JP2_Plot + geom_point() + theme(legend.position = "none"))
  
  list("Elast_1L" = Elast.1L, "StabStage_1L" = StabStage.1L, "Lambda_1L" = lambda.1L, "logr" = logr, "extinct" = extinct, "APlot" = APlot, "PPlot" = popplot, "output" = alldata, "summary" = popsummary)
} # this is the end of the function


###ONE SCENARIO EXAMPLE with a defined network
net1adj <- matrix(nrow = 4, ncol = 4, 0)
net1adj[1,2] <- 1
net1adj[1,3] <- 1
net1adj[1,4] <- 4
net1adj[2,1] <- 1
net1adj[2,3] <- 2
net1adj[3,1] <- 1
net1adj[3,2] <- 2
net1adj[4,1] <- 4

DM<-net1adj

plot(graph_from_adjacency_matrix(net1adj, mode='directed', weighted=TRUE))

#with variance on vital rates
d1 <- Amphib.hydro.1L.pop(Pej1=c(0.03515,0.02481), Pj1j2=c(0.3775,0.08847), Pj2a=c(0.2490,0.1535), Paa=c(0.08847,0.04871), Fj2=c(83.261,62.3), Fa=c(40.5645,13.4284), k=1000000000000000000000, nYears=100, npop=4, PerDisp = 1, DM=net1adj)
d1$APlot

#no variance on vital rates
d2 <- Amphib.hydro.1L.pop(Pej1=c(0.03515,0), Pj1j2=c(0.3775,0), Pj2a=c(0.2490,0), Paa=c(0.08847,0), Fj2=c(83.261,0), Fa=c(40.5645,0), k=1000000000000000000000, nYears=100, npop=4, PerDisp = 1, DM=net1adj)
d2$APlot

rtnorm(n=1, mean=0.3775, sd=0, lower=0.001, upper=Inf)



dset21 <- Amphib.hydro.1L.pop(Pej1=c(0.03515,0.02481), Pj1j2=c(0.3775,0.08847), Pj2a=c(0.2490,0.1535), Paa=c(0.08847,0.04871), Fj2=c(83.261,62.3), Fa=c(40.5645,13.4284), k=1000000000000000000000, npop=16, h1=rtnorm(n=16, mean=0.5, sd=0.2, lower=0.2, upper=0.8), PerDisp = 0, nYears=500, write = F, filename="dset21.csv")
dset21$APlot
dset21$summary
dset21$logr
dset21$extinct
dset21$output


###ONE SCENARIO REPLICATED EXAMPLE###

stime <- Sys.time()
simresults <- setNames(data.frame(matrix(nrow = 100, ncol = 8)), c("Network", "nTraps", "TSev", "Attract", "Disperse", "MaxA" , "metaR", "ExtinctTime"))
simresults[,1] <- "Full" #network type
simresults[,2] <- 5 #number of traps
simresults[,3] <- 10 #trap suit reduction
simresults[,4] <- 0 #trap attraction
simresults[,5] <- 0 #dispersal

sims <- vector("list", 100) #empty list of simulations
for (i in 1:100){
  net1 <- make_full_graph(50) #Full Network
  #net1 <- erdos.renyi.game(50,100, type="gnm",directed = F) #Erdos-Renyi
  #net1 <- make_tree(50, 2, mode = "undirected")
  net1adj <- as.matrix(as_adjacency_matrix(net1, sparse=T))
  
  #introduce traps
  traps <- rep(1,50)
  t.loc <- sample(1:50, 5)
  traps[t.loc] <- 0.9
  
  #run the simulation and extract values
  d1 <- Amphib.hydro.1L.pop(Pej1=c(0.03515,0.02481), Pj1j2=c(0.3775,0.08847), Pj2a=c(0.2490,0.1535), Paa=c(0.08847,0.04871), Fj2=c(83.261,62.3), Fa=c(40.5645,13.4284), k=1000000000000000000000, npop=50, suit=traps, h1=rtnorm(n=50, mean=0.4, sd=0.2, lower=0.2, upper=0.8), PerDisp = 0, nYears=500, DM=net1adj)
  simresults[i,6] <- max(d1$summary$TotalAdult)
  simresults[i,7] <- d1$logr
  simresults[i,8] <- d1$extinct
  sims[[i]] <- d1$summary #saves each of the simulations
}
etime <- Sys.time()
etime-stime

#simresults
mean(simresults$metaR)

write.csv(simresults, file="Full_10Random_10_0_0.csv")


###RUN A SET OF REPLICATED SIMULATIONS###

stime <- Sys.time()
simresults <- setNames(data.frame(matrix(nrow = 32400, ncol = 8)), c("Network", "nTraps", "TSev", "Attract", "Disperse", "MaxA" , "metaR", "ExtinctTime"))
#sims <- vector("list", 100) #empty list of simulations

simresults[,1] <- "Full" #network type

k <- 1 #counter variable
for (d in c(0, 0.1, 0.25, 0.5, 0.75, 0.99)){# d cycles through all the changes to trap attractiveness, or set the right amount
  for (b in c(0, 0.1, 0.25, 0.5, 0.75, 1)){ #b rotates through all the dispersal amounts, or set to the right amount
    for (a in c(1, 0.9, 0.75, 0.5, 0.25, 0)){ #a rotates through all fitness penalties
      for (j in c(0,5,12,25,37,50)){ #j rotates through all number of traps
        for (i in 1:25){
          simresults[k,2] <- j #number of traps
          simresults[k,3] <- (1-a)*100 #trap suit reduction
          simresults[k,4] <- d*100 #trap attraction
          simresults[k,5] <- b*100 #dispersal
          net1 <- make_full_graph(50) #Full Network
          #net1 <- erdos.renyi.game(50,100, type="gnm",directed = F) #Erdos-Renyi
          #net1 <- make_tree(50, 2, mode = "undirected")
          net1adj <- as.matrix(as_adjacency_matrix(net1, sparse=T))
      
          #introduce traps
          traps <- rep(1,50)
          t.loc <- sample(1:50, j)
          traps[t.loc] <- a
          
          #alter attraction of traps
          net1adj[, t.loc] <- apply(net1adj[, t.loc], 2, function(x) ifelse(x > 0, x - d, x))
      
          #run the simulation and extract values
          d1 <- Amphib.hydro.1L.pop(Pej1=c(0.03515,0.02481), Pj1j2=c(0.3775,0.08847), Pj2a=c(0.2490,0.1535), Paa=c(0.08847,0.04871), Fj2=c(83.261,62.3), Fa=c(40.5645,13.4284), k=1000000000000000000000, npop=50, suit=traps, h1=rtnorm(n=50, mean=0.4, sd=0.2, lower=0.2, upper=0.8), PerDisp = b, nYears=500, DM=net1adj)
          simresults[k,6] <- max(d1$summary$TotalAdult)
          simresults[k,7] <- d1$logr
          simresults[k,8] <- d1$extinct
          #sims[[i]] <- d1$summary #saves each of the simulations
          k <- k + 1
        }
      }
    }
  }
}
etime <- Sys.time()
etime-stime

#simresults
summarres <- simresults %>%
  group_by(nTraps, TSev, Disperse, Attract) %>%
  summarise(meanR = mean(metaR))
summarres <- as.data.frame(summarres)
summarres

write.csv(simresults, file="Full_Random_all_0_all.csv")


###LATIN SQUARE###
stime <- Sys.time() #start system timer

simresults <- setNames(data.frame(matrix(nrow = 10, ncol = 8)), c("Network", "nDisturbed", "Sev", "Attract", "Disperse", "MaxA" , "metaR", "ExtinctTime")) #Array to store results, nrow is the number of simulations to run
simresults[,1] <- "Erdos" #network type

LTsq <- LHS(factors= c("nDisturbed", "Sev", "Attract", "Dispersal"), N=10) #create simulation values, N is the number of simulations to run

for (i in 1:5000){#i is the number of simuations
  print(i)
  simresults[i,2] <- round(60 * LTsq$data[i,1]) #number of traps
  simresults[i,3] <- (1-LTsq$data[i,2])*100 #fitness reduction
  simresults[i,4] <- LTsq$data[i,3]*100 #trap attraction
  simresults[i,5] <- LTsq$data[i,4]*100 #dispersal
  #net1 <- make_full_graph(60) #Full Network
  net1 <- erdos.renyi.game(60,180, type="gnm",directed = F) #Erdos-Renyi
  #net1 <- make_tree(60, 2, mode = "undirected")
  net1adj <- as.matrix(as_adjacency_matrix(net1, sparse=T))
  
  #introduce disturbance
  Dist <- rep(1,60) #fitness in each patch is 1 by default
  
  #random selection of disturbed patches
  #d.loc <- sample(1:60, simresults[i,2]) #select disturbance location
  
  #top node degree used to select patches in Erdos network
  node.degree <- cbind(1:60,degree(net1, mode="total"), sample(1:60, 60)) #add a random number at the end to randomize equal degrees
  node.degree <- node.degree[order(-node.degree[,2], node.degree[,3]),]
  d.loc <- node.degree[1:simresults[i,2], 1]
  
  #top eccentricity used to selection patches for tree network
  #node.eccentricity <- cbind(1:60, eccentricity(net1), sample(1:60, 60)) #add a random number at the end to randomize equal degrees
  #node.eccentricity <- node.eccentricity[order(node.eccentricity[,2], node.eccentricity[,3]),]
  #d.loc <- node.eccentricity[1:simresults[i,2], 1]
  
  Dist[d.loc] <- LTsq$data[i,2] #set fitness in disturbed
  
  #alter attraction of traps
  net1adj[, d.loc] <- ifelse(net1adj[, d.loc] > 0, net1adj[, d.loc] * (1 - LTsq$data[i,3]), net1adj[, d.loc]) #changed from V6
  
  #run the simulation and extract values
  d1 <- Amphib.hydro.1L.pop(Pej1=c(0.03515,0.02481), Pj1j2=c(0.3775,0.08847), Pj2a=c(0.2490,0.1535), Paa=c(0.08847,0.04871), Fj2=c(83.261,62.3), Fa=c(40.5645,13.4284), k=1000000000000000000000, npop=60, suit=Dist, h1=rtnorm(n=60, mean=0.4, sd=0.2, lower=0.2, upper=0.8), PerDisp = LTsq$data[i,4], nYears=500, DM=net1adj)
  simresults[i,6] <- max(d1$summary$TotalAdult) #max  population size
  simresults[i,7] <- d1$logr #mean growth rate
  simresults[i,8] <- d1$extinct #time to extinction
}

etime <- Sys.time()
etime-stime
write.csv(simresults, file="Erdos_TopCon_LatinSQ.csv")


###IMPORTANCE OF Dispersal

stime <- Sys.time() #start system timer

simresults <- setNames(data.frame(matrix(nrow = 1000, ncol = 8)), c("Network", "nDisturbed", "Sev", "Attract", "Disperse", "MaxA" , "metaR", "ExtinctTime")) #Array to store results, nrow is the number of simulations to run

simresults[,1] <- "Erdos" #network type

#LTsq <- LHS(factors= c("nDisturbed", "Sev", "Attract", "Dispersal"), N=5000) #create simulation values, N is the number of simulations to run

for (i in 1:1000){#i is the number of simuations
  print(i)
  simresults[i,2] <- 0 #number of traps
  simresults[i,3] <- 0 #fitness reduction
  simresults[i,4] <- 0 #trap attraction
  simresults[i,5] <- sample(1:100, 1) #dispersal
  #net1 <- make_full_graph(60) #Full Network
  net1 <- erdos.renyi.game(60,180, type="gnm",directed = F) #Erdos-Renyi
  #net1 <- make_tree(60, 2, mode = "undirected")
  net1adj <- as.matrix(as_adjacency_matrix(net1, sparse=T))
  
  #introduce disturbance
  Dist <- rep(1,60) #fitness in each patch is 1 by default
  
  #random selection of disturbed patches
  d.loc <- sample(1:60, simresults[i,2]) #select disturbance location
  
  Dist[d.loc] <- simresults[i,3]/100 #set fitness in disturbed
  
  #alter attraction of traps
  net1adj[, d.loc] <- ifelse(net1adj[, d.loc] > 0, net1adj[, d.loc] - simresults[i,4]/100, net1adj[, d.loc]) #changed from V6
  
  #run the simulation and extract values
  d1 <- Amphib.hydro.1L.pop(Pej1=c(0.03515,0.02481), Pj1j2=c(0.3775,0.08847), Pj2a=c(0.2490,0.1535), Paa=c(0.08847,0.04871), Fj2=c(83.261,62.3), Fa=c(40.5645,13.4284), k=1000000000000000000000, npop=60, suit=Dist, h1=rtnorm(n=60, mean=0.4, sd=0.2, lower=0.2, upper=0.8), PerDisp = simresults[i,5]/100, nYears=500, DM=net1adj)
  simresults[i,6] <- max(d1$summary$TotalAdult) #max  population size
  simresults[i,7] <- d1$logr #mean growth rate
  simresults[i,8] <- d1$extinct #time to extinction
}


etime <- Sys.time()
etime-stime
write.csv(simresults, file="Erdos_Disturbance.csv")



###INTERACTION BETWEEN NUMBER OF DISTURBED PATCHES AND SEVERITY

stime <- Sys.time() #start system timer

simresults <- setNames(data.frame(matrix(nrow = 1000, ncol = 8)), c("Network", "nDisturbed", "Sev", "Attract", "Disperse", "MaxA" , "metaR", "ExtinctTime")) #Array to store results, nrow is the number of simulations to run

simresults[,1] <- "Erdos" #network type

LTsq <- LHS(factors= c("nDisturbed", "Sev"), N=1000) #create simulation values, N is the number of simulations to run

for (i in 1:1000){#i is the number of simuations
  print(i)
  simresults[i,2] <- round(60 * LTsq$data[i,1]) #number of traps
  simresults[i,3] <- (1-LTsq$data[i,2])*100 #fitness reduction
  simresults[i,4] <- 0 #trap attraction
  simresults[i,5] <- 25 #dispersal
  #net1 <- make_full_graph(60) #Full Network
  net1 <- erdos.renyi.game(60,180, type="gnm",directed = F) #Erdos-Renyi
  #net1 <- make_tree(60, 2, mode = "undirected")
  net1adj <- as.matrix(as_adjacency_matrix(net1, sparse=T))
  
  #introduce disturbance
  Dist <- rep(1,60) #fitness in each patch is 1 by default
  
  #random selection of disturbed patches
  d.loc <- sample(1:60, simresults[i,2]) #select disturbance location
  
  Dist[d.loc] <- LTsq$data[i,2] #set fitness in disturbed
  
  #alter attraction of traps
  net1adj[, d.loc] <- ifelse(net1adj[, d.loc] > 0, net1adj[, d.loc] - simresults[i,4]/100, net1adj[, d.loc]) #changed from V6
  
  #run the simulation and extract values
  d1 <- Amphib.hydro.1L.pop(Pej1=c(0.03515,0.02481), Pj1j2=c(0.3775,0.08847), Pj2a=c(0.2490,0.1535), Paa=c(0.08847,0.04871), Fj2=c(83.261,62.3), Fa=c(40.5645,13.4284), k=1000000000000000000000, npop=60, suit=Dist, h1=rtnorm(n=60, mean=0.4, sd=0.2, lower=0.2, upper=0.8), PerDisp = simresults[i,5]/100, nYears=500, DM=net1adj)
  simresults[i,6] <- max(d1$summary$TotalAdult) #max  population size
  simresults[i,7] <- d1$logr #mean growth rate
  simresults[i,8] <- d1$extinct #time to extinction
}

etime <- Sys.time()
etime-stime
write.csv(simresults, file="Erdos_ndistxsev.csv")
