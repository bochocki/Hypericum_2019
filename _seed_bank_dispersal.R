#################################################################
# Pedro F Quintana Ascencio
# reviewed 10/04/2018
#################################################################

#################################################################
### CORE PROGRAM FOR DISPERSAL AND SEED BANK SIMULATIONS 
#################################################################

#################################################################
## Part I: Read the Hypericum data from files
#################################################################
rm(list=ls())
#getwd()
#setwd("C:/Users/pquintan/Documents/Data/Hypericum/Analysis 2015")

dat <- read.table("hcdem_ABS_2015_b.txt", header=T)
#names(dat)

dat45 <- read.table("hc103burnedb.txt", header=T)

hc_data <- reshape(dat,varying=list(names(dat)[6:26],names(dat)[27:47],names(dat)[48:68],names(dat)[69:89],
                                    names(dat)[90:110],names(dat)[111:131],names(dat)[132:152]),direction="long",
                   v.names=c("surv_init","init_height","rep","stems","surv_fin","fin_height","rep_fin"),times=2014:1994)
head(hc_data)

hc_data45 <- reshape(dat45,varying=list(names(dat45)[6:13],names(dat45)[14:21],names(dat45)[22:29],names(dat45)[30:37],
                                        names(dat45)[38:45],names(dat45)[46:53],names(dat45)[54:61]),direction="long",
                     v.names=c("surv_init","init_height","rep","stems","surv_fin","fin_height","rep_fin"),times=2008:2001)


#################################################################
##  Part II: Merging files 
#################################################################

hc_data <- rbind(hc_data,hc_data45)
hc_data$gap[hc_data$site==103 & hc_data$gap==2] <- 5197
hc_data$gap[hc_data$site==103 & hc_data$gap==1] <- 5212

#################################################################
##  spliting Site 1 
#################################################################

hc_data$site[hc_data$site==1 & hc_data$gap >6] <-101
hc_data$site[hc_data$site==1 & hc_data$gap <6] <-102
hc_data <- subset(hc_data,site!=1)

#################################################################
##  Plots for log height
#################################################################

dd <- data.frame(x=hc_data$init_height,y=hc_data$fin_height)


#################################################################
#### calculating log height
#################################################################

hc_data$ln.hc <- log(hc_data$init_height)  
hc_data$ln.fin_hc <- log(hc_data$fin_height)  
hc_data$ln.hc2 <- hc_data$ln.hc^2

#################################################################
## Part V: Stage definitions
#################################################################

hc_data$nstage[hc_data$surv_init==1] <- "adult"
hc_data$nstage[hc_data$surv_init==0] <- "rip"
hc_data$nstage[hc_data$surv_init==9] <- "prevdead"
hc_data$nstage[hc_data$surv_init==2] <- "missing"
hc_data$nstage[hc_data$surv_init==5] <- "yearling"
hc_data$nstage[hc_data$surv_init==3] <- "adult"
hc_data$nstage[hc_data$stems == 1 & hc_data$surv_init==3 & hc_data$init_height < 20 ] <- "yearling"

hc_data$nstagef[hc_data$surv_fin==1] <- "adult"
hc_data$nstagef[hc_data$surv_fin==0] <- "rip"
hc_data$nstagef[hc_data$surv_fin==9] <- "prevdead"
hc_data$nstagef[hc_data$surv_fin==2] <- "missing"
hc_data$nstagef[hc_data$surv_fin==5] <- "yearling"
hc_data$nstagef[hc_data$surv_fin==3] <- "adult"
hc_data$nstagef[hc_data$stems == 1 & hc_data$surv_fin==3 & hc_data$fin_height < 20 ] <- "yearling"


#############################################
##### TSF ####
#############################################

site_info <- read.table("Balds_for_simulation_o1.txt", header=T) ####
dim.site <- dim(site_info)

itt <- array(0, c(92,54))
z <- seq(0,62,1)

for (si in 1:92) {  
  habitats <- site_info[si,8:dim.site[2]]
  habitats <- c(rep(0,4),as.numeric(habitats))
  
  if(sum(habitats)==0) {
    itt[si,] <- (1:length(habitats))+8 
    fire_yr<- -100}
  if(sum(habitats)>0) {
    fire_yr <- which(habitats==1)
    for (u in 1:length(fire_yr)){
      if (length(fire_yr)==1) {
        ittw <- 1:(fire_yr[u]-1)
        itt[si,] <- c( (ittw+8),( fire_yr[u]:length(habitats-1))-fire_yr[u] +1 ) }
      
      if (length(fire_yr) > 1) {
        if (u == 1) { ittw <- (1:(fire_yr[u]-1)+8) } #+8
        if (u > 1 & u < length(fire_yr) ) 
        {ittw <- c( (ittw),( fire_yr[u-1]:(fire_yr[u]-1))-(fire_yr[u-1]-1) )} 
        if (u > 1 & u == length(fire_yr) ) 
        { ittw <- c( (ittw),( fire_yr[u-1]:(fire_yr[u]-1))-(fire_yr[u-1]-1)  )
        itt[si,] <- c( ittw,( fire_yr[u]:length(habitats-1))-fire_yr[u] +1 ) }                       
      }
    }
  }
}
#print(t(itt[1:92,]))


  
  c_vec.a <- seq(min(hc_data$ln.hc[!is.na(hc_data$ln.hc) & hc_data$nstage=="adult"]),
                 max(hc_data$ln.hc[!is.na(hc_data$ln.hc) & hc_data$nstage=="adult"]),length.out=151) 
  ad.plants <- hist(hc_data$ln.hc[!is.na(hc_data$ln.hc) & hc_data$nstage=="adult"],breaks=c_vec.a,plot=FALSE)$counts/2
  c_vec.s <- seq(min(hc_data$ln.hc[!is.na(hc_data$ln.hc) & hc_data$nstage=="yearling"]),
                 max(hc_data$ln.hc[!is.na(hc_data$ln.hc) & hc_data$nstage=="yearling"]),length.out=151) 
  ad.yearli <- hist(hc_data$ln.hc[!is.na(hc_data$ln.hc) & hc_data$nstage=="yearling"],breaks=c_vec.s,plot=FALSE)$counts/2
  
  vec <- c(20*sum(ad.yearli,ad.plants),ad.yearli,ad.plants)
  sum(vec)
  sum(ad.yearli)
  sum(ad.plants)
  #length(vec)
  n_vec  <- matrix(rep(vec,92),c(301,92))
  r <- Nseeds <- N_plants <-N_adults<- N_yearlings<-array(0,c(dim.site[1],50+4))
  
  
  #################################################### 
  #####  Dispersal
  ####################################################
  mat.dis <- read.table("bald_dispersal_0.txt")
  # Matrices with dispersal probabilities 
  m.seeds <- array(0,c(92,92))
  #################################################### 
  #####  Seed bank
  ####################################################
  sbk <- 1 # 1, 0.8, 1.2
  hc_data$tsf <- hc_data$time - hc_data$burn_yr   
  z1 <- (z- mean(hc_data$tsf[!is.na(hc_data$tsf)]))/sd(hc_data$tsf[!is.na(hc_data$tsf)])
  
  #################################################### 
  #####  Start post simulations
  ####################################################
  nb <- 1 #92
  all.mat <- array(0, c(301,301,57,10))
  
  for(t in 1:54) {    
    for (si in 1:nb) { 

      if (t == 1){
      matname <- paste("matrices_",si,".RData",sep="")
      load(matname, .GlobalEnv,verbose=TRUE) 
      all.mat[,,,si] <- sim
      #####################################
      zz <- z1[itt[si,t]]
      S.dorm <- -0.1852*zz^2 + 0.0387*zz +0.471
      if (zz > 0.40){ S.dorm <- 0.38}
      sim[1,,itt[si,t]] <-(sim[1,,itt[si,t]]/S.dorm)*(S.dorm*sbk)
      sim[1,2:301,itt[si,t]] <- (sim[1,2:301,itt[si,t]]/(1-S.dorm))*(1-S.dorm*sbk)
      n_vec[,si] <- sim[,,itt[t]]%*%n_vec[,si]  
      Nseeds[si,t] <- n_vec[1,si]
      m.seeds[si,1:92] <- Nseeds[si,t]*mat.dis[,si]
      N_plants[si,t] <- sum(n_vec[2:301,si])
      N_adults[si,t] <- sum(n_vec[152:301,si])
      N_yearlings[si,t] <- sum(n_vec[2:151,si])
      #####################################
      } 
      if (t > 1){
        #####################################
        n_vec[1,si] <- sum(m.seeds[si,1:92])
        zz <- z1[itt[si,t]]
        S.dorm <- -0.1852*zz^2 + 0.0387*zz +0.471
        if (zz > 0.40){ S.dorm <- 0.38}
        if (si==2 & itt[si,t] > 57) {itt[si,t] <-57 }
        all.mat[1,,itt[si,t],si] <-(all.mat[1,,itt[si,t],si]/S.dorm)*(S.dorm*sbk)
        all.mat[2:150,1:301,itt[si,t],si] <- (all.mat[2:150,1:301,itt[si,t],si]/(1-S.dorm))*(1-S.dorm*sbk)
        ##################################### 
        n_vec[,si] <- all.mat[,,itt[si,t],si]%*%n_vec[,si]   	# project 1 year ahead,
        Nseeds[si,t] <- n_vec[1,si]   			# sum n to get new total pop. size,
        m.seeds[si,1:92] <- Nseeds[si,t]*mat.dis[,si]
    if (itt[si,t]==1) {
      n_vec[2:151,si] <- rep(0,150)
      if (sum(n_vec[152:301,si]) < 0){
        n_vec[152:301,si] <- rep(0,150)
      }
    }
    N_plants[si,t] <- sum(n_vec[2:301,si])
    N_adults[si,t] <- sum(n_vec[152:301,si])
    N_yearlings[si,t] <- sum(n_vec[2:151,si])
    #r[si,t]=log(N_plants[si,t])			      # calculate log growth rate, and
  } # si

  } # if > t
} # t



#NumCores <- as.numeric(Sys.getenv("SLURM_NPROCS"))
#JobNumber <- as.numeric(Sys.getenv("SLURM_JOB_ID"))

#NumCores <- as.numeric(Sys.getenv("SLURM_NPROCS"))
#JobNumber <- as.numeric(Sys.getenv("SLURM_JOB_ID"))

##put your desired output name file here
##put your desired output name file here
#Outname1 <- paste("Base for adults_32",ppp,".txt",sep="")
#Outname2 <- paste("Base for seeds_32",ppp,".txt",sep="")

## cat output name + directory
#OutDir1 <- paste(JobNumber, "/", Outname1, sep="")
#write.table(N_plants[ppp,], OutDir1 ,sep="\t")
#OutDir2 <- paste(JobNumber, "/", Outname2, sep="")
#write.table(Nseeds[ppp,], OutDir2 ,sep="\t")
#OutDir3 <- paste(JobNumber, "/", Outname3, sep="")
#save(sim, file = OutDir3 , ascii=T)


## cat output name + directory
  #write.table(N_plants[ppp,],Outname1 ,sep="\t")
  #write.table(Nseeds[ppp,], Outname2 ,sep="\t")
  t(N_plants[1:nb,])
  t(Nseeds[1:nb,])
  