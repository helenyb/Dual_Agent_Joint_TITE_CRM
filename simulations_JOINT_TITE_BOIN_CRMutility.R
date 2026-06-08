#joint tite boin simulations
library(doParallel)
library(mvtnorm)
registerDoParallel(cores=20)
nsims<-1000

source("JOINT_TITE_BOIN_CRMutility.R")


#define the labelling of doses
dose.ind.mat<-which(matrix(c(1:10),nrow=2)>0,arr.ind=T)
#       row col
# [1,]   1   1
# [2,]   2   1
# [3,]   1   2
# [4,]   2   2
# [5,]   1   3
# [6,]   2   3
# [7,]   1   4
# [8,]   2   4
# [9,]   1   5
# [10,]  2   5

#definition of scenarios (prob for ALL CYCLES)
#ordering corresponds to dose.ind.matrix
tox_scen1<-c(0.03,0.05,0.07,0.09,0.11,0.13,0.15,0.25,0.2,0.3)
tox_scen2<-c(0.1,0.45,0.15,0.5,0.2,0.55,0.3,0.6,0.4,0.6)
tox_scen3<-c(0.05,0.1,0.08,0.12,0.15,0.3,0.2,0.4,0.45,0.5)
tox_scen4<-c(0.1,0.3,0.2,0.45,0.4,0.55,0.5,0.6,0.6,0.6)
tox_scen5<-c(0.3,0.4,0.45,0.5,0.5,0.55,0.55,0.6,0.6,0.6)
tox_scen6<-c(0.4,0.4,0.4,0.4,0.5,0.5,0.5,0.5,0.6,0.6)

eff_scen1<-c(0.2,0.25,0.3,0.4,0.35,0.45,0.5,0.6,0.55,0.65)
eff_scen2<-c(0.3,0.34,0.32,0.36,0.38,0.42,0.4,0.44,0.46,0.48)
eff_scen3<-c(0.06,0.1,0.08,0.15,0.12,0.25,0.2,0.35,0.3,0.4)
eff_scen4<-c(0.05,0.1,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55)
eff_scen5<-c(0.1,0.2,0.12,0.3,0.14,0.4,0.16,0.5,0.18,0.6)
eff_scen6<-c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2)
#plateau
eff_scen7<-c(0.2,0.25,0.3,0.35,0.4,0.45,0.4,0.45,0.4,0.45)
#bell
eff_scen8<-c(0.2,0.3,0.3,0.35,0.4,0.45,0.3,0.35,0.2,0.25)

ndoses<-10

#how to split the activity across cycles for data generation
eff.pattern1<-rep(1,3)/3
# eff.pattern2<-c(1:3)/6
# eff.pattern3<-c(3:1)/6
eff.pattern<-1
ncycles<-3



#calculate parameters for data generation
for(pattern in 1){
  for(eff_scen in 1:8){
    parMat<-matrix(nrow=2,ncol=ndoses)
    cycleMat<-     cyc_func_eff_v2(cyc_all_vec = get(paste(c("eff_scen",eff_scen),collapse="")), 
                                   split_vec = get(paste(c("eff.pattern",pattern),collapse="")))
    for(j in 1:ncol(cycleMat)){
      parMat[,j]<- find_lognormal_parms3(p1=cycleMat[1,j],p3=sum(cycleMat[,j]),int2=seq(0.01,10,0.01))
    }
    assign(paste(c("eff_scen",eff_scen,".",pattern,"pars"),collapse=""),
           parMat)
  }
}


for(tox in 1:6){
  parMat<-matrix(nrow=2,ncol=ndoses)
  cycleMat<-cyc_func_tox_v2(cyc_all_vec=get(paste(c("tox_scen",tox),collapse="")),cyc1_prop=0.75)
  for(j in 1:ncol(cycleMat)){
    parMat[,j]<-  find_lognormal_parms3(p1=cycleMat[1,j],p3=cycleMat[2,j],int2=seq(0.01,10,0.01))
  }
  assign(paste(c("tox_scen",tox,".pars"),collapse=""),
         parMat)
}





#number of simulations
nsims<-1000 



#default ordering
def.order<-c(1:10)


##boin utility definition
u00<-40 #no tox no act
u01<-100 #no tox yes act
u10<-0 #yes tox no act
u11<-60 #yes tox yes act
utility_mat<-matrix(c(u00,u10,u01,u11),nrow=2)


for(eff.scen.index in 1:8){
  eff<-eff.scen.index
  for(tox.scen.index in 1:6){
    
    tox<-tox.scen.index
    
    assign(paste(c("JOINTTITEBOIN.eff",eff,".",eff.pattern,"_tox",tox,"st4"),collapse=""),
          # foreach(i=1:nsims, combine = list) %dopar% {
              foreach(i=1:nsims, combine = list) %do% {
             ##function
             
             JOINT.TITE.BOIN_CRMutility(seed=i,tru.E.pars = get(paste(c("eff_scen",eff,".",pattern,"pars"),collapse="")),
                             tru.T.pars=get(paste(c("tox_scen",tox,".pars"),collapse="")),tru.corET=-0.5,
                             co_size=3,ncohorts=20 ,target=0.3, targetE=0.2,
                             ncycles=ncycles,
                             sufficient.information=T,sufficient.information.lim=29,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                             safety.stopping.high.toosafe=T,initial.one.cycle=T,
                             C_eff=0.2,C_tox=0.2,effbound=0.2,toxbound=0.3,
                             backfill=F,TITE=T,dose.indices=dose.ind.mat,
                             default.order=def.order,start.dose=4,utility_mat =utility_mat,pause=T, a.stop.bound = 61 )
             
           }
           
           
    )
    
    save.image(paste(c("JOINT_TITE_BOIN_sims.RData"),collapse=""))
    print(timestamp())
    
  }
}






