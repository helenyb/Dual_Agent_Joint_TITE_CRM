#joint tite pocrm simulations

library(doParallel)
library(mvtnorm)
library(rjags)
library(cubature)
library(pocrm)
registerDoParallel(cores=20)
nsims<-1000

source("JOINT_TITE_POCRM.R")




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

ndoses<-10

#how to split the activity across cycles for data generation
eff.pattern1<-rep(1,3)/3
# eff.pattern2<-c(1:3)/6
# eff.pattern3<-c(3:1)/6
eff.pattern<-1
ncycles<-3



#calculate parameters for data generation
for(pattern in 1){
  for(eff_scen in 1:6){
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


#prior hyper-parameters for part 1 
prior_hyp1<-c()
#tox
prior_hyp1[1]<- 0#mu_alphaT
prior_hyp1[2]<-1/1.34^2 #tau_alphaT

#activity
prior_hyp1[3]<-0 #mu_alphaE
prior_hyp1[4]<-1/1.34^2  #tau_alphaE



prior_hyp2<-c()
#tox
prior_hyp2[1]<- (-4)#mu_alphaT
prior_hyp2[2]<-1 #tau_alphaT
prior_hyp2[3]<- 2.5 #mu_betaT
prior_hyp2[4]<-4 #tau_betaT
#activity
prior_hyp2[5]<-(-4) #mu_alphaE
prior_hyp2[6]<-1/1.5 #tau_alphaE
prior_hyp2[7]<-2 #mu_betaE
prior_hyp2[8]<- 1/0.5 #tau_betaE
#relationship
prior_hyp2[9]<-0 #mu_phi
prior_hyp2[10]<-0.01 #tau_phi


#skeletons
skeletonE<-getprior(halfwidth = 0.05,target=0.3,nu=1,nlevel=10)
skeletonT<-getprior(halfwidth = 0.02,target=0.3,nu=2,nlevel=10)




#orderings 
#numbering corresponds to rows of dose.ind.matrix
#e.g order 1 (1,1)->(1,2)->(1,3) etc

order1<-c(1,3,5,7,9,2,4,6,8,10) 

order2<-c(1,2,3,4,5,6,7,8,9,10)

order3<-c(1,3,2,5,4,7,6,9,8,10)

order4<-c(1,3,2,4,5,7,6,8,9,10)

order5<-c(1,2,3,5,4,7,6,9,8,10)


#working models
wmE<-getwm(matrix(c(order1,order2,order3,order4,order5),nrow=5, byrow=T),skeletonE)
wmT<-getwm(matrix(c(order1,order2,order3,order4,order5),nrow=5, byrow=T),skeletonT)

#prior probs assigned to working models for part 1
prior.probsE<-rep(1/nrow(wmE),nrow(wmE))
prior.probsT<-rep(1/nrow(wmT),nrow(wmT))


#number of simulations
nsims<-1000 


#escalation only allowed in 1 direction at a time
ds_rule<-"ON"

#default ordering
def.order<-c(1:10)

for(eff.scen.index in 1:6){
  eff<-eff.scen.index
  for(tox.scen.index in 1:6){
    
    tox<-tox.scen.index
    
    
    
    assign(paste(c("JOINTTITEPOCRM.eff",eff,".",eff.pattern,"_tox",tox),collapse=""),
           foreach(i=1:nsims, combine = list) %dopar% {
             ##function
             
             JOINT.TITE.POCRM(seed=i,tru.E.pars = get(paste(c("eff_scen",eff,".",pattern,"pars"),collapse="")),
                                  tru.T.pars=get(paste(c("tox_scen",tox,".pars"),collapse="")),tru.corET=-0.5,
                                  co_size=3,ncohorts=20 ,target=0.3,
                                  ncycles=ncycles,dose.skipping.rule=ds_rule, 
                                  prior_vec1=prior_hyp1, prior_vec2=prior_hyp2,
                                  sufficient.information=T,sufficient.information.lim=29,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                                  safety.stopping.high.toosafe=T,initial.one.cycle=T,
                                  w1=0.33, w2=1.09, upper.tox=0.301,C_eff=0.2,C_tox=0.2,effbound=0.2,toxbound=0.3,
                                  backfill=F,TITE=T,dose.indices=dose.ind.mat,
                                  wmE=wmE,wmT=wmT,prior.o.T=prior.probsT,prior.o.E=prior.probsE,
                                  default.order=def.order,start.dose=4,pause=T,a.stop.bound = 0)
             
           }
           
           
    )
    
    
    save.image(paste(c("JOINT_TITE_POCRM_sims.RData"),collapse=""))
    print(timestamp())
  }
}










