source("JOINT_TITE_BLRM.R")


library(doParallel)
library(mvtnorm)
library(rjags)
library(cubature)
registerDoParallel(cores=20)

#define the labelling of doses
dose.ind.mat<-which(matrix(c(1:10),nrow=2)>0,arr.ind=T)

#scens

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


##picking scens to calibrate prior over
#T1.A3
tox_scenP1<-tox_scen1
eff_scenP1<-eff_scen3

#T2.A4
tox_scenP2<-tox_scen2
eff_scenP2<-eff_scen4

#T3.A5
tox_scenP3<-tox_scen3
eff_scenP3<-eff_scen5

#T3.A6
tox_scenP4<-tox_scen3
eff_scenP4<-eff_scen6

#T4.A1
tox_scenP5<-tox_scen4
eff_scenP5<-eff_scen1

#T4.A5
tox_scenP6<-tox_scen4
eff_scenP6<-eff_scen5

#T5.A2
tox_scenP7<-tox_scen5
eff_scenP7<-eff_scen2

#T6.A1
tox_scenP8<-tox_scen6
eff_scenP8<-eff_scen1

scens<-c()
for(sceni in 1:8){
  scens[sceni]<-paste(c("P",sceni),collapse = "")
}


ndoses<-10

#how to split the activity across cycles for data generation
eff.pattern1<-rep(1,3)/3
# eff.pattern2<-c(1:3)/6
# eff.pattern3<-c(3:1)/6
eff.pattern<-1
ncycles<-3



#calculate parameters for data generation
for(pattern in 1){
  for(eff_scen in scens){
    parMat<-matrix(nrow=2,ncol=ndoses)
    cycleMat<-     cyc_func_eff_v2(cyc_all_vec = get(paste(c("eff_scen",eff_scen),collapse="")), 
                                   split_vec = get(paste(c("eff.pattern",pattern),collapse="")))
    for(j in 1:ncol(cycleMat)){
      parMat[,j]<- find_lognormal_parms3(p1=cycleMat[1,j],p3=sum(cycleMat[,j]),int2=seq(0.01,10,0.01))
    }
    assign(paste(c("eff_scen",eff_scen,".pars"),collapse=""),
           parMat)
  }
}


for(tox in scens){
  parMat<-matrix(nrow=2,ncol=ndoses)
  cycleMat<-cyc_func_tox_v2(cyc_all_vec=get(paste(c("tox_scen",tox),collapse="")),cyc1_prop=0.75)
  for(j in 1:ncol(cycleMat)){
    parMat[,j]<-  find_lognormal_parms3(p1=cycleMat[1,j],p3=cycleMat[2,j],int2=seq(0.01,10,0.01))
  }
  assign(paste(c("tox_scen",tox,".pars"),collapse=""),
         parMat)
}


#prior hyper-parameters 
#prior_vecZ: prior hyperparameters for joint 2-parameter model (Z)
prior_vecZ<-c(
  , #mu_alphaT=prior_vecZ[1] prior mean for a_T (intercept, toxicity)
  , #tau_alphaT=prior_vecZ[2] prior precision for a_T (intercept, toxicity)
  , #mu_betaT=prior_vecZ[3] prior mean for b_T (slope, toxicity)
  ,  #tau_betaT=prior_vecZ[4] prior precision for b_T (slope, toxicity)
  ,  #mu_alphaE=prior_vecZ[5] prior mean for a_E (intercept, efficacy)
  , #tau_alphaE=prior_vecZ[6] prior precision for a_E (intercept, efficacy)
  ,  #mu_betaE=prior_vecZ[7] prior mean for b_E (slope, efficacy)
    #tau_betaE=prior_vecZ[8] prior precision for b_E (slope, efficacy)
)

#prior_vecW: prior hyperparameters for joint 2-parameter model (W)
prior_vecW<-c(
  ,  #mu_alphaT=prior_vecW[1] prior mean for a_T (intercept, toxicity)
  ,  #tau_alphaT=prior_vecW[2] prior precision for a_T (intercept, toxicity)
  , #mu_betaT=prior_vecW[3] prior mean for b_T (slope, toxicity)
  , #tau_betaT=prior_vecW[4] prior precision for b_T (slope, toxicity)
  ,  #mu_alphaE=prior_vecW[5] prior mean for a_E (intercept, efficacy)
  , #tau_alphaE=prior_vecW[6] prior precision for a_E (intercept, efficacy)
  , #mu_betaE=prior_vecW[7] prior mean for b_E (slope, efficacy)
  #tau_betaE=prior_vecW[8] prior precision for b_E (slope, efficacy)
)

#prior_vecI: prior hyperparameters for interactions
prior_vecI<-c(
  0,  #mu_phi=prior_vecI[1] prior mean for phi (Gumbel model)
  0.01,  #tau_phi=prior_vecI[2] prior precision for phi (Gumbel model)
  0,  #mu_eta_T=prior_vecI[3] prior mean for eta_T (W/Z interaction for toxicity)
  1,  #tau_eta_T=prior_vecI[4] prior precision for eta_T (W/Z interaction for toxicity)
  0, #mu_eta_E=prior_vecI[5] prior mean for eta_E (W/Z interaction for efficacy)
  1 #tau_eta_E=prior_vecI[6] prior precision for eta_E (W/Z interaction for efficacy)
)
#escalation only allowed in 1 direction at a time
ds_rule<-"ON"



#default ordering
def.order<-c(1:10)

#doses
dosesW<-c(600,1200)
dosesZ<-c(50,75,100,125,150)



#number of simulations
nsims<-1000 

#escalation only allowed in 1 direction at a time
ds_rule<-"ON"

#default ordering
def.order<-c(1:10)

####different TOX priors


##tox hyperparameters 1
prior_vecZ.1<-prior_vecZ
prior_vecZ.1[1:4]<-c(
, #mu_alphaT=prior_vecZ[1] prior mean for a_T (intercept, toxicity)
, #tau_alphaT=prior_vecZ[2] prior precision for a_T (intercept, toxicity)
, #mu_betaT=prior_vecZ[3] prior mean for b_T (slope, toxicity)
)  #tau_betaT=prior_vecZ[4] prior precision for b_T (slope, toxicity)
#
prior_vecW.1<-prior_vecW
prior_vecW.1[1:4]<-c(
  ,  #mu_alphaT=prior_vecW[1] prior mean for a_T (intercept, toxicity)
  ,  #tau_alphaT=prior_vecW[2] prior precision for a_T (intercept, toxicity)
  , #mu_betaT=prior_vecW[3] prior mean for b_T (slope, toxicity)
  ) #tau_betaT=prior_vecW[4] prior precision for b_T (slope, toxicity)

##tox hyperparameters 2
prior_vecZ.2<-prior_vecZ
prior_vecZ.2[1:4]<-c(
  , #mu_alphaT=prior_vecZ[1] prior mean for a_T (intercept, toxicity)
  , #tau_alphaT=prior_vecZ[2] prior precision for a_T (intercept, toxicity)
  , #mu_betaT=prior_vecZ[3] prior mean for b_T (slope, toxicity)
)  #tau_betaT=prior_vecZ[4] prior precision for b_T (slope, toxicity)
#
prior_vecW.2<-prior_vecW
prior_vecW.2[1:4]<-c(
  ,  #mu_alphaT=prior_vecW[1] prior mean for a_T (intercept, toxicity)
  ,  #tau_alphaT=prior_vecW[2] prior precision for a_T (intercept, toxicity)
  , #mu_betaT=prior_vecW[3] prior mean for b_T (slope, toxicity)
) #tau_betaT=prior_vecW[4] prior precision for b_T (slope, toxicity)

##tox hyperparameters 3
prior_vecZ.3<-prior_vecZ
prior_vecZ.3[1:4]<-c(
  , #mu_alphaT=prior_vecZ[1] prior mean for a_T (intercept, toxicity)
  , #tau_alphaT=prior_vecZ[2] prior precision for a_T (intercept, toxicity)
  , #mu_betaT=prior_vecZ[3] prior mean for b_T (slope, toxicity)
)  #tau_betaT=prior_vecZ[4] prior precision for b_T (slope, toxicity)
#
prior_vecW.3<-prior_vecW
prior_vecW.3[1:4]<-c(
  ,  #mu_alphaT=prior_vecW[1] prior mean for a_T (intercept, toxicity)
  ,  #tau_alphaT=prior_vecW[2] prior precision for a_T (intercept, toxicity)
  , #mu_betaT=prior_vecW[3] prior mean for b_T (slope, toxicity)
) #tau_betaT=prior_vecW[4] prior precision for b_T (slope, toxicity)


nsims<-1000



for(tox.prior.Z.par in c(1:3)){
  for(tox.prior.W.par in c(1:3)){ 
    for(scen in scens){   
      
      assign(paste(c("JOINTTITEBLRM.parT.W",tox.prior.W.par,".Z",tox.prior.Z.par,"_",scen),collapse=""),
             
                 foreach(i=1:nsims, combine = list) %dopar% {
               ##function
               
               JOINT.TITE.BLRM(seed=i,tru.E.pars = get(paste(c("eff_scen",scen,".pars"),collapse="")),
                               tru.T.pars=get(paste(c("tox_scen",scen,".pars"),collapse="")),tru.corET=-0.5,
                               co_size=3,ncohorts=20 ,target=0.3,
                               ncycles=ncycles,dose.skipping.rule=ds_rule, 
                               prior_vecZ=get(paste(c("prior_vecZ.",tox.prior.Z.par),collapse="")),
                               prior_vecW=get(paste(c("prior_vecW.",tox.prior.W.par),collapse="")), prior_vecI=prior_vecI,
                               sufficient.information=T,sufficient.information.lim=29,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                               safety.stopping.high.toosafe=T,initial.one.cycle=T,
                               w1=0.33, w2=1.09, upper.tox=0.301,C_eff=0.2,C_tox=0.2,effbound=0.2,toxbound=0.3,
                               backfill=F,TITE=T,dose.indices=dose.ind.mat,
                               dosesW=dosesW, dosesZ=dosesZ,
                               default.order=def.order,start.dose=4, pause=T,a.stop.bound = 0
               )
               
               
               
               
             }
             
             
      )
      
      
      
      
     
      
      save.image("JOINT_TITEBLRM_TCAL1.RData")
        }
  }
}
#########################################################################    


