#Joint TITE-BLRM prior calibration stage 1
source("JOINT_TITE_BLRM_priorcal.R")


library(rjags)
library(cubature)
library(mvtnorm)
library(doParallel)
registerDoParallel(cores=20)
##setting up scenarios


#mapping for generation of scenarios
scen.order1<-c(1,6,2,7,3,8,4,9,5,10)
scen.order2<-c(1:10)
scen.order3<-c(1,3,2,5,4,7,6,9,8,10)
scen.order4<-c(1,3,2,4,5,7,6,8,9,10)
scen.order5<-c(1,2,3,5,4,6,7,9,8,10)
##for calibration only:
scen.order6<-c(1,4,2,5,3,7,6,8,9,10)

#dose indices
dose.ind.mat<-which(matrix(c(1:10),nrow=2)>0,arr.ind=T)

##sequences
seq1<-seq(from=0.03,by=0.03,length.out=10) #10
seq2<-seq(from=0.02,by=0.04,length.out=10) #8
seq3<-seq(from=0.05,by=0.05,length.out=10) #6
seq4<-c(seq(from=0.06,by=0.06,length.out=8),seq(from=0.5,by=0.1,length.out=2)) #5
seq5<-c(0.2,seq(from=0.3,by=0.05,length.out=9)) #2
seq6<-seq(from=0.3,by=0.05,length.out=10) #1
seq7<-c(seq(from=0.4,by=0.05,length.out=6),rep(0.7,4)) #none

#scenarios
prior_scen1<-matrix(seq1[scen.order1],ncol=5)

prior_scen2<-matrix(seq2[scen.order4],ncol=5)

prior_scen3<-matrix(seq3[scen.order6],ncol=5)#

prior_scen4<-matrix(seq4[scen.order3],ncol=5)

prior_scen5<-matrix(seq4[scen.order2],ncol=5)

prior_scen6<-matrix(seq5[scen.order5],ncol=5)

prior_scen7<-matrix(seq6[scen.order3],ncol=5)#

prior_scen8<-matrix(seq7[scen.order4],ncol=5)#

#list of tox scenarios
prior_scen_list<-lapply(c(1:8),function(x) get(paste(c("prior_scen",x),collapse="")))

#eff not used in calibration, so fix at 1
eff_scen1<-matrix(rep(1,10),ncol=5)#

#how to split the activity across cycles for data generation
eff.pattern1<-c(1,0,0)
# eff.pattern2<-c(1:3)/6
# eff.pattern3<-c(3:1)/6

ndoses<-10
ncycles<-3


#calculate parameters for data generation
for(pattern in 1){
  for(eff_scen in 1){
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


for(tox in 1:8){
  parMat<-matrix(nrow=2,ncol=ndoses)
  cycleMat<-cyc_func_tox_v2(cyc_all_vec=get(paste(c("prior_scen",tox),collapse="")),cyc1_prop=0.75)
  for(j in 1:ncol(cycleMat)){
    parMat[,j]<-  find_lognormal_parms3(p1=cycleMat[1,j],p3=cycleMat[2,j],int2=seq(0.01,10,0.01))
  }
  assign(paste(c("prior_scen",tox,".pars"),collapse=""),
         parMat)
}


##dosea

dosesW<-c(600,1200)
dosesZ<-c(50,75,100,125,150)


#options for hyper-parameters
mu_alphaT_Z_priors<-c(-6.5,-6,-5.5)
mu_alphaT_W_priors<-c(-6.5,-6,-5.5)
tau_alphaT_priors<-c(2)
mu_betaT_Z_priors<-c(0,0.5,1)
mu_betaT_W_priors<-c(0,0.25,0.5)
tau_betaT_priors<-c(2)


prior_vecI_1<-c(0,1)

tau_alphaT_i<-1
tau_betaT_i<-1


nsims<-1000

for (mu_alphaT_W_i in 1:length(mu_alphaT_W_priors)){
  for (mu_betaT_W_i in 1:length(mu_betaT_W_priors)){
    for (mu_alphaT_Z_i in 1:length(mu_alphaT_Z_priors)){
      for (mu_betaT_Z_i in 1:length(mu_betaT_Z_priors)){
        #prior hyper-parameters
        prior_vecZ_1<-c()
        prior_vecZ_1[1]<-mu_alphaT_Z_priors[mu_alphaT_Z_i] #mu_alphaT
        prior_vecZ_1[2]<-tau_alphaT_priors[tau_alphaT_i] #tau_alphaT
        prior_vecZ_1[3]<- mu_betaT_Z_priors[mu_betaT_Z_i] #mu_betaT
        prior_vecZ_1[4]<-tau_betaT_priors[tau_betaT_i] #tau_betaT
        
        prior_vecW_1<-c()
        prior_vecW_1[1]<-mu_alphaT_W_priors[mu_alphaT_W_i] #mu_alphaT
        prior_vecW_1[2]<-tau_alphaT_priors[tau_alphaT_i] #tau_alphaT
        prior_vecW_1[3]<- mu_betaT_W_priors[mu_betaT_W_i] #mu_betaT
        prior_vecW_1[4]<-tau_betaT_priors[tau_betaT_i] #tau_betaT
        
        
        
        for(scen.index in 1:length(prior_scen_list)){
          tox<-scen.index
          # browser()
          assign(paste(c("JOINTTITEBLRM.prior_alpha_mu_Z",mu_alphaT_Z_i,"_beta_mu_Z",mu_betaT_Z_i,"_alpha_mu_W",mu_alphaT_W_i,"_beta_mu_W",mu_betaT_W_i,"_scen",tox,"_start4"),collapse=""),
                 foreach(i=1:nsims, combine = list) %dopar% { #full (clusters)
                   #   foreach(i=1:nsims, combine = list) %do% { #practice (windows)
                   ##function
                   
                   
                   JOINT.TITE.BLRM.priorcal(seed=i,tru.E.pars=eff_scen1.1pars,tru.T.pars=get(paste(c("prior_scen",tox,".pars"),collapse="")),tru.corET=-0.5,
                                      co_size=3,ncohorts=20 ,target=0.3,
                                      ncycles=ncycles,dose.skipping.rule="ON", 
                                      prior_vecZ=prior_vecZ_1, 
                                      prior_vecW=prior_vecW_1, 
                                      prior_vecI=prior_vecI_1,
                                      sufficient.information=T,sufficient.information.lim=29,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                                      safety.stopping.high.toosafe=T,initial.one.cycle=T,
                                      C_tox=0.2,toxbound=0.3,
                                      backfill=F,backfill.num=2,TITE=T,dose.indices=dose.ind.mat,
                                      dosesW=dosesW, dosesZ=dosesZ, default.order=c(1:10),start.dose=4,gs.iter=10000)
                   
                   
                 })#for assign
          
          
          save.image(paste(c("JOINT_TITE_BLRM_priorcal1.RData"),collapse=""))
        }#for scen
        
      } #for mu beta Z
    } #for mu alpha Z
  } #for mu beta W
} #for mu alpha W



