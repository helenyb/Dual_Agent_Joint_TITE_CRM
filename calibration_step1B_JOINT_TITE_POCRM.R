#2nd round of prior cals, fixed skeleton and changing prior for mu and tau
library(pocrm)
library(rjags)
library(cubature)
library(mvtnorm)
library(doParallel)
registerDoParallel(cores=20)

source("JOINT_TITE_POCRM_priorcal.R")


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

#options for skeletons
skeleton1<-getprior(halfwidth = 0.02,target=0.3,nu=1,nlevel=10)
skeleton2<-getprior(halfwidth = 0.1,target=0.3,nu=1,nlevel=10)
skeleton3<-getprior(halfwidth = 0.02,target=0.3,nu=5,nlevel=10)
skeleton4<-getprior(halfwidth = 0.1,target=0.3,nu=5,nlevel=10)
skeleton5<-getprior(halfwidth = 0.02,target=0.3,nu=8,nlevel=10)
skeleton6<-getprior(halfwidth = 0.1,target=0.3,nu=8,nlevel=10)

skeleton_mat<-matrix(c(skeleton1,skeleton2,skeleton3,skeleton4,skeleton5,skeleton6),byrow=T,nrow=6)

#how to split the activity across cycles for data generation (no activity used in priorcal, so this is set to 1,0,0)
eff.pattern1<-c(1,0,0)

ndoses<-10
ncycles<-3

nsims<-1000 #full

#calculate parameters for data generation (no activity in priorcal)
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


#prior hyper-parameters for power model are fixed
prior_hyp1<-c()
#tox
prior_hyp1[1]<- 0#mu_alphaT
prior_hyp1[2]<-1/1.34^2 #tau_alphaT



#define starting dose (row of dose.ind.mat)
starting.dose<-4


#options for hyper-parameters
mu_alphaT_priors<-c(-5,-4)
tau_alphaT_priors<-c(1,2,4)
mu_betaT_priors<-c(2,2.5)
tau_betaT_priors<-c(1,2,4)

#orderings 
#numbering corresponds to rows of dose.ind.matrix
#e.g order 1 (1,1)->(1,2)->(1,3) etc

order1<-c(1,3,5,7,9,2,4,6,8,10) 

order2<-c(1,2,3,4,5,6,7,8,9,10)

order3<-c(1,3,2,5,4,7,6,9,8,10)

order4<-c(1,3,2,4,5,7,6,8,9,10)

order5<-c(1,2,3,5,4,7,6,9,8,10)


C_tox_val<-0.2

for (tau_alphaT_i in 1:length(tau_alphaT_priors)){
  for (tau_betaT_i in 1:length(tau_betaT_priors)){
    for (mu_alphaT_i in 1:length(mu_alphaT_priors)){
      for (mu_betaT_i in 1:length(mu_betaT_priors)){
        #prior hyper-parameters
        prior_hyp2<-c()
        #tox
        prior_hyp2[1]<-mu_alphaT_priors[mu_alphaT_i] #mu_alphaT
        prior_hyp2[2]<-tau_alphaT_priors[tau_alphaT_i] #tau_alphaT
        prior_hyp2[3]<- mu_betaT_priors[mu_betaT_i] #mu_betaT
        prior_hyp2[4]<-tau_betaT_priors[tau_betaT_i] #tau_betaT
        
        #fixed skeleton
        wmT<-getwm(matrix(c(order1,order2,order3,order4,order5),nrow=5, byrow=T),skeleton_mat[3,])
        
        prior.probsT<-rep(1/nrow(wmT),nrow(wmT))
        
        
        
        for(scen.index in 1:length(prior_scen_list)){
          tox<-scen.index
          # browser()
          assign(paste(c("JOINTTITEPOCRM.prior2_sk",3,"_alpha_mu",mu_alphaT_i,"_alpha_tau",tau_alphaT_i,"_beta_mu",mu_betaT_i,"_beta_tau",tau_betaT_i,"_scen",tox,"_start4"),collapse=""),
                 foreach(i=1:nsims, combine = list) %dopar% { #full (clusters)
                   #foreach(i=1:nsims, combine = list) %do% { #practice (windows)
                   ##function
                   
                   
                   JOINT.TITE.POCRM.priorcal(seed=i,tru.E.pars=eff_scen1.1pars,tru.T.pars=get(paste(c("prior_scen",tox,".pars"),collapse="")),tru.corET=-0.5,
                                       co_size=3,ncohorts=20 ,target=0.3,
                                       ncycles=ncycles,dose.skipping.rule="ON", 
                                       prior_vec1=prior_hyp1, prior_vec2=prior_hyp2,
                                       sufficient.information=T,sufficient.information.lim=29,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                                       safety.stopping.high.toosafe=T,initial.one.cycle=T,
                                       C_tox=C_tox_val,toxbound=0.3,
                                       backfill=F,backfill.num=2,TITE=T,dose.indices=dose.ind.mat,
                                       wmT=wmT,prior.o.T=prior.probsT,default.order=c(1:10),start.dose=4)
                 })#for assign
          
          
          save.image(paste(c("JOINT_TITE_POCRM_priorcal_st4_2.RData"),collapse=""))
        }#for scen
        
      } #for mu_beta
    } #for mu_alpha
  } #for tau_beta
} #for tau_alpha