##pocrm prior calibration results extractions stage 2
load("JOINT_TITE_POCRM_priorcal_st4_2.RData")



##correct dose pairs
##prior calibrations
MTD_function<-function(truth,target){
  if(all(truth>target+1e-20)){
    return("UNSAFE")
  }else if(all(truth<target-1e-20)){
    return("TOOSAFE")
  }else{
    difs<-target-truth
    difs[difs<(-1e-10)]<-10
    MTD_out<-max(which(difs==min(difs)))  
    return(MTD_out)
  }
}

correct_list<-unlist(lapply(prior_scen_list, function(x) MTD_function(x,0.3)))


#correct_in<-unlist(lapply(prior_scen_list, function(x) which(abs(x-0.3)<1e-10)))

#empty results array
result_array<-array(NA,dim=c(length(prior_scen_list),
                             length(mu_alphaT_priors),
                             length(mu_betaT_priors),
                             length(tau_alphaT_priors),
                             length(tau_betaT_priors) 
))      


skeleton_i<-3

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
        
        wmT<-getwm(matrix(c(order1,order2,order3,order4,order5),nrow=5, byrow=T),skeleton_mat[skeleton_i,])
        
        prior.probsT<-rep(1/nrow(wmT),nrow(wmT))
        
        
        
        for(scen.index in 1:length(prior_scen_list)){
          tox<-scen.index
          
          # browser()
          res_list<- get(paste(c("JOINTTITEPOCRM.prior2_sk",3,"_alpha_mu",mu_alphaT_i,"_alpha_tau",tau_alphaT_i,"_beta_mu",mu_betaT_i,"_beta_tau",tau_betaT_i,"_scen",tox,"_start4"),collapse=""))
          
          nsims<-length(res_list)
          dose_rec_vector<-unlist(lapply(res_list,'[[','dose.rec'))
          MTD<-correct_list[scen.index]
          if(MTD=="UNSAFE"){
            stop_list<-(lapply(res_list,'[[','stop.code'))
            stop_mat<-do.call(rbind,stop_list)
            #total of safety
            safety<-(rowSums(stop_mat[,c(1,5,6)],na.rm=T)>0)
            result_array[scen.index,mu_alphaT_i,mu_betaT_i,tau_alphaT_i,tau_betaT_i]<-100*sum(safety)/nsims
          }else if(MTD=="TOOSAFE"){
            stop_list<-(lapply(res_list,'[[','stop.code'))
            stop_mat<-do.call(rbind,stop_list)
            #total of safety
            safety1<-(stop_mat[,7]>0)
            result_array[scen.index,mu_alphaT_i,mu_betaT_i,tau_alphaT_i,tau_betaT_i]<-100*sum(safety1,na.rm=T)/nsims
            
          }else{
            result_array[scen.index,mu_alphaT_i,mu_betaT_i,tau_alphaT_i,tau_betaT_i]<-100*sum(dose_rec_vector==MTD,na.rm=T)/nsims
          }               
          
          
          
        }#for scen
      } #for mu beta
    } #for mu alpha
  } #for tau beta
} #for tau alpha


#geometric mean
geomean_array<-exp(colMeans(log(result_array)))
#choice of parameters (index)
geomean_choice<-which(geomean_array==max(geomean_array),arr.ind = T)
#choice of parameters (parameters)
par_choices<-c(mu_alphaT_priors[geomean_choice[1]],# mean alpha
               mu_betaT_priors[geomean_choice[2]],# mean beta
               tau_alphaT_priors[geomean_choice[3]],# prec alpha
               tau_betaT_priors[geomean_choice[4]])# prec beta) 


plot(c(1:length(geomean_array)),geomean_array,ylim=c(0,35))

