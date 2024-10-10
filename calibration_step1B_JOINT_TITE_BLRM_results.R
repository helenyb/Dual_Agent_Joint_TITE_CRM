##blrm prior calibration results extractions stage 2

#start at 4
#all scens
load("JOINT_TITE_BLRM_priorcal2.RData")

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
                             length(mu_alphaT_W_priors),
                             length(mu_betaT_W_priors),
                             length(mu_alphaT_Z_priors),
                             length(mu_betaT_Z_priors),
                             length(tau_alphaT_priors),
                             length(tau_betaT_priors)
))



for (mu_alphaT_W_i in 1:length(mu_alphaT_W_priors)){
  for (mu_betaT_W_i in 1:length(mu_betaT_W_priors)){
    for (mu_alphaT_Z_i in 1:length(mu_alphaT_Z_priors)){
      for (mu_betaT_Z_i in 1:length(mu_betaT_Z_priors)){
        for (tau_alphaT_i in 1:length(tau_alphaT_priors)){
          for (tau_betaT_i in 1:length(tau_betaT_priors)){
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
              
              
              res_list<- get(paste(c("JOINTTITEBLRM.prior_","alpha_tau",tau_alphaT_i,"_beta_tau",tau_betaT_i,
                                     "alpha_mu_Z",mu_alphaT_Z_i,"_beta_mu_Z",mu_betaT_Z_i,"_alpha_mu_W",mu_alphaT_W_i,"_beta_mu_W",mu_betaT_W_i,"_scen",tox,"_start4"),collapse=""))
              
              nsims<-length(res_list)
              dose_rec_vector<-unlist(lapply(res_list,'[[','dose.rec'))
              MTD<-correct_list[scen.index]
              if(MTD=="UNSAFE"){
                stop_list<-(lapply(res_list,'[[','stop.code'))
                stop_mat<-do.call(rbind,stop_list)
                #total of safety
                safety<-(rowSums(stop_mat[,c(1,5,6)],na.rm=T)>0)
                result_array[scen.index,mu_alphaT_W_i,mu_betaT_W_i,mu_alphaT_Z_i,mu_betaT_Z_i,tau_alphaT_i,tau_betaT_i]<-100*sum(safety)/nsims
              }else if(MTD=="TOOSAFE"){
                stop_list<-(lapply(res_list,'[[','stop.code'))
                stop_mat<-do.call(rbind,stop_list)
                #total of safety
                safety1<-(stop_mat[,7]>0)
                result_array[scen.index,mu_alphaT_W_i,mu_betaT_W_i,mu_alphaT_Z_i,mu_betaT_Z_i,tau_alphaT_i,tau_betaT_i]<-100*sum(safety1,na.rm=T)/nsims
                
              }else{
                result_array[scen.index,mu_alphaT_W_i,mu_betaT_W_i,mu_alphaT_Z_i,mu_betaT_Z_i,tau_alphaT_i,tau_betaT_i]<-100*sum(dose_rec_vector==MTD,na.rm=T)/nsims
              }               
              
              
              
            }#for scen
            
          } #for tau beta
        } #for tau alpha
      } #for mu betaZ
    } #for mu alphaZ
  } #for mu betaW
} #for mu alphaW




#geometric mean
geomean_array<-exp(colMeans(log(result_array)))
#choice of parameters (index)
geomean_choice<-which(geomean_array==max(geomean_array),arr.ind = T)
#choice of parameters (parameters)
par_choices<-c(mu_alphaT_W_priors[geomean_choice[1]],
               mu_betaT_W_priors[geomean_choice[2]],
               mu_alphaT_Z_priors[geomean_choice[3]],
               mu_betaT_Z_priors[geomean_choice[4]],
               tau_alphaT_priors[geomean_choice[5]],
               tau_betaT_priors[geomean_choice[6]] )

#plotting to check spread
plot(c(1:length(geomean_array)),geomean_array,ylim=c(0,35))


