##prior results extractions 


correct_func<-function(prior_scens,good_tox=0.3,good_eff=0.2,w1=0.33,w2=1.09){
  correct_vec<-c()
  for(sc.i in 1:length(prior_scens)){
    sc<- prior_scens[sc.i]
    eff_probs<-get(paste(c("eff_scen",sc),collapse=""))
    tox_probs<-get(paste(c("tox_scen",sc),collapse=""))
    utilities<-eff_probs-w1*tox_probs-w2*tox_probs*as.numeric(tox_probs>0.3)
    if(all(tox_probs>good_tox)){
      correct_vec[sc.i]<-NA
    }else if(all(tox_probs<good_tox)){
      correct_vec[sc.i]<-11
    }else if(all(eff_probs<good_eff)){
      
      correct_vec[sc.i]<-NA
    }else{
      corrects<-(tox_probs<=good_tox)&(eff_probs>=good_eff)
      if(sum(corrects)>0){
        #   browser()
        correct_vec[sc.i]<-which(corrects)[which.max(utilities[corrects])]
      }else{
        correct_vec[sc.i]<-NA
      }
      
    }
    
    
  }
  
  return(correct_vec)
}


correct_vector<-correct_func(scens)

#POCRM
load("JOINT_TITEPOCRM_TCAL1.RData")
#empty results array
result_array<-array(NA,dim=c(length(scens),
                             
                             3,#pars
                             4 #skeleton
))      



for(prior.par in c(1:3)){
  for(prior.sk in c(1:4)){ 
    
    for(sc.i in 1:length(scens)){
      scen<- scens[sc.i]
      
      
      res_list<-get(paste(c("JOINTTITEPOCRM.par",prior.par,".sk",prior.sk,"_",scen),collapse=""))
      
      nsims<-length(res_list)
      dose_rec_vector<-unlist(lapply(res_list,'[[','dose.rec'))
      OBD<-correct_vector[sc.i]
      if(is.na(OBD)){
        stop_list<-(lapply(res_list,'[[','stop.code'))
        stop_mat<-do.call(rbind,stop_list)
        #total of stops
        stops<-(rowSums(stop_mat[,c(1,5,6)],na.rm=T)>0)
        result_array[sc.i,prior.par,prior.sk]<-100*sum(stops)/nsims
      }else{
        result_array[sc.i,prior.par,prior.sk]<-100*sum(dose_rec_vector==OBD,na.rm=T)/nsims
        
      }
    }
  }
}



#geometric mean
geomean_array<-exp(colMeans(log(result_array)))
#choice of parameters (index)
geomean_choice<-which(geomean_array==max(geomean_array),arr.ind = T)
#plot geomean
plot(c(1:length(geomean_array)),geomean_array)


#BLRM
load("JOINT_TITEBLRM_TCAL1.RData")


#empty results array
result_array<-array(NA,dim=c(length(scens),
                             
                             3,#Z pars
                             3 #W pars
))      



for(prior.Z.par in c(1:3)){
  for(prior.W.par in c(1:3)){ 
    
    for(sc.i in 1:length(scens)){
      scen<- scens[sc.i]
      
      res_list<-get( paste(c("JOINTTITEBLRM.parT.W",prior.W.par,".Z",prior.Z.par,"_",scen),collapse=""))
      
      nsims<-length(res_list)
      dose_rec_vector<-unlist(lapply(res_list,'[[','dose.rec'))
      OBD<-correct_vector[sc.i]
      if(is.na(OBD)){
        stop_list<-(lapply(res_list,'[[','stop.code'))
        stop_mat<-do.call(rbind,stop_list)
        #total of stops
        stops<-(rowSums(stop_mat[,c(1,5,6)],na.rm=T)>0)
        result_array[sc.i,prior.Z.par,prior.W.par]<-100*sum(stops)/nsims
      }else{
        result_array[sc.i,prior.Z.par,prior.W.par]<-100*sum(dose_rec_vector==OBD,na.rm=T)/nsims
        
      }
    }
  }
}




#geometric mean
geomean_array<-exp(colMeans(log(result_array)))
#choice of parameters (index)
geomean_choice<-which(geomean_array==max(geomean_array),arr.ind = T)
#plot geomean
plot(c(1:length(geomean_array)),geomean_array)
