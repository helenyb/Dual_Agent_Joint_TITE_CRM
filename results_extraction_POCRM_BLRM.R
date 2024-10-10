##For utility based on CRM:

#load:
load("JOINT_TITE_BLRM_sims.RData")
load("JOINT_TITE_POCRM_sims.RData")

methods_vec<-c("JOINTTITEPOCRM","JOINTTITEBLRM")

dose.names.vec<-c()
for(i in 1:nrow(dose.ind.mat)){
  dose.names.vec[i]<-paste(c("[",dose.ind.mat[i,1],",",dose.ind.mat[i,2],"]"),collapse="")
}



dose_rec_frame<-function(method,eff_sc,tox_sc,eff_pattern){
  res_list<-get(paste(c(method,".eff",eff_sc,".",eff.pattern,"_tox",tox_sc ),collapse=""))
  
  dose_rec_vector<-unlist(lapply(res_list,'[[','dose.rec'))
  for(i in 1:length(dose_rec_vector)){
    if(res_list[[i]]$stop.code[7]==1){
      dose_rec_vector[i]<-11
    }
    if(res_list[[i]]$stop.code[1]==1){
      dose_rec_vector[i]<-NA
    }
  }
  
  
  out<-data.frame(c("N",dose.names.vec,"H"),c(sum(is.na(dose_rec_vector)),tabulate(dose_rec_vector,nbins=ndoses+1))/(length(dose_rec_vector)/100),method)
  names(out)<-c("Dose","Rec","Method")    
  out$Dose<-factor(out$Dose,levels=unique(c("N",dose.names.vec,"H")))
  
  return(out)
}



#dose.rec matrix last column is (1) all too safe (2) no admissables
dose_rec_matrix<-function(method,eff_sc,tox_sc,eff_pattern){
  res_list<-get(paste(c(method,".eff",eff_sc,".",eff.pattern,"_tox",tox_sc ),collapse=""))
  
  dose_rec_vector<-unlist(lapply(res_list,'[[','dose.rec'))
  for(i in 1:length(dose_rec_vector)){
    if(res_list[[i]]$stop.code[7]==1){
      dose_rec_vector[i]<-11
    }
    if(res_list[[i]]$stop.code[1]==1){
      dose_rec_vector[i]<-NA
    }
  }
  
  
  out<-data.frame(c("N",dose.names.vec,"H"),c(sum(is.na(dose_rec_vector)),tabulate(dose_rec_vector,nbins=ndoses+1))/(length(dose_rec_vector)/100),method)
  out.matrix<-matrix(c(tabulate(dose_rec_vector,nbins=ndoses+1),sum(is.na(dose_rec_vector)))/(length(dose_rec_vector)/100),nrow=2,ncol=6)
  return(out.matrix)
}


#dose assignment matrix 
dose_assi_matrix<-function(method,eff_sc,tox_sc,eff_pattern){
  res_list<-get(paste(c(method,".eff",eff_sc,".",eff.pattern,"_tox",tox_sc ),collapse=""))
  
  dose_assi_matall<-matrix(3*unlist(lapply(res_list,'[[','dose.ass')),nrow=10)
  out_mat<-matrix(rowMeans(dose_assi_matall),nrow=2)
  return(out_mat)
}

for(method in methods_vec){
  for(eff.pattern in 1){
    for(eff in 1:6){
      for(tox in 1:6){
        assign(paste(c("rec_frame_",method,"_eff",eff,".",eff.pattern,"_tox",tox),collapse=""),
               dose_rec_frame(method=method,eff_sc=eff,tox_sc=tox,eff_pattern=eff.pattern))
        
      }
    }
  }
}

#pcs
correct_mat_func<-function(eff_sc_vec,tox_sc_vec,good_tox=0.3,good_eff=0.2,w1=0.33,w2=1.09){
  correct_mat<-matrix(nrow=length(eff_sc_vec),ncol=length(tox_sc_vec))
  for(eff_sc in eff_sc_vec){
    for(tox_sc in tox_sc_vec){
      eff_probs<-get(paste(c("eff_scen",eff_sc),collapse=""))
      tox_probs<-get(paste(c("tox_scen",tox_sc),collapse=""))
      utilities<-eff_probs-w1*tox_probs-w2*tox_probs*as.numeric(tox_probs>0.3)
      if(all(tox_probs>good_tox)){
        correct_mat[eff_sc,tox_sc]<-NA
      }else if(all(tox_probs<good_tox)){
        correct_mat[eff_sc,tox_sc]<-11
      }else if(all(eff_probs<good_eff)){
        
        correct_mat[eff_sc,tox_sc]<-NA
      }else{
        corrects<-(tox_probs<=good_tox)&(eff_probs>=good_eff)
        if(sum(corrects)>0){
          #   browser()
          correct_mat[eff_sc,tox_sc]<-which(corrects)[which.max(utilities[corrects])]
        }else{
          correct_mat[eff_sc,tox_sc]<-NA
        }
        
      }
      
    }
  }
  
  return(correct_mat)
}

#pas
acceptable_list_func<-function(eff_sc_vec,tox_sc_vec,upper_tox=0.33,ok_eff=0.2,good_tox=0.3){
  acceptable_list<-list()
  
  for(eff_sc in eff_sc_vec){
    acceptable_list[[eff_sc]]<-list()
    for(tox_sc in tox_sc_vec){
      eff_probs<-get(paste(c("eff_scen",eff_sc),collapse=""))
      tox_probs<-get(paste(c("tox_scen",tox_sc),collapse=""))
      if(all(tox_probs>upper_tox)){
        acceptable_list[[eff_sc]][[tox_sc]]<-NA
      }else if(all(eff_probs<ok_eff)){
        acceptable_list[[eff_sc]][[tox_sc]]<-NA
      }else{
        acceptables<-(tox_probs<=upper_tox)&(eff_probs>=ok_eff)
        
        if(sum(acceptables)>0){
          #   browser()
          which.acceptables<-which(acceptables)
          if(all(tox_probs<good_tox)){
            which.acceptables<-c(which.acceptables,11)
          }
          acceptable_list[[eff_sc]][[tox_sc]]<-which.acceptables
        }else{
          acceptable_list[[eff_sc]][[tox_sc]]<-NA
        }
        
      }
      
    }
    
  }
  return(acceptable_list)
}


#good
good_list_func<-function(eff_sc_vec,tox_sc_vec,upper_tox=0.33,ok_eff=0.2,good_tox=0.3,
                         good_eff=0.2,w1=0.33,w2=1.09,utility.margin=0.15){
  
  
  correct_mat<-matrix(nrow=length(eff_sc_vec),ncol=length(tox_sc_vec))
  acceptable_list<-list()
  good_list<-list()
  for(eff_sc in eff_sc_vec){
    acceptable_list[[eff_sc]]<-list()
    good_list[[eff_sc]]<-list()
    for(tox_sc in tox_sc_vec){
      eff_probs<-get(paste(c("eff_scen",eff_sc),collapse=""))
      tox_probs<-get(paste(c("tox_scen",tox_sc),collapse=""))
      utilities<-eff_probs-w1*tox_probs-w2*tox_probs*as.numeric(tox_probs>0.3)
      
      if(all(tox_probs>upper_tox)){
        acceptable_list[[eff_sc]][[tox_sc]]<-NA
        good_list[[eff_sc]][[tox_sc]]<-NA
      }else if(all(eff_probs<ok_eff)){
        acceptable_list[[eff_sc]][[tox_sc]]<-NA
        good_list[[eff_sc]][[tox_sc]]<-NA
      }else{
        acceptables<-(tox_probs<=upper_tox)&(eff_probs>=ok_eff)
        
        if(sum(acceptables)>0){
          #   browser()
          which.acceptables<-which(acceptables)
          if(all(tox_probs<good_tox)){
            which.acceptables<-c(which.acceptables,11)
          }
          acceptable_list[[eff_sc]][[tox_sc]]<-which.acceptables
          
          max.utility<-max(utilities[which.acceptables])
          utility.margin.boundary<-max.utility-utility.margin
          acceptable.utilities<-rep(-100,length(eff_probs))
          acceptable.utilities[which.acceptables]<-utilities[which.acceptables]
          which.good<-which(acceptable.utilities>=utility.margin.boundary)
          good_list[[eff_sc]][[tox_sc]]<-which.good
          
        }else{
          acceptable_list[[eff_sc]][[tox_sc]]<-NA
          good_list[[eff_sc]][[tox_sc]]<-NA
        }
        
      }
      
      
      if(all(tox_probs>good_tox)){
        correct_mat[eff_sc,tox_sc]<-NA
      }else if(all(tox_probs<good_tox)){
        correct_mat[eff_sc,tox_sc]<-11
      }else if(all(eff_probs<good_eff)){
        
        correct_mat[eff_sc,tox_sc]<-NA
      }else{
        corrects<-(tox_probs<=good_tox)&(eff_probs>=good_eff)
        if(sum(corrects)>0){
          #   browser()
          correct_mat[eff_sc,tox_sc]<-which(corrects)[which.max(utilities[corrects])]
        }else{
          correct_mat[eff_sc,tox_sc]<-NA
        }
        
      }
      
    }
  }
  
  
  
  
  return(good_list)
}




correct_acceptable_good_func<-function(method,eff_sc_vec,tox_sc_vec,eff_pattern_vec,w1=0.33,w2=1.09,utility.margin=0.1){
  correct_perc_array<-array(dim = c(length(eff_sc_vec),length(tox_sc_vec),length(eff_pattern_vec)))
  acceptable_perc_array<-array(dim = c(length(eff_sc_vec),length(tox_sc_vec),length(eff_pattern_vec)))
  good_perc_array<-array(dim = c(length(eff_sc_vec),length(tox_sc_vec),length(eff_pattern_vec)))
  
  acceptable_list<-acceptable_list_func(eff_sc_vec = eff_sc_vec,tox_sc_vec = tox_sc_vec)
  correct_mat<-correct_mat_func(eff_sc_vec = eff_sc_vec,tox_sc_vec = tox_sc_vec,w1=w1,w2=w2)
  good_list<-good_list_func(eff_sc_vec = eff_sc_vec,tox_sc_vec = tox_sc_vec,w1=w1,w2=w2,utility.margin=utility.margin)
  
  for(eff.pattern in eff_pattern_vec){
    for(eff_sc in eff_sc_vec){
      for(tox_sc in tox_sc_vec){   
        res_list<-get(paste(c(method,".eff",eff_sc,".",eff.pattern,"_tox",tox_sc ),collapse=""))
        #        res_list<-get(paste(c(method,".eff",eff_sc,".",eff.pattern,"_tox",tox_sc),collapse=""))
        dose_rec_vector<-unlist(lapply(res_list,'[[','dose.rec'))
        for(i in 1:length(dose_rec_vector)){
          if(res_list[[i]]$stop.code[7]==1){
            dose_rec_vector[i]<-11
          }
          if(res_list[[i]]$stop.code[1]==1){
            dose_rec_vector[i]<-NA
          }
        }
        
        
        if(is.na(correct_mat[eff_sc,tox_sc])){
          correct_perc_array[eff_sc,tox_sc,eff.pattern]<-100*mean(is.na(dose_rec_vector))
        }else{
          # browser()
          correct_perc_array[eff_sc,tox_sc,eff.pattern]<-100*sum(dose_rec_vector==correct_mat[eff_sc,tox_sc],na.rm=T)/length(dose_rec_vector)
          
        }
        
        if(any(is.na(acceptable_list[[eff_sc]][[tox_sc]]))){
          acceptable_perc_array[eff_sc,tox_sc,eff.pattern]<-100*mean(is.na(dose_rec_vector))
          good_perc_array[eff_sc,tox_sc,eff.pattern]<-100*mean(is.na(dose_rec_vector))
          
        }else{
          acceptable_perc_array[eff_sc,tox_sc,eff.pattern]<-100*mean(dose_rec_vector %in% (acceptable_list[[eff_sc]][[tox_sc]]))
          good_perc_array[eff_sc,tox_sc,eff.pattern]<-100*mean(dose_rec_vector %in% (good_list[[eff_sc]][[tox_sc]]))
          
        }
      }
    }
  }
  #transposed to match excel sheet (one efficacy pattern only)
  return(list(correct_perc_array=t(correct_perc_array[,,1]),
              acceptable_perc_array=t(acceptable_perc_array[,,1]),
              good_perc_array=t(good_perc_array[,,1])))
}   




for(method in methods_vec){
  assign(paste(c("correct_acc_mt_",method),collapse=""),
         correct_acceptable_good_func(method=method,eff_sc_vec=c(1:6),tox_sc_vec=c(1:6),eff_pattern_vec=1,w1=0.33,w2=1.09,utility.margin = 0.1)
  )
  
  
}


######Sample Size
SampleSize_func<-function(method,eff_sc,tox_sc,eff_pattern){
  res_list<-get(paste(c(method,".eff",eff_sc,".",eff.pattern,"_tox",tox_sc ),collapse=""))
  
  num_pat_vector<-unlist(lapply(res_list,'[[','num.pat'))
  
  return(mean(num_pat_vector))
  
}

SampleSize_func_all<-function(method,eff_sc_vec,tox_sc_vec,eff_pattern){
  out.mat<-matrix(ncol=max(eff_sc_vec),nrow=max(tox_sc_vec))
  
  for(eff in eff_sc_vec){
    for(tox in tox_sc_vec){
      
      
      res_list<-get(paste(c(method,".eff",eff,".",eff.pattern,"_tox",tox ),collapse=""))
      
      out.mat[tox,eff]<-mean(unlist(lapply(res_list,'[[','num.pat')))
    }
  }
  
  return(out.mat)
  
}



for(method in methods_vec){
  
  assign(paste(c("SampleSize_mt_",method),collapse=""),
         SampleSize_func_all(method=method,eff_sc_vec=c(1:6),tox_sc_vec=c(1:6),eff_pattern=1)
  )
  
  
}


##unsafe dose assignment
unsafe_list_func<-function(eff_sc_vec,tox_sc_vec,upper_tox=0.3){
  unsafe_list<-list()
  
  for(eff_sc in eff_sc_vec){
    unsafe_list[[eff_sc]]<-list()
    for(tox_sc in tox_sc_vec){
      eff_probs<-get(paste(c("eff_scen",eff_sc),collapse=""))
      tox_probs<-get(paste(c("tox_scen",tox_sc),collapse=""))
      
      unsafe<-(tox_probs>upper_tox)
      
      if(sum(unsafe)>0){
        #   browser()
        which.unsafe<-which(unsafe)
      }else{
        which.unsafe<-0
      }
      unsafe_list[[eff_sc]][[tox_sc]]<-which.unsafe
      
      
    }
    
  }
  
  
  return(unsafe_list)
}



sum_unsafe<-function(dose_vector,unsafe_doses,co_size,ndoses){
  all.assigned<-rep(dose_vector,each=co_size)
  sum.assigned<-c()
  for(dose.level in 1:ndoses){
    sum.assigned[dose.level]<-sum(all.assigned==dose.level,na.rm=T)
  }
  if(sum(unsafe_doses)==0){
    tot.sum<-0
  }else{
    tot.sum<-sum(sum.assigned[unsafe_doses])
  }
  return(tot.sum)
}

unsafe_func<-function(method,eff_sc_vec,tox_sc_vec,eff_pattern_vec=1,co_size=3,ndoses){
  unsafe_perc_array<-matrix(ncol=length(eff_sc_vec),nrow=length(tox_sc_vec))
  
  unsafe_list<-unsafe_list_func(eff_sc_vec = eff_sc_vec,tox_sc_vec = tox_sc_vec)
  for(eff.pattern in eff_pattern_vec){
    for(eff in eff_sc_vec){
      for(tox in tox_sc_vec){   
        res_list<-get(paste(c(method,".eff",eff,".",eff.pattern,"_tox",tox ),collapse=""))
        
        dose_ass_list<-(lapply(res_list,'[[','dosevec'))
        
        unsafe_perc_array[tox,eff]<-mean(unlist(lapply(dose_ass_list,sum_unsafe,unsafe_doses=unsafe_list[[eff]][[tox]],co_size=co_size,ndoses=ndoses)))
        #browser()
      }
    }
  }
  return(unsafe_perc_array=unsafe_perc_array)
}   

for(method in methods_vec){
  
  assign(paste(c("Unsafe_assigment_mt_",method),collapse=""),
         unsafe_func(method=method,eff_sc_vec=c(1:6),tox_sc_vec=c(1:6),eff_pattern_vec=1,co_size=3,ndoses=10)
  )
  
  
}

##########
