#TITE-BOIN12 COMBO

source("data_generation_TTE_v2.R")

##pre-amble functions


#logit & expit
logit<-function(p){
  return(log(p/(1-p)))
}
expit<-function (x) {
  return(exp(x)/(1 + exp(x)))
}

##this function breaks ties by choosing lowest as default (for comparability in sims), specify anything else in break.ties arg for random
which.is.max.func<-function(x,break.ties="min"){
  if(break.ties=="min"){
    if(is.vector(x)){
      y<-seq_along(x)[x==max(x)]
      if (length(y) > 1L) 
        return(y[1])
      else return(y)
    }else if(is.matrix(x)){
      y<-which(x==max(x),arr.ind = T)
      if(nrow(y)>1L){
        return(y[1,])
      }else{
        return(y)
      }
    }
  }else{
    if(is.vector(x)){
      y<-seq_along(x)[x==max(x)]
      if (length(y) > 1L) 
        return(sample(y, 1L))
      else return(y)
    }else if(is.matrix(x)){
      y<-which(x==max(x),arr.ind = T)
      if(nrow(y)>1L){
        return(y[sample(c(1:nrow(y)),1L),])
      }else{
        return(y)
      }
    }
  }
}





current_patient_data_frame<-function(current_time,patient.dataframe,follow_up){
  #input the whole data and convert to format for TITECRM
  #each row is a patient, with entries c("entry.time","dose.level","DLT","DLT.time","Eff","Eff.time") 
  num_patients<-max(patient.dataframe$patient_ID)
  patient.dataframe<-patient.dataframe[patient.dataframe$time_of<=current_time,]
  patient.dataframe<-patient.dataframe[patient.dataframe$cycle_num<=follow_up,]
  
  
  patient_dataframe<-matrix(NA,nrow=num_patients,ncol=6)
  for (i in 1:num_patients){
    patient_data_ind<-patient.dataframe[patient.dataframe$patient_ID==i,]
    patient_DLT<-max(patient_data_ind$DLT)
    patient_Eff<-max(patient_data_ind$Eff)
    patient_entry.time<-patient_data_ind$entry_time[1]
    patient_dose.level<-patient_data_ind$dose_level[1]
    if(patient_DLT==1){
      patient_DLT.time<-max(patient_data_ind$DLT.time,na.rm=T)
    }else{
      patient_DLT.time<-NA
    }
    
    if(patient_Eff==1){
      patient_Eff.time<-max(patient_data_ind$Eff.time,na.rm=T)
    }else{
      patient_Eff.time<-NA
    }
    patient_dataframe[i,]<-c(patient_entry.time,patient_dose.level,patient_DLT,patient_DLT.time,patient_Eff,patient_Eff.time)
  }
  
  patient_dataframe<-data.frame(patient_dataframe)
  names(patient_dataframe)<-c("entry.time","dose.level","DLT","DLT.time","Eff","Eff.time") 
  follow_up_time<-function(x) min(follow_up,x)
  patient_follow_up<- sapply(current_time-patient_dataframe$entry.time,follow_up_time)
  patient_weights_tox<-patient_follow_up/follow_up
  patient_weights_eff<-patient_follow_up/follow_up
  
  #only DLT if we have seen it at the current time
  current.DLT<-patient_dataframe$DLT
  current.DLT[which(current.DLT==1)]<-patient_dataframe$DLT.time[which(current.DLT==1)]<=current_time
  
  #only DLT time observed if current DLT is true
  current.DLT.time<-patient_dataframe$DLT.time
  current.DLT.time[which(current.DLT==0)]<-NA
  
  #only Eff if we have seen it at the current time
  current.Eff<-patient_dataframe$Eff
  current.Eff[which(current.Eff==1)]<-patient_dataframe$Eff.time[which(current.Eff==1)]<=current_time
  
  #only Eff time observed if current Eff is true
  current.Eff.time<-patient_dataframe$Eff.time
  current.Eff.time[which(current.Eff==0)]<-NA
  
  entry.time<-patient_dataframe$entry.time
  dose.level<-patient_dataframe$dose.level
  
  patient_weights_tox[which(current.DLT==1)]<-1
  patient_weights_eff[which(current.DLT==1)]<-(current.DLT.time[which(current.DLT==1)]-entry.time[which(current.DLT==1)])/follow_up
  patient_weights_eff[which(current.Eff==1)]<-1
  
  
  
  return(data.frame(entry.time,dose.level,current.DLT,current.DLT.time,current.Eff,current.Eff.time,patient_weights_tox,patient_weights_eff,patient_follow_up))
  
}


#create hard safety matrix

hard.safety.mat.function<-function(perc,co.size,max.cohorts){
  hard.safety.mat<-matrix(nrow=2,ncol=max.cohorts)
  hard.safety.mat[2,]<-seq(from=co.size,to=max.cohorts*co.size, by=co.size)
  for(i in 1:max.cohorts){
    numi<- hard.safety.mat[2,i]
    try_vec<-c(1:numi)
    probs_vec<-1-pbeta(0.3,1+try_vec,1+numi-try_vec)
    hard.safety.mat[1,i]<- min( which(100*probs_vec>perc))
  }                     
  return(hard.safety.mat)
}





##INPUT:
#seed=seed for reproducability
#tru.T.pars=parameters for data generation of DLT times 
#tru.E.pars=parameters for data generation of efficacy times 
#tru.corET=correlation between DLT times and efficacy times in data generation
#co_size=cohort size
#ncohorts=maximum number of cohorts in study
#target=target toxicity (all cycles)
#ncycles=number of cycles
#dose.skipping.rule= type of dose.skipping rule:
#"ON.diag.allowed" is same as "ON" but also allows (a+1,b+1)
#sufficient.information==enforce stopping for sufficient information? no more than 9 patients per dose. T=enforce stopping for sufficient information (default)
#sufficient.information.lim= number of patients needed before stopping when the next assignment is the same. 
#hard.safety.rule=percentage for hard safety rule based on Beta(1.1)? 
#e.g. : 85= 2/3,3/6,4/9. 90=2/3,4/6,5/9, 95=3/3,4/6,5/9, <50 means no hard safety enforced
#safety.stopping.low.unsafe= (T=stop when P(p1>0.3)>0.8 (cycle 1))
#safety.stopping.high.toosafe= (T= stop when P(pJ>0.3)>0.8 (cycle 1))
#initial.one.cycle: Is the initial period based on one cycle at a time? (F=wait until all cycles completed before next dose in initial period)
#C_eff: Dose is "admissible in efficacy if P(P(efficacy)>effbound)>C_eff
#C_tox: Dose is "admissible in safety if P(P(DLT)<toxbound)>C_tox
#effbound: Dose is "admissible in efficacy if P(P(efficacy)>effbound)>C_eff
#toxbound: Dose is "admissible in safety if P(P(DLT)<toxbound)>C_tox
#pause: Recruitment paused if no admissible doses, to allow observation of current patients more cycles
#a.stop.bound: number of patients required before activity is allowed to stop trial



#backfill: should doses deemed safe be backfilled? (Default FALSE)
#backfill.num: How many cohorts to add as backfilling?
#dose.indices: matrix with 2 col and ndoses rows, of indices for doses (i.e row 2 tells us dose "2" is [2,1] on the grid)
#default.order: Order in which we escalate if we see no activity/toxicity
#start.dose: The starting dose, as defined in the default order



##OUTPUT:
# dose.rec: dose recommendation
# num.pat: number of patients
# dose.ass: number of cohorts per dose (vector)
# stop.code: stopping reason 
#1= No admissible doses
#2= Precision (N/A for pocrm)
#3= Max patients
#4= Sufficient information
#5= Lowest dose fails hard safety
#6= Lowest dose unsafe
#7= Highest Dose too safe
# num.DLT: number of DLTS (total)
# DLT.mat: 2xnumdoses matrix. rows for grades 1=no DLT, 2=DLT. cols for doses.
# Eff.mat: 2xnumdoses matrix. rows for grades 1=noEff, 2=Eff. cols for doses.
# duration: total trial duration (until all recruited patients are fully observed)
# max_admissable max admissable dose at the end of the trial according to hard safety only (has hard safety eliminated any?)
# dosevec: Dose assignment of cohorts in sequence order
# DLT.vec: Binary sequence of DLT outcomes (all cycles)
# EFF.vec: Binary sequence of efficacy outcomes (all cycles)
# MTD.rec: MTD recommendation (note this is not subject to safety rules, will always output a dose)


JOINT.TITE.BOIN<-function(seed,tru.E.pars,tru.T.pars,tru.corET,co_size,ncohorts ,target,targetE,
                             ncycles,dose.skipping.rule="ON", 
                             sufficient.information=T,sufficient.information.lim=29,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                             safety.stopping.high.toosafe=T,initial.one.cycle=T,
                             C_eff=0.2,C_tox=0.2,effbound=0.2,toxbound=0.3,
                             backfill=F,backfill.num=2,TITE=T,dose.indices,
                             default.order,start.dose,utility_mat,  lower_tox_lim=0.25,
                             upper_tox_lim=0.35,pause=T,a.stop.bound){
  set.seed(seed)
  stop_vec<-rep(0,7)
  ndoses<-nrow(dose.indices)
  dosevec<-c()
  patient_ID1<-1
  current.time<-0
  
  excluded<-rep(0,ndoses)
  
  
  #check that the data generation inputs tally with the number of doses and cycles
  if(!all(dim(tru.T.pars)==dim(tru.E.pars))){
    stop("Incompatible truth matrices for toxicity and efficacy")
  }
  
  
  if(ndoses!=ncol(tru.T.pars)){
    stop("Incompatable truth matrix with number of doses")
  }
  
  if((TITE==F)&(ncycles!=1)){
    stop("Incompatable TITE and ncycles input: ncycles must equal 1 for TITE=F")
  }
  
  
  #safety for BOIN
  
  
  delta_L<-target-(log((1-lower_tox_lim)/(1-target)))/(log((target*(1-lower_tox_lim))/(lower_tox_lim*target)))
  delta_U<-target+(log((1-target)/(1-upper_tox_lim)))/(log((upper_tox_lim*(1-target))/(target*upper_tox_lim)))
  
  
  WT<-0.11
  # u00<-utility_mat[1,1]
  # u01<-utility_mat[1,2]
  # u10<-utility_mat[2,1]
  # u11<-utility_mat[2,2]
  # u_bar<-u00*(1-targetE)*(1-target)+u01*(targetE)*(1-target)+u11*(targetE)*(target)
  #benchmark
  u_bar<-targetE-WT*target
  #u_b<-(u_bar+(100-u_bar)/2)/100
  u_b<-((0.5*u_bar +0.5)+WT)/(1+WT)
  
  
  ##definitions of allowable escalations for dose.skipping.rule
  #define the matrix
  dose.mat<-matrix(nrow=max(dose.indices[,1]),ncol=max(dose.indices[,2]))
  for(dose_i in 1:nrow(dose.indices)){
    dose.mat[dose.indices[dose_i,1],dose.indices[dose_i,2]]<-dose_i
  }
  
  #define the allowable.escalations list
  if(dose.skipping.rule=="ON.diag.allowed"){
    allowable.under.int<-list()
    allowable.in.int<-list()
    allowable.over.int<-list()
    for(dose_i in 1:nrow(dose.indices)){
      
      dose_index<-dose.indices[dose_i,]
      all.below<-which((dose.indices[,1]<=dose_index[1])&(dose.indices[,2]<=dose_index[2]))
      
      #de-escalation
      x_minus_1_y<-which((dose.indices[,1]==(dose_index[1]-1))&(dose.indices[,2]==dose_index[2]))
      x_y_minus_1<-which((dose.indices[,1]==dose_index[1])&(dose.indices[,2]==(dose_index[2]-1)))
      
      #off-diagonal
      x_plus_1_y_minus_1<-which((dose.indices[,1]==(dose_index[1]+1))&(dose.indices[,2]==(dose_index[2]-1)))
      x_minus_1_y_pluus_1<-which((dose.indices[,1]==(dose_index[1]-1))&(dose.indices[,2]==(dose_index[2]+1)))
      
      #escalation
      x_plus_1_y<-which((dose.indices[,1]==(dose_index[1]+1))&(dose.indices[,2]==dose_index[2]))
      x_y_plus_1<-which((dose.indices[,1]==dose_index[1])&(dose.indices[,2]==(dose_index[2]+1)))
      x_plus_1_y_plus_1<-which((dose.indices[,1]==(dose_index[1]+1))&(dose.indices[,2]==(dose_index[2]+1)))
      
      
      
      allowable.under.int[[dose_i]]<-unique(c(x_minus_1_y,x_y_minus_1,
                                              x_plus_1_y_minus_1,x_minus_1_y_pluus_1,
                                              x_plus_1_y,x_y_plus_1,x_plus_1_y_plus_1,
                                              dose_i))
      allowable.in.int[[dose_i]]<-unique(c(x_minus_1_y,x_y_minus_1,
                                           dose_i))
      allowable.over.int[[dose_i]]<-unique(c(x_minus_1_y,x_y_minus_1))
    }
  }else{
    allowable.under.int<-list()
    allowable.in.int<-list()
    allowable.over.int<-list()
    for(dose_i in 1:nrow(dose.indices)){
      
      dose_index<-dose.indices[dose_i,]
      all.below<-which((dose.indices[,1]<=dose_index[1])&(dose.indices[,2]<=dose_index[2]))
      
      #de-escalation
      x_minus_1_y<-which((dose.indices[,1]==(dose_index[1]-1))&(dose.indices[,2]==dose_index[2]))
      x_y_minus_1<-which((dose.indices[,1]==dose_index[1])&(dose.indices[,2]==(dose_index[2]-1)))
      
      #off-diagonal
      x_plus_1_y_minus_1<-which((dose.indices[,1]==(dose_index[1]+1))&(dose.indices[,2]==(dose_index[2]-1)))
      x_minus_1_y_pluus_1<-which((dose.indices[,1]==(dose_index[1]-1))&(dose.indices[,2]==(dose_index[2]+1)))
      
      #escalation
      x_plus_1_y<-which((dose.indices[,1]==(dose_index[1]+1))&(dose.indices[,2]==dose_index[2]))
      x_y_plus_1<-which((dose.indices[,1]==dose_index[1])&(dose.indices[,2]==(dose_index[2]+1)))
      
      
      
      allowable.under.int[[dose_i]]<-unique(c(x_minus_1_y,x_y_minus_1,
                                              x_plus_1_y_minus_1,x_minus_1_y_pluus_1,
                                              x_plus_1_y,x_y_plus_1,dose_i))
      allowable.in.int[[dose_i]]<-unique(c(x_minus_1_y,x_y_minus_1,
                                           dose_i))
      allowable.over.int[[dose_i]]<-unique(c(x_minus_1_y,x_y_minus_1))
    }
    
    
  }
  
  
  
  
  #define all doses as admissable before any are dropped for safety
  max_admissable<-ndoses
  
  if(hard.safety.rule>50){
    hard.safety<-T
    hard.safety.mat<-hard.safety.mat.function(perc=hard.safety.rule,co.size = co_size,max.cohorts = ncohorts)
  }else{
    hard.safety<-F
  }
  
  
  expand_vec<-rep(0, ndoses)
  dose_rec<-NA
  nextdose<-start.dose
  stop<-0
  
  
  initial<-1
  
  while(stop==0){
    
    
    #DATA GENERATION
    
    if(current.time==0){
      #first cohort
      all.data<-multiple_patient_generation_TTE(patient_ID1=patient_ID1,efficacy_pars =tru.E.pars,tox_pars = tru.T.pars,corET=tru.corET,ncycles=ncycles ,thenextdose=nextdose,
                                                entry_time=current.time,num_patients = co_size)
      current.time<-current.time+1
      patient_ID1<-max(all.data$patient_ID)+1
      current.data<-all.data[all.data$time_of<=current.time,]
      dosevec[current.time]<-nextdose
    }else{
      
      #subsequent cohorts
      #escalation
      all.data<-rbind(all.data,multiple_patient_generation_TTE(patient_ID1=patient_ID1,efficacy_pars =tru.E.pars,tox_pars = tru.T.pars,corET=tru.corET,ncycles=ncycles ,thenextdose=nextdose,
                                                               entry_time=current.time,num_patients = co_size)
      )
      patient_ID1<-max(all.data$patient_ID)+1
      if(backfill==T){
        #expansion
        if(nextdose>1){ #if there are doses below the next dose
          below.doses<-c(1:(nextdose-1)) #which doses are below? 
          expand.below<- expand_vec[below.doses] #have the below doses been expanded?
          
          # which.expand<-which[expand.below==0] #which doses should now be expanded?
          for(expand.doses in below.doses){ #for each dose that needs expanding
            if(expand.below[expand.doses]==0){
              #generate data for 2 cohorts on the dose
              all.data<-rbind(all.data,multiple_patient_generation_TTE(patient_ID1=patient_ID1,efficacy_pars =tru.E.pars,tox_pars = tru.T.pars,corET=tru.corET,ncycles=ncycles ,thenextdose=expand.doses,
                                                                       entry_time=current.time,num_patients = backfill.num*co_size))
              #update patient ID for next assignment
              patient_ID1<-max(all.data$patient_ID)+1
              #update the expansion vector to say this dose has been expanded
              expand_vec[expand.doses]<-1
            }
          }
        }
        
      }
      
      current.time<-current.time+1
      
      current.data<-all.data[all.data$time_of<=current.time,]
      dosevec[current.time]<-nextdose
      
    }
    #  patient.data<-patient_data_frame(all.data)
    
    current.patient.data<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=ncycles)
    
    
    ##TITE CRM has an initial period where we keep escalating until we see a DLT
    
    if(initial==1){
      
      
      
      #next dose is next in the default order
      initial.dose.seq.def<-c()
      initial.dose.seq<-current.patient.data$dose.level
      for(i in 1:length(initial.dose.seq)){
        initial.dose.seq.def[i]<-which(default.order==initial.dose.seq[i])
      }
      nextdose.def<-min(max(initial.dose.seq.def)+1,ndoses)
      nextdose<-default.order[nextdose.def]
      
      
      
      
      ##sufficient information 
      if(sufficient.information==T){
        npats_doses<-sapply(c(1:ndoses), function(x) sum(current.patient.data$dose.level==x))
        num_ass<-tabulate(dosevec,nbins = ndoses) #number of assignments
        
        if(npats_doses[nextdose]>=sufficient.information.lim){
          stop<-4
          stop_vec[4]<-1
          dose_rec<-nextdose
          break
        }
      }
      
      if(initial.one.cycle==T){
        if((sum(current.patient.data$current.DLT)>0)|(sum(current.patient.data$current.Eff)>0)){
          initial<-0
        }
        
      }else{
        
        for (cyc in 1:(ncycles-1)){
          current.patient.data<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=ncycles)
          current.time<-current.time+1
          
          if((sum(current.patient.data$current.DLT)>0)|(sum(current.patient.data$current.Eff)>0)){
            initial<-0
            current.time<-current.time-1
          }
          
        }
      }
      
    }
    if(initial==0){
      #posterior 
      y_E<-current.patient.data$current.Eff
      y_T<-current.patient.data$current.DLT
      patdoses<-current.patient.data$dose.level
      # print(patdoses)
      if(TITE==T){
        weightE<-current.patient.data$patient_weights_eff
        weightT<-current.patient.data$patient_weights_tox
      }else{
        weightE<-rep(1,length(current.patient.data$patient_weights_eff))
        weightT<-rep(1,length(current.patient.data$patient_weights_tox))
      }
      
      
      
      
      prob_over_eff<-prob_under_tox<-P_UB<-pi_T_hat<-pi_E_hat<-c()
      admiss<-eff_obs<-tox_obs<-x_obs<-rep(0,ndoses)
      
      for(dose_i in 1:ndoses){
        
        if(sum(patdoses==dose_i)>0){
          weightT_i<-weightT[patdoses==dose_i]
          weightE_i<-weightE[patdoses==dose_i]
          y_E_i<-y_E[patdoses==dose_i]
          y_T_i<-y_T[patdoses==dose_i]
          
          
          pi_T_hat[dose_i]<-sum((weightT_i==1)*y_T_i)/sum(((weightT_i==1)*y_T_i)+((weightT_i==1)*(1-y_T_i))+(weightT_i*(1-(weightT_i==1))))
          
          pi_E_hat[dose_i]<-sum((weightE_i==1)*y_E_i)/sum(((weightE_i==1)*y_E_i)+((weightE_i==1)*(1-y_E_i))+(weightE_i*(1-(weightE_i==1))))
          part.eff<-part.tox<-part.ut<-pyT1<-pyE1<-c()
          for(pat_i in 1:length(y_T_i)){
            #estimates based on partial obs
            if(weightT_i[pat_i]!=1){
              pyT1[pat_i]<-pi_T_hat[dose_i]*(1-weightT_i[pat_i])/(1-pi_T_hat[dose_i]*weightT_i[pat_i])
            } 
            if(weightE_i[pat_i]!=1){
              pyE1[pat_i]<-pi_E_hat[dose_i]*(1-weightE_i[pat_i])/(1-pi_E_hat[dose_i]*weightE_i[pat_i])
            }
            
            
            
            
            #full utility on full obs
            if((weightE_i[pat_i]==1)&(weightT_i[pat_i]==1)){
              # part.ut[pat_i]<-utility_mat[y_T_i[pat_i]+1,y_E_i[pat_i]+1]
              part.ut[pat_i]<-y_E_i[pat_i]+WT*(1-y_T_i[pat_i])
            }
            #partial utility on partial obs
            if((weightE_i[pat_i]!=1)&(weightT_i[pat_i]==1)){
              # part.ut[pat_i]<-sum(matrix(c((1-y_T_i[pat_i])*(1-pyE1[pat_i]),
              #                               (y_T_i[pat_i])*(1-pyE1[pat_i]),
              #                               (1-y_T_i[pat_i])*(pyE1[pat_i]),
              #                               (y_T_i[pat_i])*(pyE1[pat_i])),nrow=2)*utility_mat)
              part.ut[pat_i]<-pyE1[pat_i]+WT*(1-y_T_i[pat_i])
              
            }
            
            if((weightE_i[pat_i]==1)&(weightT_i[pat_i]!=1)){
              # part.ut[pat_i]<-sum(matrix(c((1-pyT1[pat_i])*(1-y_E_i[pat_i]),
              #                              (pyT1[pat_i])*(1-y_E_i[pat_i]),
              #                              (1-pyT1[pat_i])*(y_E_i[pat_i]),
              #                              (pyT1[pat_i])*(y_E_i[pat_i])),nrow=2)*utility_mat)
              part.ut[pat_i]<-y_E_i[pat_i]+WT*(1-pyT1[pat_i])
              
            }
            if((weightE_i[pat_i]!=1)&(weightT_i[pat_i]!=1)){
              # part.ut[pat_i]<-sum(matrix(c((1-pyT1[pat_i])*(1-pyE1[pat_i]),
              #                              (pyT1[pat_i])*(1-pyE1[pat_i]),
              #                              (1-pyT1[pat_i])*(pyE1[pat_i]),
              #                              (pyT1[pat_i])*(pyE1[pat_i])),nrow=2)*utility_mat)
              part.ut[pat_i]<-pyE1[pat_i]+WT*(1-pyT1[pat_i])
              
            }
            
            ####safety only
            #full utility on full obs
            if((weightT_i[pat_i]==1)){
              part.tox[pat_i]<-y_T_i[pat_i]
            }
            #partial utility on partial obs
            if((weightT_i[pat_i]!=1)){
              part.tox[pat_i]<-pyT1[pat_i]
            }
            
            ####activity only
            #full utility on full obs
            if((weightE_i[pat_i]==1)){
              part.eff[pat_i]<-y_E_i[pat_i]
            }
            #partial utility on partial obs
            if((weightE_i[pat_i]!=1)){
              part.eff[pat_i]<-pyE1[pat_i]
            }  
            
          } # for pat_i
          
          
          
          
          #based on follow up observed (not event obs)
          x_obs[dose_i]<-sum(part.ut)/(1+WT)
          tox_obs[dose_i]<-sum(part.tox)
          prob_under_tox[dose_i]<-pbeta(toxbound,1+tox_obs[dose_i],1+length(y_T_i)-tox_obs[dose_i])
          
          eff_obs[dose_i]<-sum(part.eff)
          prob_over_eff[dose_i]<-1-pbeta(effbound,1+eff_obs[dose_i],1+length(y_E_i)-eff_obs[dose_i])
          
          
          #P(utility>benchmark)
          P_UB[dose_i]<-1-pbeta(u_b,1+x_obs[dose_i],1+length(y_T_i)-x_obs[dose_i])
          
          
        }else{  # ELSE if dose_i has no patients
          P_UB[dose_i]<-1-pbeta(u_b,1,1)
          prob_over_eff[dose_i]<-1-pbeta(effbound,1,1)
          prob_under_tox[dose_i]<-pbeta(toxbound,1,1)
        }
        
        #define admissible
        
        admiss[dose_i]<-  (prob_over_eff[dose_i]>C_eff)&(prob_under_tox[dose_i]>C_tox)
        
        
        
      } # for dose_i
      admiss.doses<-which(admiss==1)
      
      current_dose<-dosevec[current.time]
      #browser()
      #current dose: is pitT under/in/over interval?
      #if safety above interval: de-escalation
      if(pi_T_hat[current_dose]>delta_U){
        allowable_next<-allowable.over.int[[current_dose]]
        
      }
      
      #if safety in interval: stay or de-escalate
      if((pi_T_hat[current_dose]>=delta_L)&(pi_T_hat[current_dose]<=delta_U)){
        allowable_next<-allowable.in.int[[current_dose]]
        
      }
      #if safety below interval: de-escalate, stay or escalate
      if(pi_T_hat[current_dose]<delta_L){
        allowable_next<-allowable.under.int[[current_dose]]
        
      }
      #define set of doses to go to (this includes are they admissable?)
      allowable_next_adm2<-(c(1:ndoses)%in%allowable_next)+(c(1:ndoses)%in%admiss.doses)
      allowable_next_adm<-which(allowable_next_adm2==2)
      
      #choose admissible dose with max prob of being over benchmark utility
      if(length(allowable_next_adm)>0){ #if there are doses admissiable in the set
        
        #maximises prob of being greater than benchmark
        nextdose<-allowable_next_adm[which(P_UB[allowable_next_adm]==max(P_UB[allowable_next_adm]))]
        
        #choose randomly if a draw
        if(length(nextdose)>1){
          nextdose<-sample(nextdose)[1]
        }
      }else{ ##if no admissable doses
        
        if(pause==T){
          pause.time<-current.time
          while((length(allowable_next_adm)==0)&(current.time<pause.time+ncycles-1)){
            current.time<-current.time+1
            current.data<-all.data[all.data$time_of<=current.time,]
            
            current.patient.data<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=ncycles)
            
            
            
            
            y_E<-current.patient.data$current.Eff
            y_T<-current.patient.data$current.DLT
            
            patdoses<-current.patient.data$dose.level 
            
            if(TITE==T){
              weightE<-current.patient.data$patient_weights_eff
              weightT<-current.patient.data$patient_weights_tox
            }else{
              weightE<-rep(1,length(current.patient.data$patient_weights_eff))
              weightT<-rep(1,length(current.patient.data$patient_weights_tox))
            }
            
            
            
            prob_over_eff<-prob_under_tox<-P_UB<-pi_T_hat<-pi_E_hat<-c()
            admiss<-eff_obs<-tox_obs<-x_obs<-rep(0,ndoses)
            
            for(dose_i in 1:ndoses){
              
              if(sum(patdoses==dose_i)>0){
                weightT_i<-weightT[patdoses==dose_i]
                weightE_i<-weightE[patdoses==dose_i]
                y_E_i<-y_E[patdoses==dose_i]
                y_T_i<-y_T[patdoses==dose_i]
                
                
                pi_T_hat[dose_i]<-sum((weightT_i==1)*y_T_i)/sum(((weightT_i==1)*y_T_i)+((weightT_i==1)*(1-y_T_i))+(weightT_i*(1-(weightT_i==1))))
                
                pi_E_hat[dose_i]<-sum((weightE_i==1)*y_E_i)/sum(((weightE_i==1)*y_E_i)+((weightE_i==1)*(1-y_E_i))+(weightE_i*(1-(weightE_i==1))))
                part.eff<-part.tox<-part.ut<-pyT1<-pyE1<-c()
                for(pat_i in 1:length(y_T_i)){
                  #estimates based on partial obs
                  if(weightT_i[pat_i]!=1){
                    pyT1[pat_i]<-pi_T_hat[dose_i]*(1-weightT_i[pat_i])/(1-pi_T_hat[dose_i]*weightT_i[pat_i])
                  } 
                  if(weightE_i[pat_i]!=1){
                    pyE1[pat_i]<-pi_E_hat[dose_i]*(1-weightE_i[pat_i])/(1-pi_E_hat[dose_i]*weightE_i[pat_i])
                  }
                  #full utility on full obs
                  if((weightE_i[pat_i]==1)&(weightT_i[pat_i]==1)){
                    part.ut[pat_i]<-utility_mat[y_T_i[pat_i]+1,y_E_i[pat_i]+1]
                  }
                  #partial utility on partial obs
                  if((weightE_i[pat_i]!=1)&(weightT_i[pat_i]==1)){
                    part.ut[pat_i]<-sum(matrix(c((1-y_T_i[pat_i])*(1-pyE1[pat_i]),
                                                 (y_T_i[pat_i])*(1-pyE1[pat_i]),
                                                 (1-y_T_i[pat_i])*(pyE1[pat_i]),
                                                 (y_T_i[pat_i])*(pyE1[pat_i])),nrow=2)*utility_mat)
                  }
                  
                  if((weightE_i[pat_i]==1)&(weightT_i[pat_i]!=1)){
                    part.ut[pat_i]<-sum(matrix(c((1-pyT1[pat_i])*(1-y_E_i[pat_i]),
                                                 (pyT1[pat_i])*(1-y_E_i[pat_i]),
                                                 (1-pyT1[pat_i])*(y_E_i[pat_i]),
                                                 (pyT1[pat_i])*(y_E_i[pat_i])),nrow=2)*utility_mat)
                  }
                  if((weightE_i[pat_i]!=1)&(weightT_i[pat_i]!=1)){
                    part.ut[pat_i]<-sum(matrix(c((1-pyT1[pat_i])*(1-pyE1[pat_i]),
                                                 (pyT1[pat_i])*(1-pyE1[pat_i]),
                                                 (1-pyT1[pat_i])*(pyE1[pat_i]),
                                                 (pyT1[pat_i])*(pyE1[pat_i])),nrow=2)*utility_mat)
                  }
                  
                  ####safety only
                  #full utility on full obs
                  if((weightT_i[pat_i]==1)){
                    part.tox[pat_i]<-y_T_i[pat_i]
                  }
                  #partial utility on partial obs
                  if((weightT_i[pat_i]!=1)){
                    part.tox[pat_i]<-pyT1[pat_i]
                  }
                  
                  ####activity only
                  #full utility on full obs
                  if((weightE_i[pat_i]==1)){
                    part.eff[pat_i]<-y_E_i[pat_i]
                  }
                  #partial utility on partial obs
                  if((weightE_i[pat_i]!=1)){
                    part.eff[pat_i]<-pyE1[pat_i]
                  }  
                  
                } # for pat_i
                
                
                
                
                #based on follow up observed (not event obs)
                x_obs[dose_i]<-0.01*sum(part.ut) 
                tox_obs[dose_i]<-sum(part.tox)
                prob_under_tox[dose_i]<-pbeta(toxbound,1+tox_obs[dose_i],1+length(y_T_i)-tox_obs[dose_i])
                
                eff_obs[dose_i]<-sum(part.eff)
                prob_over_eff[dose_i]<-1-pbeta(effbound,1+eff_obs[dose_i],1+length(y_E_i)-eff_obs[dose_i])
                
                
                #P(utility>benchmark)
                P_UB[dose_i]<-1-pbeta(u_b,1+x_obs[dose_i],1+length(y_T_i)-x_obs[dose_i])
                
                
              }else{  # ELSE if dose_i has no patients
                P_UB[dose_i]<-1-pbeta(u_b,1,1)
                prob_over_eff[dose_i]<-1-pbeta(effbound,1,1)
                prob_under_tox[dose_i]<-pbeta(toxbound,1,1)
              }
              
              #define admissible
              
              admiss[dose_i]<-  (prob_over_eff[dose_i]>C_eff)&(prob_under_tox[dose_i]>C_tox)
              
              
              
            } # for dose_i
            admiss.doses<-which(admiss==1)
            
            current_dose<-dosevec[pause.time]
            #browser()
            #current dose: is pitT under/in/over interval?
            #if safety above interval: de-escalation
            if(pi_T_hat[current_dose]>delta_U){
              allowable_next<-allowable.over.int[[current_dose]]
              
            }
            
            #if safety in interval: stay or de-escalate
            if((pi_T_hat[current_dose]>=delta_L)&(pi_T_hat[current_dose]<=delta_U)){
              allowable_next<-allowable.in.int[[current_dose]]
              
            }
            #if safety below interval: de-escalate, stay or escalate
            if(pi_T_hat[current_dose]<delta_L){
              allowable_next<-allowable.under.int[[current_dose]]
              
            }
            #define set of doses to go to (this includes are they admissable?)
            allowable_next_adm2<-(c(1:ndoses)%in%allowable_next)+(c(1:ndoses)%in%admiss.doses)
            allowable_next_adm<-which(allowable_next_adm2==2)
            
            #choose admissible dose with max prob of being over benchmark utility
            if(length(allowable_next_adm)>0){ #if there are doses admissiable in the set
              
              #maximises prob of being greater than benchmark
              nextdose<-allowable_next_adm[which(P_UB[allowable_next_adm]==max(P_UB[allowable_next_adm]))]
              
              #choose randomly if a draw
              if(length(nextdose)>1){
                nextdose<-sample(nextdose)[1]
              }
              
              
              
              
              
            }
          }#for while
          
          
          
        }#for pause
        if(length(allowable_next_adm)==0){
          
          if(length(y_T)>=a.stop.bound){
            stop<-1
            stop_vec[1]<-1
            dose_rec<-NA
            break
          }else{
            #choose admissibility based on tox only
            
            for(dose_i in 1:ndoses){
              admiss[dose_i]<-  (prob_under_tox[dose_i]>C_tox)
              
              
              
            } # for dose_i
            admiss.doses<-which(admiss==1)
            
            current_dose<-dosevec[pause.time]
            #browser()
            #current dose: is pitT under/in/over interval?
            #if safety above interval: de-escalation
            if(pi_T_hat[current_dose]>delta_U){
              allowable_next<-allowable.over.int[[current_dose]]
              
            }
            
            #if safety in interval: stay or de-escalate
            if((pi_T_hat[current_dose]>=delta_L)&(pi_T_hat[current_dose]<=delta_U)){
              allowable_next<-allowable.in.int[[current_dose]]
              
            }
            #if safety below interval: de-escalate, stay or escalate
            if(pi_T_hat[current_dose]<delta_L){
              allowable_next<-allowable.under.int[[current_dose]]
              
            }
            #define set of doses to go to (this includes are they admissable?)
            allowable_next_adm2<-(c(1:ndoses)%in%allowable_next)+(c(1:ndoses)%in%admiss.doses)
            allowable_next_adm<-which(allowable_next_adm2==2)
            
            #choose admissible dose with max prob of being over benchmark utility
            if(length(allowable_next_adm)>0){ #if there are doses admissiable in the set
              
              #maximises prob of being greater than benchmark
              nextdose<-allowable_next_adm[which(P_UB[allowable_next_adm]==max(P_UB[allowable_next_adm]))]
              
              #choose randomly if a draw
              if(length(nextdose)>1){
                nextdose<-sample(nextdose)[1]
              }
              
            }
            
            if(length(allowable_next_adm)==0){
              stop<-1
              stop_vec[1]<-1
              dose_rec<-NA
              break
            }
            
          }
        }
      }#for if no admissible
    }
    #num pats
    npats_doses<-sapply(c(1:ndoses), function(x) sum(current.patient.data$dose.level==x))
    #number of DLTs per dose level
    nDLTs_doses<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)]))
    nDLTs_doses_c1<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)&(current.data$cycle==1)]))
    
    
    #stopping rules
    if(((safety.stopping.low.unsafe==T)&(npats_doses[1]>0))|((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0))){
      
      
      current.patient.data_cyc1<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=1)
      
      
      weightE1<-weightT1<-c()
      patdoses1<-current.patient.data_cyc1$dose.level
      
      
      #
      if(TITE==T){
        
        for(obser in 1:nrow(current.patient.data_cyc1)){
          weightE1[obser]<-current.patient.data_cyc1$patient_weights_eff[obser]
          weightT1[obser]<-current.patient.data_cyc1$patient_weights_tox[obser]
          
        }
      }else{
        
        weightE1<-rep(1,nrow(current.patient.data_cyc1))
        weightT1<-rep(1,nrow(current.patient.data_cyc1))
        
      }
      
      
      
      
      
      
      #only cycle 1 
      
      if((safety.stopping.low.unsafe==T)&(npats_doses[1]>0)){
        nDLTs1<- nDLTs_doses_c1[1]
        cyc1_0.3g<-1-pbeta(0.3,1+nDLTs1,1+npats_doses[1]-nDLTs1)
        
        
        if(cyc1_0.3g>0.8){
          stop<-6
          stop_vec[6]<-1
          nextdose<-NA
          
        }
      }
      
      if((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0)){
        
        nDLTsN<- nDLTs_doses_c1[ndoses]
        
        cycJ_0.3l<-pbeta(0.3,1+nDLTsN,1+npats_doses[ndoses]-nDLTsN)
        
        if(cycJ_0.3l>0.8){
          
          stop<-7
          stop_vec[7]<-1
          nextdose<-NA
          
        }
      }
      
    }
    
    
    
    
    
    
    ##hard safety
    if(hard.safety==T){
      
      explored<-c(1:ndoses)[npats_doses>0]
      
      for(do in explored){
        
        if(nDLTs_doses_c1[do]>=hard.safety.mat[1, which(hard.safety.mat[2,]==npats_doses[do])]){
          
          min.over.indices.hold<-dose.indices[do,]
          excluded.hold<-as.numeric((dose.indices[1]>=min.over.indices.hold[1])&(dose.indices[2]>=min.over.indices.hold[2]))
          excluded<-excluded+excluded.hold
          
          
        }
      }
      if(!is.na(nextdose)){
        if(excluded[nextdose]>0){#if the next dose is in the exclusion zone
          
          #maximises prob of being greater than benchmark
          nextdose<-allowable_next_adm[which(P_UB[allowable_next_adm]==max(P_UB[allowable_next_adm]))]
          
          
          included.P_UB<-P_UB
          included.P_UB[excluded>0]<--100
          included.P_UB[-allowable_next_adm]<--100
          nextdose<-which.max(included.P_UB)
        }
      }
      if(all(excluded>0)){
        stop<-5
        stop_vec[5]<-1
        dose_rec<-NA
        nextdose<-NA
        
      }
      
      
    }
    
    
    
    
    ##sufficient information 
    if((sufficient.information==T)&(!is.na(nextdose))){
      if(npats_doses[nextdose]>=sufficient.information.lim){
        
        
        
        if(prob_over_eff[dose_i]<C_eff){
          stop<-1
          stop_vec[c(1,4)]<-1
          dose_rec<-NA
          break
        }else{
          stop<-4
          stop_vec[4]<-1
          dose_rec<-nextdose
        }
        
        
        
      }
    }
    
    
    
    
    if(nrow(current.patient.data)==(ncohorts*co_size)){ #max patients reached
      ##choose final rec
      
      current.patient.data<-current_patient_data_frame(current_time=(current.time+ncycles),patient.dataframe=all.data,follow_up = ncycles)
      patdoses<-current.patient.data$dose.level
      y_E<-current.patient.data$current.Eff
      y_T<-current.patient.data$current.DLT
      
      if(TITE==T){
        weightE<-current.patient.data$patient_weights_eff
        weightT<-current.patient.data$patient_weights_tox
      }else{
        weightE<-rep(1,length(current.patient.data$patient_weights_eff))
        weightT<-rep(1,length(current.patient.data$patient_weights_tox))
      }
      
      ###########################
      
      
      prob_over_eff<-prob_under_tox<-P_UB<-pi_T_hat<-pi_E_hat<-c()
      admiss<-eff_obs<-tox_obs<-x_obs<-rep(0,ndoses)
      
      for(dose_i in 1:ndoses){
        
        if(sum(patdoses==dose_i)>0){
          weightT_i<-weightT[patdoses==dose_i]
          weightE_i<-weightE[patdoses==dose_i]
          y_E_i<-y_E[patdoses==dose_i]
          y_T_i<-y_T[patdoses==dose_i]
          
          
          pi_T_hat[dose_i]<-sum((weightT_i==1)*y_T_i)/sum(((weightT_i==1)*y_T_i)+((weightT_i==1)*(1-y_T_i))+(weightT_i*(1-(weightT_i==1))))
          
          pi_E_hat[dose_i]<-sum((weightE_i==1)*y_E_i)/sum(((weightE_i==1)*y_E_i)+((weightE_i==1)*(1-y_E_i))+(weightE_i*(1-(weightE_i==1))))
          part.eff<-part.tox<-part.ut<-pyT1<-pyE1<-c()
          for(pat_i in 1:length(y_T_i)){
            #estimates based on partial obs
            if(weightT_i[pat_i]!=1){
              pyT1[pat_i]<-pi_T_hat[dose_i]*(1-weightT_i[pat_i])/(1-pi_T_hat[dose_i]*weightT_i[pat_i])
            } 
            if(weightE_i[pat_i]!=1){
              pyE1[pat_i]<-pi_E_hat[dose_i]*(1-weightE_i[pat_i])/(1-pi_E_hat[dose_i]*weightE_i[pat_i])
            }
            #full utility on full obs
            if((weightE_i[pat_i]==1)&(weightT_i[pat_i]==1)){
              # part.ut[pat_i]<-utility_mat[y_T_i[pat_i]+1,y_E_i[pat_i]+1]
              part.ut[pat_i]<-y_E_i[pat_i]+WT*(1-y_T_i[pat_i])
            }
            #partial utility on partial obs
            if((weightE_i[pat_i]!=1)&(weightT_i[pat_i]==1)){
              # part.ut[pat_i]<-sum(matrix(c((1-y_T_i[pat_i])*(1-pyE1[pat_i]),
              #                               (y_T_i[pat_i])*(1-pyE1[pat_i]),
              #                               (1-y_T_i[pat_i])*(pyE1[pat_i]),
              #                               (y_T_i[pat_i])*(pyE1[pat_i])),nrow=2)*utility_mat)
              part.ut[pat_i]<-pyE1[pat_i]+WT*(1-y_T_i[pat_i])
              
            }
            
            if((weightE_i[pat_i]==1)&(weightT_i[pat_i]!=1)){
              # part.ut[pat_i]<-sum(matrix(c((1-pyT1[pat_i])*(1-y_E_i[pat_i]),
              #                              (pyT1[pat_i])*(1-y_E_i[pat_i]),
              #                              (1-pyT1[pat_i])*(y_E_i[pat_i]),
              #                              (pyT1[pat_i])*(y_E_i[pat_i])),nrow=2)*utility_mat)
              part.ut[pat_i]<-y_E_i[pat_i]+WT*(1-pyT1[pat_i])
              
            }
            if((weightE_i[pat_i]!=1)&(weightT_i[pat_i]!=1)){
              # part.ut[pat_i]<-sum(matrix(c((1-pyT1[pat_i])*(1-pyE1[pat_i]),
              #                              (pyT1[pat_i])*(1-pyE1[pat_i]),
              #                              (1-pyT1[pat_i])*(pyE1[pat_i]),
              #                              (pyT1[pat_i])*(pyE1[pat_i])),nrow=2)*utility_mat)
              part.ut[pat_i]<-pyE1[pat_i]+WT*(1-pyT1[pat_i])
              
            }
            
            ####safety only
            #full utility on full obs
            if((weightT_i[pat_i]==1)){
              part.tox[pat_i]<-y_T_i[pat_i]
            }
            #partial utility on partial obs
            if((weightT_i[pat_i]!=1)){
              part.tox[pat_i]<-pyT1[pat_i]
            }
            
            ####activity only
            #full utility on full obs
            if((weightE_i[pat_i]==1)){
              part.eff[pat_i]<-y_E_i[pat_i]
            }
            #partial utility on partial obs
            if((weightE_i[pat_i]!=1)){
              part.eff[pat_i]<-pyE1[pat_i]
            }  
            
          } # for pat_i
          
          
          
          
          #based on follow up observed (not event obs)
          x_obs[dose_i]<-0.01*sum(part.ut) 
          tox_obs[dose_i]<-sum(part.tox)
          prob_under_tox[dose_i]<-pbeta(toxbound,1+tox_obs[dose_i],1+length(y_T_i)-tox_obs[dose_i])
          
          eff_obs[dose_i]<-sum(part.eff)
          prob_over_eff[dose_i]<-1-pbeta(effbound,1+eff_obs[dose_i],1+length(y_E_i)-eff_obs[dose_i])
          
          
          #P(utility>benchmark)
          P_UB[dose_i]<-1-pbeta(u_b,1+x_obs[dose_i],1+length(y_T_i)-x_obs[dose_i])
          
          
        }else{  # ELSE if dose_i has patients
          P_UB[dose_i]<-1-pbeta(u_b,1,1)
          prob_over_eff[dose_i]<-1-pbeta(effbound,1,1)
          prob_under_tox[dose_i]<-pbeta(toxbound,1,1)
        }
        
        #define admissible
        
        admiss[dose_i]<-  (prob_over_eff[dose_i]>C_eff)&(prob_under_tox[dose_i]>C_tox)
        
        
        
      } # for dose_i
      admiss.doses<-which(admiss==1)
      #final rec
      final_safe<-which(pi_T_hat<=target)
      allowable_final2<-(c(1:ndoses)%in%final_safe)+(c(1:ndoses)%in%admiss.doses)
      allowable_final<-which(allowable_final2==2)
      
      #choose admissible dose with max prob of being over benchmark utility
      if(length(allowable_final)>0){ #if there are doses admissiable in the set
        
        #maximises prob of being greater than benchmark
        dose_rec<-allowable_final[which(P_UB[allowable_final]==max(P_UB[allowable_final]))]
        
        #choose randomly if a draw
        if(length(dose_rec)>1){
          dose_rec<-sample(dose_rec)[1]
        }
        
        stop<-3
        stop_vec[3]<-1
        
      }else{ 
        stop<-1
        stop_vec[1]<-1
        dose_rec<-NA
        break
      }    ##if no admissable doses
      
      ##########################
      
      
    }
    
  }
  #follow up for all patients
  current.data<-all.data
  current.time<-max(current.data$time_of)
  current.patient.data<-current_patient_data_frame(current_time = current.time,patient.dataframe = all.data,follow_up = ncycles)
  
  
  
  #grade matrix out
  DLT.matrix.out<-matrix(0,ncol=ndoses,nrow=2)
  for (pat in 1:(max(current.data$patient))){
    pat_max_grade<-max(current.data[current.data$patient_ID==pat,]$DLT)
    pat_dose<-current.data[current.data$patient_ID==pat,]$dose_level[1]
    DLT.matrix.out[pat_max_grade+1,pat_dose]<- DLT.matrix.out[pat_max_grade+1,pat_dose]+1
  }
  
  #grade matrix out
  Eff.matrix.out<-matrix(0,ncol=ndoses,nrow=2)
  for (pat in 1:(max(current.data$patient))){
    pat_max_grade<-max(current.data[current.data$patient_ID==pat,]$Eff)
    pat_dose<-current.data[current.data$patient_ID==pat,]$dose_level[1]
    Eff.matrix.out[pat_max_grade+1,pat_dose]<- Eff.matrix.out[pat_max_grade+1,pat_dose]+1
  }
  
  #MTD.rec<-min(c(max(dosevec),max_admissable,which.min(abs(target-expit(mean(gibbs_out$alphaT)+doses*mean(gibbs_out$betaT))))))
  
  
  #output list
  output<-list(
    dose.rec=dose_rec ,#dose recommendation
    num.pat=max(current.data$patient) ,#: number of patients
    dose.ass=tabulate(dosevec,nbins = ndoses) ,# : number of cohorts per dose (vector)
    stop.code=stop_vec, #: stopping reason 
    num.DLT= sum(current.data$DLT[current.data$patient_ID>0]),#: number of DLTS (total)
    DLT.mat=DLT.matrix.out ,#: 2xnumdoses matrix. rows for grades 1=no DLT, 2=DLT. cols for doses.
    Eff.mat=Eff.matrix.out ,#: 2xnumdoses matrix. rows for grades 1=noEff, 2=Eff. cols for doses.
    duration= current.time,#: total trial duration (until all recruited patients are fully observed)
    max_admissable=max_admissable, # max admissable dose at the end of the trial (has hard safety eliminated any?)
    dosevec=dosevec,
    DLT.vec=current.patient.data$current.DLT,
    EFF.vec=current.patient.data$current.Eff,
    all.data=current.patient.data
  )
  
  return(output)
  
}

