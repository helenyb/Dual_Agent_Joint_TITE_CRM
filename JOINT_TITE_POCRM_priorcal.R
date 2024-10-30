#JOINT TITE-POCRM function function tox only (for calibration)

#part 2 is 2 parameter joint
#part 1 has 1 parameter independent


library(rjags)
library(cubature)
library(mvtnorm)
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


#posterior as function of parameters, for integrations (independent, 1 par)
post.fun.ind.1par<-function(parameters,y,patdoses,wm_skeleton,weight,prior_hyp){
  ##parameter
  alpha<-parameters
  
  #prior hyper-parameters
  
  mu_alpha<-prior_hyp[1]
  tau_alpha<-prior_hyp[2]
  
  
  
  G<-c()
  pi<-c()
  ndoses<-length(wm_skeleton)
  for(j in 1:ndoses){
    pi[j] <- (wm_skeleton[j])^exp(alpha)
    
  }
  
  
  for(i in 1:length(patdoses)){
    G[i]<-weight[i]*pi[patdoses[i]]
    
  }
  
  prod((G^y)*((1-G)^(1-y)))*
    dnorm(alpha, mean=mu_alpha,sd=sqrt(1/tau_alpha))
  
}
#vectorises for integration function
V.post.fun.ind.1par<-function(parameters,y,patdoses,wm_skeleton,weight,prior_hyp){
  sapply(parameters,post.fun.ind.1par,y,patdoses,wm_skeleton,weight,prior_hyp)
  
}



#get working model
#orders is a matrix, each row is a different ordering
getwm<-function (orders, skeleton) 
{
  d <- ncol(orders) #num of dose combinations
  s <- nrow(orders) #num of orderings
  alpha <- matrix(0, nrow = s, ncol = d)
  for (j in 1:s) {
    alpha[j, ] <- skeleton[order(orders[j, ])]
  }
  alpha
}


#model for rjags
model.crm.joint.string <-"
model {




for(j in 1:ndoses){

logit(piT[j]) = alphaT + betaT*dosesT[j]

}

for(i in 1:length(yDLT)){

G_piT[i]=weightT[i]*piT[patdoses[i]]
yDLT[i] ~ dbinom(G_piT[i],1)

}



log_betaT ~ dnorm(mu_betaT, tau_betaT)
betaT = exp(log_betaT)
alphaT ~ dnorm(mu_alphaT, tau_alphaT)


}
"

##function for gibbs sampler
gibbs_sampler.crm.combo<-function(data_list,iter,model_string){
  model1.spec<-textConnection(model_string)
  jags <- jags.model(model1.spec,data =data_list,n.chains=1,n.adapt=1000,quiet=T)
  update(jags, 1000,progress.bar="none")
  tt<-jags.samples(jags,c('alphaT','betaT'),iter,progress.bar="none")
  return(tt)
}


#defines category



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


#function that takes the allowables list, the vector of doses assigned so far, and gives all allowable doses
allowable_func<-function(allowable.list.in,doses.so.far){
  
  all.allowable<-unique(unlist(allowable.list.in[doses.so.far]))
  return(all.allowable)
}


##INPUT:
#seed=seed for reproducibility
#tru.E.pars=parameters for data generation of efficacy times 
#tru.T.pars=parameters for data generation of DLT times 
#tru.corET=correlation between DLT times and efficacy times in data generation
#co_size=cohort size
#ncohorts=maximum number of cohorts in study
#target=target toxicity (all cycles)
#ncycles=number of cycles
#dose.skipping.rule= type of dose.skipping rule:
#"OFF" is no restrictions
#"ON" is once dose(a,b) is explored, all doses(x,y) such that x<=a+1 and y<=b are allowable, and x<=a and y<=b+1
#"ON.diag.allowed" is same as "ON" but also allows (a+1,b+1)
#prior_vec1: prior hyperparameters for 1 parameter independent model
#prior_vec1[1]: mean for toxicity
#prior_vec1[2]: precision for toxicity
#prior_vec2: prior hyperparameters for joint 2-parameter model
#mu_alphaT=prior_vec2[1] prior mean for a_T (intercept, toxicity)
#tau_alphaT=prior_vec2[2] prior precision for a_T (intercept, toxicity)
#mu_betaT=prior_vec2[3] prior mean for b_T (slope, toxicity)
#tau_betaT=prior_vec2[4] prior precision for b_T (slope, toxicity)
#sufficient.information==enforce stopping for sufficient information? no more than 9 patients per dose. T=enforce stopping for sufficient information (default)
#sufficient.information.lim= number of patients needed before stopping when the next assignment is the same. 
#hard.safety.rule=percentage for hard safety rule based on Beta(1.1)? 
#e.g. : 85= 2/3,3/6,4/9. 90=2/3,4/6,5/9, 95=3/3,4/6,5/9, <50 means no hard safety enforced
#safety.stopping.low.unsafe= (T=stop when P(p1>0.3)>0.8 (cycle 1))
#safety.stopping.high.toosafe= (T= stop when P(pJ>0.3)>0.8 (cycle 1))
#initial.one.cycle: Is the initial period based on one cycle at a time? (F=wait until all cycles completed before next dose in initial period)
#C_tox: Dose is "admissible in safety if P(P(DLT)<toxbound)>C_tox
#toxbound: Dose is "admissible in safety if P(P(DLT)<toxbound)>C_tox
#backfill: should doses deemed safe be backfilled? (Default FALSE) *NOTE backfill only uses default order*
#backfill.num: How many cohorts to add as backfilling?
#TITE: Is TITE-CRM used? (Default=T, F implies CRM only and is only compatible with ncycles=1)
#dose.indices: matrix with 2 col and ndoses rows, of indices for doses (i.e row 2 tells us dose "2" is [2,1] on the grid)
#wmT: matrix of working model (skeletons) for toxicity
#prior.o.T: prior probability of the orderings for toxicity
#default.order: Order in which we escalate if we see no activity/toxicity
#start.dose: The starting dose, as defined 


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
# all.data: All of the trial data

JOINT.TITE.POCRM.priorcal<-function(seed,tru.E.pars,tru.T.pars,tru.corET,co_size,ncohorts ,target,
                              ncycles,dose.skipping.rule, prior_vec1, prior_vec2,
                              sufficient.information=T,sufficient.information.lim=29,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                              safety.stopping.high.toosafe=T,initial.one.cycle=T,
                              C_tox=0.2,toxbound=0.3,
                              backfill=F,backfill.num=2,TITE=T,dose.indices,
                              wmT,prior.o.T,default.order,start.dose){
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
  
  
  
  
  if(!dose.skipping.rule%in%c("ON","ON.diag.allowed","OFF")){
    stop("Dose skipping rule undefined")
  }
  
  
  ##definitions of allowable escalations for dose.skipping.rule
  #define the matrix
  dose.mat<-matrix(nrow=max(dose.indices[,1]),ncol=max(dose.indices[,2]))
  for(dose_i in 1:nrow(dose.indices)){
    dose.mat[dose.indices[dose_i,1],dose.indices[dose_i,2]]<-dose_i
  }
  
  #define the allowable.escalations list
  if(dose.skipping.rule=="ON"){
    allowable.escalations<-list()
    for(dose_i in 1:nrow(dose.indices)){
      
      dose_index<-dose.indices[dose_i,]
      all.below<-which((dose.indices[,1]<=dose_index[1])&(dose.indices[,2]<=dose_index[2]))
      a_plus_1<-which((dose.indices[,1]<=(dose_index[1]+1))&(dose.indices[,2]<=dose_index[2]))
      b_plus_1<-which((dose.indices[,1]<=dose_index[1])&(dose.indices[,2]<=(dose_index[2])+1))
      
      allowable.escalations[[dose_i]]<-unique(c(all.below,a_plus_1,b_plus_1))
    }
    
  }else if(dose.skipping.rule=="ON.diag.allowed"){
    allowable.escalations<-list()
    for(dose_i in 1:nrow(dose.indices)){
      
      dose_index<-dose.indices[dose_i,]
      
      allowable.escalations[[dose_i]]<-which((dose.indices[,1]<=(dose_index[1]+1))&(dose.indices[,2]<=(dose_index[2])+1))
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
      
      patdoses<-current.patient.data$dose.level
      # print(patdoses)
      if(TITE==T){
        weightE<-current.patient.data$patient_weights_eff
        weightT<-current.patient.data$patient_weights_tox
      }else{
        weightE<-rep(1,length(current.patient.data$patient_weights_eff))
        weightT<-rep(1,length(current.patient.data$patient_weights_tox))
      }
      #formatting for gibbs sampler
      #yOUT 
      #1= no E no T
      #2= no E yes T
      #3= yes E no T
      #4= yes E yes T
      y_E<-current.patient.data$current.Eff
      y_T<-current.patient.data$current.DLT
      ########################################  
      ##decisions##
      ###################################  
      
      
      
      ######### 1 par ind
      
      margT1<-c()
      
      
      
      
      
      
      
      for (kT in 1:length(prior.o.T)) { 
        # browser()
        
        margT1[kT]<-integrate(V.post.fun.ind.1par,lower=-20,upper=20,y=y_T,patdoses=patdoses,wm_skeleton=wmT[kT,],
                              weight=weightT,prior_hyp=prior_vec1[1:2])$value
        
        
        #print(kT)
      } # for kT
      
      #browser()
      #posterior
      
      pordT1 <- (margT1 * prior.o.T)/sum(margT1 * prior.o.T)
      #maximise the posterior - which ordering?
      
      ordT1.I<-which.is.max.func(pordT1)
      
      
      
      
      
      
      #calculate posterior for that ordering
      currentdata_pocrm1.I<-list(ndoses=ncol(wmT),weightT=weightT,yDLT=y_T,
                                 dosesT=wmT[ordT1.I,],patdoses=patdoses,
                                 mu_alphaT=prior_vec2[1],mu_betaT=prior_vec2[3],
                                 tau_alphaT=prior_vec2[2],tau_betaT=prior_vec2[4])
      
      #browser()
      
      
      gibbs_out_pocrm1.I<-gibbs_sampler.crm.combo(data_list=currentdata_pocrm1.I,iter=10000,model_string=model.crm.joint.string)
      #browser()
      #choose next dose
      
      
      #R value for the best ordering
      rpredT1.I <- expit(mean(gibbs_out_pocrm1.I$alphaT)+wmT[ordT1.I,]*mean(gibbs_out_pocrm1.I$betaT))
      
      #browser()
      
      # choose the next dose level
      
      #admissable doses
      admiss<-probtox<-c()
      for(dose.level in 1:ndoses){
        probtox[dose.level]<-mean(expit((gibbs_out_pocrm1.I$alphaT)+wmT[ordT1.I,dose.level]*(gibbs_out_pocrm1.I$betaT))<toxbound)
        
        admiss[dose.level]<-  probtox[dose.level]>C_tox
      }
      admiss.doses<-which(admiss)
      
      #dose skipping
      if((dose.skipping.rule=="ON")|(dose.skipping.rule=="ON.diag.allowed")){
        allowable<-allowable_func(allowable.list.in = allowable.escalations,doses.so.far = dosevec)
        
        admiss.doses.only<-admiss.doses
        
        admiss.doses<-admiss.doses[admiss.doses%in%allowable]
        
      }
      
      #browser()
      utility1.I<-p.tox.est<-c()
      for(dose.level in admiss.doses){
        p.tox.est[dose.level]<-expit(mean(gibbs_out_pocrm1.I$alphaT)+wmT[ordT1.I,dose.level]*mean(gibbs_out_pocrm1.I$betaT))
        
        utility1.I[dose.level]<- -abs(p.tox.est[dose.level] - target)
        
      }
      
      
      ##choosing the next dose:
      
      nextdose<-which.max(utility1.I)
      
      utility<-utility1.I
      # if(ordT1.J!=ordT2){
      #   browser()
      # }
      
      ##if no admissable doses
      if(length(nextdose)==0){
        stop<-1
        stop_vec[1]<-1
        dose_rec<-NA
        break
      }    
      
   
      
      npats_doses<-sapply(c(1:ndoses), function(x) sum(current.patient.data$dose.level==x))
      #stopping rules
      if(((safety.stopping.low.unsafe==T)&(npats_doses[1]>0))|((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0))){
        
        
        current.patient.data_cyc1<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=1)
        
        
        
        weightT1<-c()
        patdoses1<-current.patient.data_cyc1$dose.level
        yDLT1<-current.patient.data_cyc1$current.DLT
        
        
        #
        if(TITE==T){
          
          for(obser in 1:nrow(current.patient.data_cyc1)){
            weightT1[obser]<-current.patient.data_cyc1$patient_weights_tox[obser]
            
          }
        }else{
          
          weightT1<-rep(1,nrow(current.patient.data_cyc1))
          
        }
        
        
        
        
        
        
        #only cycle 1 
        currentdata_pocrm_1cyc<-list(ndoses=ncol(wmT),yDLT=yDLT1,
                                     dosesT=wmT[ordT1.I,],
                                     patdoses=patdoses1,
                                     mu_alphaT=prior_vec2[1],mu_betaT=prior_vec2[3],
                                     tau_alphaT=prior_vec2[2],tau_betaT=prior_vec2[4],
                                     weightT=weightT1)
        
        
        
        
        gibbs_out1<-gibbs_sampler.crm.combo(data_list=currentdata_pocrm_1cyc,iter=5000,model_string=model.crm.joint.string)
        
        if((safety.stopping.low.unsafe==T)&(npats_doses[1]>0)){
          posterior1<- expit((gibbs_out1$alphaT)+wmT[ordT1.I,1]*(gibbs_out1$betaT))
          
          posterior1<-posterior1[1:length(posterior1)]
          cyc1_0.3g<-mean(posterior1>0.3,na.rm=T)
          
          if(cyc1_0.3g>0.8){
            stop<-6
            stop_vec[6]<-1
            nextdose<-NA
            
          }
        }
        
        if((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0)){
          posteriorJ<-expit((gibbs_out1$alphaT)+wmT[ordT1.I,ncol(wmT)]*(gibbs_out1$betaT))
          posteriorJ<-posteriorJ[1:length(posteriorJ)]
          cycJ_0.3l<-mean(posteriorJ<0.3,na.rm=T)
          if(cycJ_0.3l>0.8){
            
            stop<-7
            stop_vec[7]<-1
            nextdose<-NA
            
          }
        }
        
      }
      
      
      
      
      #number of DLTs per dose level
      nDLTs_doses<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)]))
      nDLTs_doses_c1<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)&(current.data$cycle==1)]))
      
      
      ##hard safety
      if(hard.safety==T){
        
        explored<-c(1:ndoses)[npats_doses>0]
        
        for(do in explored){
          
          if(nDLTs_doses_c1[do]>=hard.safety.mat[1, which(hard.safety.mat[2,]==npats_doses[do])]){
            
            min.over.indices.hold<-dose.indices[do,]
            excluded.hold<-as.numeric((dose.indices[,1]>=min.over.indices.hold[1])&(dose.indices[,2]>=min.over.indices.hold[2]))
            excluded<-excluded+excluded.hold
            
            
          }
        }
        if(!is.na(nextdose)){
          if(excluded[nextdose]>0){#if the next dose is in the exclusion zone
            included.utility<-utility
            included.utility[excluded>0]<--100
            included.utility[-admiss.doses]<--100
            nextdose<-which.max(included.utility)
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
          
          stop<-4
          stop_vec[4]<-1
          dose_rec<-nextdose
          
        }
      }
      
      
      
    }
    if(nrow(current.patient.data)==(ncohorts*co_size)){ #max patients reached
      
      
      current.patient.data<-current_patient_data_frame(current_time=(current.time+ncycles),patient.dataframe=all.data,follow_up = ncycles)
      patdoses<-current.patient.data$dose.level
      
      patdoses<-current.patient.data$dose.level
      
      if(TITE==T){
        
        weightT<-current.patient.data$patient_weights_tox
      }else{
        weightT<-rep(1,length(current.patient.data$patient_weights_tox))
      }
      #formatting for gibbs sampler
      #yOUT 
      #1= no E no T
      #2= no E yes T
      #3= yes E no T
      #4= yes E yes T
      
      y_T<-current.patient.data$current.DLT
      ########################################  
      
      ######### 1 par ind
      
      margT1<-c()
      
      
      
      
      
      
      
      for (kT in 1:length(prior.o.T)) { 
        # browser()
        
        margT1[kT]<-integrate(V.post.fun.ind.1par,lower=-20,upper=20,y=y_T,patdoses=patdoses,wm_skeleton=wmT[kT,],
                              weight=weightT,prior_hyp=prior_vec1[1:2])$value
        
        
        #print(kT)
      } # for kT
      
      #browser()
      #posterior
      
      pordT1 <- (margT1 * prior.o.T)/sum(margT1 * prior.o.T)
      #maximise the posterior - which ordering?
      
      ordT1.I<-which.is.max.func(pordT1)
      
      
      
      
      
      
      #calculate posterior for that ordering
      currentdata_pocrm1.I<-list(ndoses=ncol(wmT),weightT=weightT,yDLT=y_T,
                                 dosesT=wmT[ordT1.I,],patdoses=patdoses,
                                 mu_alphaT=prior_vec2[1],mu_betaT=prior_vec2[3],
                                 tau_alphaT=prior_vec2[2],tau_betaT=prior_vec2[4])
      
      #browser()
      
      
      gibbs_out_pocrm1.I<-gibbs_sampler.crm.combo(data_list=currentdata_pocrm1.I,iter=10000,model_string=model.crm.joint.string)
      #browser()
      #choose next dose
      
      
      #R value for the best ordering
      rpredT1.I <- expit(mean(gibbs_out_pocrm1.I$alphaT)+wmT[ordT1.I,]*mean(gibbs_out_pocrm1.I$betaT))
      
      #browser()
      
      # choose the next dose level
      
      #admissable doses
      admiss<-probeff<-probtox<-c()
      for(dose.level in 1:ndoses){
        probtox[dose.level]<-mean(expit((gibbs_out_pocrm1.I$alphaT)+wmT[ordT1.I,dose.level]*(gibbs_out_pocrm1.I$betaT))<toxbound)
        
        admiss[dose.level]<-  probtox[dose.level]>C_tox
      }
      admiss.doses<-which(admiss)
      #browser()
      utility1.I<-p.tox.est<-p.eff.est<-c()
      for(dose.level in admiss.doses){
        p.tox.est[dose.level]<-expit(mean(gibbs_out_pocrm1.I$alphaT)+wmT[ordT1.I,dose.level]*mean(gibbs_out_pocrm1.I$betaT))
        
        utility1.I[dose.level]<- -abs(p.tox.est[dose.level] - target)
        
      }
      
      
      ##choosing the next dose:
      
      nextdose1.I<-which.max(utility1.I)
      
      
      
      nextdose<-nextdose1.I
      utility<-utility1.I
      ##choosing the next dose:
      
      dose_rec<-which.max(utility)
      stop<-3
      stop_vec[3]<-1
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

