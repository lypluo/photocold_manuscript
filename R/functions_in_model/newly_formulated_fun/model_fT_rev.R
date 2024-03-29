#################################################
##considering the temperature modifier f_Ts(s stands for State of acclimation)
##refer Mekela et al., 2004 and Mekela et al., 2008
##also add a scaler to improve the GPP in the mid-season
#################################################
#
f_Ts_rev<-function(df,
     par=c("tau"=5,"X0"=-10,"Smax"=5,"k"=1),
     plot=FALSE){

  ##parameters:
  #tau:represents the speed of response of the current acclimation status to
  #changes in temperature (Tian et al., 2021)
  #X0:lower temperature limit parameter(degreeC) above which f_Ts>0 and on the
     #state of acclimation
  #Smax:minimum temperature threshold parameter at which canopy photosynthesis is
     #not limited by low temperature

  ## data frame df must contain columns:
  ## 'tmin': daily minimum temperature
  ## "doy"': day of the year
  #initial setup:
  f_stress<- rep(NA, nrow(df)) #f_stress:modifying factor by Temperature
  # S<-rep(NA,nrow(df)) #S:state of acclimation
  X<-rep(NA,nrow(df)) #X:the intermediate variable that used to calculate S


  #
  for (idx in seq(nrow(df))){
    #---------
    #1):X-->calcuating the intermediate values:
    #---------
    if(idx==1){
      X[idx]=df$tmin[idx]
    }
    #
    if(idx>1){
      X[idx]=X[idx-1]+c(df$tmin[idx]-X[idx-1])/par["tau"]
    }

    #---------
    #2):S-->state of acclimation
    #---------
    # S[idx]<-max(0, (X[idx] - par["X0"]))
    S<-max(0, (X[idx] - par["X0"]))
    #-----------
    #3):f_stress:modifying factor for low temperature
    #-----------
    # f_stress[idx]<-min(S[idx]/par["Smax"],1)
    f_stress_temp<-min(S/par["Smax"],1)
    #---------
    #(4):stress factor multiplied by a scalar to
    ## allow for higher mid-season GPP after down-scaling early season GPP
    ## set different scaler at different period
    #between 4-01(from April air temperature starts >0) and 9-22(autumn equinox)
    #---------
    if(df$doy[idx]>121 & df$doy[idx]<266){
      f_stress[idx] <- f_stress_temp*par["k"]
    }else
    { f_stress[idx] <- f_stress_temp}
  }

  #plotting for evaluation:
  if(plot){
    plot(f_stress)
  }
  #
  return(f_stress)
}


#update: add additional parameter m in 2022-07:
f_Ts_rev_2pars<-function(df,
                   par=c("tau"=5,"X0"=-10,"Smax"=5,"k"=1,"m"=1),
                   plot=FALSE){
  ## data frame df must contain columns:
  ## 'tmin': daily minimum temperature
  ## "doy"': day of the year
  #initial setup:
  f_stress<- rep(NA, nrow(df)) #f_stress:modifying factor by Temperature
  # S<-rep(NA,nrow(df)) #S:state of acclimation
  X<-rep(NA,nrow(df)) #X:the intermediate variable that used to calculate S
  
  
  #
  for (idx in seq(nrow(df))){
    #---------
    #1):X-->calcuating the intermediate values:
    #---------
    if(idx==1){
      X[idx]=df$tmin[idx]
    }
    #
    if(idx>1){
      X[idx]=X[idx-1]+c(df$tmin[idx]-X[idx-1])/par["tau"]
    }
    
    #---------
    #2):S-->state of acclimation
    #---------
    # S[idx]<-max(0, (X[idx] - par["X0"]))
    S<-max(0, (X[idx] - par["X0"]))
    #-----------
    #3):f_stress:modifying factor for low temperature
    #-----------
    # f_stress[idx]<-min(S[idx]/par["Smax"],1)
    f_stress_temp<-min(S/par["Smax"],1)
    #---------
    #(4):stress factor multiplied by a scalar to
    ## allow for higher mid-season GPP after down-scaling early season GPP
    ## set different scaler at different period
    #between 4-01(from April air temperature starts >0) and 9-22(autumn equinox)
    #:k1-->no T limitation
    #other period:k0
    #---------
    if(df$doy[idx]>121 & df$doy[idx]<266){
      f_stress[idx] <- f_stress_temp*par["k"]
    }else
    { f_stress[idx] <- f_stress_temp*par["m"]}
  }
  
  #plotting for evaluation:
  if(plot){
    plot(f_stress)
  }
  #
  return(f_stress)
}
