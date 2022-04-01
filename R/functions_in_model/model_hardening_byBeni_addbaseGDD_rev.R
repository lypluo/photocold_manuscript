#################################################
##Based on Beni's function
#################################################
#function for the hardening
f_hardening <- function(temp, par){
  #changing the function according to the logistic function in weikipeidia:
  #https://en.wikipedia.org/wiki/Logistic_function
  xx <- temp # * ppfd
  xx <- (-1)*par["b"] * (xx - par["a"])
  yy <- 1 / (1 + exp(xx))
  return(yy)
}

#test plot
# ggplot() +
#   geom_function(fun = f_hardening, args = list(par = c("a" = 0, "b" = 20))) +
#   xlim(-10, 30)

#function for the dehardening
f_dehardening <- function(temp, par){
  #changing the function according to the logistic function in weikipeidia:
  #https://en.wikipedia.org/wiki/Logistic_function
  xx <- temp # * ppfd
  xx <- (-1) *par["d"] * (xx - par["c"])
  yy <- 1 / (1 + exp(xx))
  return(yy)
}
#test plot
# ggplot() +
#   geom_function(fun = f_dehardening, args = list(par = c("c" = 50, "d" = 0.1))) +
#   xlim(-10, 100)

#one parameters(e) for scaling
model_hardening_1par <- function(df,
                                 par = c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1,"k" = 5),
                                 plot = FALSE){

  ## data frame df must contain columns:
  ## 'temp': daily mean temperature
  ## 'tmin': daily minimum temperature

  level_hard <- 1.0  # start without hardening
  gdd <- 0 #growing degree day
  f_stress <- rep(NA, nrow(df))

  for (idx in seq(nrow(df))){

    ## determine hardening level - responds instantaneously to minimum temperature
    level_hard_new <-  f_hardening(df$tmin[idx], par)

    if (level_hard_new < level_hard){

      ## entering deeper hardening
      level_hard <- level_hard_new

      ## re-start recovery
      gdd <- 0

      # print(paste("Hardening to", level_hard, "on", df$date[idx]))
    }

    ## accumulate growing degree days (GDD)
    gdd <- gdd + max(0, (df$temp[idx] - par["k"])) #set parameter k for base GDD

    ## de-harden based on GDD. f_stress = 1: no stress
    level_hard <- level_hard + (1-level_hard) * f_dehardening(gdd, par)

    ## stress function is hardening level multiplied by a scalar to
    ## allow for higher mid-season GPP after down-scaling early season GPP
    ## set different scaler at different period
    #between 3-22(spring equinox) and 9-22(autumn equinox)
    if(df$doy[idx]>81 & df$doy[idx]<266){
      f_stress[idx] <- level_hard*par["e"]
    }else
    { f_stress[idx] <- level_hard}

  }
  if(plot){
    plot(f_stress)
  }

  return(f_stress)
}

###two parameters(e,f) for scaling
model_hardening_2par <- function(df,
    par = c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1,"f"=1,"k" = 5),
    plot = FALSE){

  ## data frame df must contain columns:
  ## 'temp': daily mean temperature
  ## 'tmin': daily minimum temperature

  level_hard <- 1.0  # start without hardening
  gdd <- 0 #growing degree day
  f_stress <- rep(NA, nrow(df))

  for (idx in seq(nrow(df))){

    ## determine hardening level - responds instantaneously to minimum temperature
    level_hard_new <-  f_hardening(df$tmin[idx], par)

    if (level_hard_new < level_hard){

      ## entering deeper hardening
      level_hard <- level_hard_new

      ## re-start recovery
      gdd <- 0

      # print(paste("Hardening to", level_hard, "on", df$date[idx]))
    }

    ## accumulate growing degree days (GDD)
    gdd <- gdd + max(0, (df$temp[idx] - par["k"])) #set parameter k for base GDD

    ## de-harden based on GDD. f_stress = 1: no stress
    level_hard <- level_hard + (1-level_hard) * f_dehardening(gdd, par)

    ## stress function is hardening level multiplied by a scalar to
    ## allow for higher mid-season GPP after down-scaling early season GPP
    ## set different scaler at different period
    #between 3-22(spring equinox) and 9-22(autumn equinox)
    if(df$doy[idx]>81 & df$doy[idx]<266){
      f_stress[idx] <- level_hard*par["e"]
    }else
    { f_stress[idx] <- level_hard*par["f"]}

  }
  if(plot){
    plot(f_stress)
  }

  return(f_stress)
}

