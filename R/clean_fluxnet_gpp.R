#modified based on functions from Beni's package-->ingestr:
#https://github.com/stineb/ingestr/blob/master/R/get_obs_bysite_fluxnet.R
#
clean_fluxnet_gpp <- function(df, remove_neg = FALSE, filter_ntdt){
  ##--------------------------------------------------------------------
  ## Cleans daily data using criteria 1-4 as documented in Tramontana et al., 2016
  ## gpp_nt: based on nighttime flux decomposition ("NT")
  ## gpp_dt: based on daytime flux decomposition ("DT")
  ##--------------------------------------------------------------------
  
  replace_with_na_qc <- function(gpp, qc, threshold){
    gpp[which(qc < threshold)] <- NA
    return(gpp)
  }
  replace_with_na_neg <- function(gpp){
    gpp[which(gpp<0)] <- NA
    return(gpp)
  }
  replace_with_na_res <- function(gpp, res, q025, q975){
    gpp[ res > q975 | res < q025  ] <- NA
    return(gpp)
  }
  #commented by YP: do not use this:
  # df <- df %>%
  #   mutate(GPP_NT_VUT_REF = replace_with_na_qc(GPP_NT_VUT_REF, NEE_VUT_REF_NIGHT_QC, threshold),
  #          GPP_DT_VUT_REF = replace_with_na_qc(GPP_DT_VUT_REF, NEE_VUT_REF_DAY_QC,   threshold))
  
  # ## Remove data points that are based on too much gap-filled data in the underlying half-hourly data
  # gpp_nt[ which(qflag_nt < threshold) ] <- NA  ## based on fraction of data based on gap-filled half-hourly
  # gpp_dt[ which(qflag_dt < threshold) ] <- NA  ## based on fraction of data based on gap-filled half-hourly
  
  if (filter_ntdt){
    ## Remove data points where the two flux decompositions are inconsistent,
    ## i.e. where the residual of their regression is above the 97.5% or below the 2.5% quantile.
    df <- df %>%
      mutate(res = GPP_NT_VUT_REF - GPP_DT_VUT_REF)
    
    q025 <- quantile( df$res, probs = 0.025, na.rm=TRUE )
    q975 <- quantile( df$res, probs = 0.975, na.rm=TRUE )
    
    
    ## remove data outside the quartiles of the residuals between the DT and NT estimates
    df <- df %>%
      mutate(GPP_NT_VUT_REF = replace_with_na_res(GPP_NT_VUT_REF, res, q025, q975),
             GPP_DT_VUT_REF = replace_with_na_res(GPP_DT_VUT_REF, res, q025, q975)
      )
  }
  ## remove outliers
  df <- df %>%
    mutate(GPP_NT_VUT_REF = remove_outliers(GPP_NT_VUT_REF, coef = 1.5),
           GPP_DT_VUT_REF = remove_outliers(GPP_DT_VUT_REF, coef = 1.5)
    )
  
  ## remove negative GPP
  if (remove_neg){
    df <- df %>%
      mutate(GPP_NT_VUT_REF = replace_with_na_neg(GPP_NT_VUT_REF),
             GPP_DT_VUT_REF = replace_with_na_neg(GPP_DT_VUT_REF)
      )
  }
  return(df)
}
