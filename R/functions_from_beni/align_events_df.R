#' Aligns data by events-->code mainly adopt from Beni
#'
#' Uses a vectory specifying whether data falls into an event to reshape data, aligning by the onset of the event
#' 
#' @param df A data frame containing all data continuously along time. The data frame must contain one column of type logical
#' named \code{"isevent"}, specifying whether respective dates satisfy a user-defined condition that is used to define events. 
#' Events are determined by the function based on consecutive dates, where this condition is satisfied (minimum length for 
#' defining an event is given by \code{leng_threshold}).
#' @param dovars A vector of character strings specifying which columns of \code{df} to re-arrange.
#' @param before An integer specifying the number of days before the event onset to be retained in re-arranged data 
#' @param after An integer specifying the number of days after the event onset to be retained in re-arranged data 
#' @param do_norm A logical specifying whether re-arranged data is to be normalised with respect to its value before the drought onset

#'
#' @return An aligned data frame
#' @export
#'
#' @examples df_alg <- align_events( df, before=30, after=300 )
#' 

# source the function:
# source.path<-'C:/Users/yluo/Desktop/R_testcode/PhotoCold/Second_round_of_code/R/Functions/functions_from_beni/'
# source(file=paste0(source.path,"get_consecutive.R"))
#
q33 <- function( vec, ... ){
  quantile( vec, 0.33, ...)
}

q66 <- function( vec, ... ){
  quantile( vec, 0.66, ...)
}

align_events <- function( df,dovars,leng_threshold, before, after, nbins, do_norm=FALSE ){
  #!comment by YPL:test using some df
  # df<-df_events.years
  # dovars<-c("gpp")
  # leng_threshold<-5
  # before=60
  # after=0
  # nbins=10
  # do_norm=FALSE
  
  require(plyr)
  require( dplyr )
  require( tidyr )
  
  #!comment by YPL:here changing the "isevent" to "is_event"
  if (!("is_event" %in% names(df))){   
    rlang::abort("align_events(): Column named isevent is missing in data frame df.")
  }
  

  ## merge df_isevent into df
  df <- df %>% mutate( idx_df = 1:nrow(df) )
  
  ##--------------------------------------------------------
  ## Identify events ()
  ##--------------------------------------------------------
  events <- get_consecutive(
    df$is_event,
    leng_threshold = leng_threshold,
    do_merge       = TRUE,
    merge_len = 10
  )
  
  ##--------------------------------------------------------
  ## Re-arrange data, aligning by beginning of events
  ## Creates data frame where not all rows are retained from df
  ## and columns added for 'dday' (number of day relative to onset of event)
  ## and 'iinst' number of event to which row belongs.
  ##--------------------------------------------------------
  if(nrow(events)==1){
    df_dday <- c()
    
    after_inst<-events$len 
    dday <- seq( from=-before, to=after_inst, by=1 )
    idxs <- dday + events$idx_start
    #
    # pos_event<-c(pos_event,idxs)
    drophead <- which( idxs < 1 )  #comment by YP:keep the data when the idxs>1
    if (length(drophead)>0){
      idxs <- idxs[ -drophead ]
      dday <- dday[ -drophead ]
    }
    addrows <- df %>% slice( idxs ) %>% mutate( dday=dday, inst=1 )
    df_dday <- rbind(df_dday,addrows)  
    
    ## Bins for different variables XXX a bit weird with default values
    #changed by YP for Bins
    if(after==0){
      bins  <- seq( from=-before, to=after, by=(after+before)/nbins )
    }
    if(after>0){
      bins<-seq( from=-before, to=max(df_dday$dday), by=c(max(df_dday$dday)+before)/nbins )
    }
    
  }
  if (nrow(events)>1){
    
    df_dday <- c()
    pos_event<-c()   ##added by YP-->the position of "events"
    
    for ( iinst in 1:nrow(events) ){
      #after_inst <- min( after, events$len[iinst] )
      after_inst<-events$len[iinst]  #revised by YP:selected the data until end of events
      dday <- seq( from=-before, to=after_inst, by=1 )
      idxs <- dday + events$idx_start[iinst]
      #
      pos_event<-c(pos_event,idxs)
      drophead <- which( idxs < 1 )  #comment by YP:keep the data when the idxs>1
      if (length(drophead)>0){
        idxs <- idxs[ -drophead ]
        dday <- dday[ -drophead ]
      }
      addrows <- df %>% slice( idxs ) %>% mutate( dday=dday, inst=iinst )
      df_dday <- rbind(df_dday,addrows)             
    }
    norm<-""
    
    ## Bins for different variables XXX a bit weird with default values
    #changed by YP for Bins
    if(after==0){
      bins  <- seq( from=-before, to=after, by=(after+before)/nbins )
    }
    if(after>0){
      bins<-seq( from=-before, to=max(df_dday$dday), by=c(max(df_dday$dday)+before)/nbins )
    }
    
    ##--------------------------------------------------------
    ## Normalise re-arranged data relative to a certain bin's median for each site
    ##--------------------------------------------------------
    #commented by YP:
    if (do_norm){
      ## add column for bin
      df_dday <- df_dday %>% mutate( inbin  = cut( as.numeric(dday), breaks = bins ) )
      
      ## Normalise by median value in dday-bin before drought onset ("zero-bin")
      ## Get median in zero-bin
      #added by YP-->do not use sdovars, only use dovars and dsdovars,add medium var called meddovars
      # sdovars  <- paste0("s",  dovars)
      meddovas<-paste0("med",dovars)
      dsdovars <- paste0("ds", dovars)
      
      # ## Add median in zero-bin (dsdovars), separate for each site, aggregated across instances (drought events)
      # df_dday <- df_dday %>% group_by( site, inbin ) %>%
      #   summarise_at( vars(one_of(sdovars)), funs(median( ., na.rm=TRUE )) ) %>%
      #   filter( !is.na(inbin) ) %>% 
      #   filter( grepl(",0]", inbin) ) %>% 
      #   setNames( c( "site", "inbin", paste0("d", sdovars) ) ) %>% 
      #   select(-inbin) %>% 
      #   right_join(df_dday, by="site") %>% 
      #   ungroup()
      
      #changed by YP
      #-->change "site" to "sitename" below in order to process the data
      norm <- df_dday %>% 
        group_by( sitename, inbin ) %>%
        summarise_at( vars(one_of(dovars)), funs(median( ., na.rm=TRUE )) ) %>%
        dplyr::filter( !is.na(inbin) ) %>% 
        dplyr::filter( grepl(",0]", inbin) ) %>% 
        setNames( c( "sitename", "inbin", paste0("med", dovars) ) ) %>% 
        dplyr::select(-inbin)
      
      df_dday <- df_dday %>% 
        left_join(norm, by="sitename") %>% 
        ungroup()
      
      ## Divide by median in zero-bin([..,0])
      # added by YP:make the median value as the center-->i.e. for istance for c(1,2,3,4,5)-->
      # normalized by (1-3)/3,(2-3)/3....
      # -->hence,changing the beni's function below
      get_dsdovar <- function(df, dovar){
        # df<-df_dday
        # dovar<-"gpp_obs"

        meddovar<-paste0("med",dovar)
        dsdovar <- paste0("ds", dovar)
        # df[[dsdovar]] <- df[[dovar]] / df[[dsdovar]]
        #changed by YP-->normalize the data as making the "dsdovar" as the standard
        df[[dsdovar]] <- c(df[[dovar]] - df[[meddovar]])/abs(df[[meddovar]])
        return(dplyr::select(df, dsdovar))
      }
      df_dday <- purrr::map_dfc(as.list(dovars), ~get_dsdovar(df_dday, .)) %>% 
        bind_cols(dplyr::select(df_dday, -one_of(dsdovars)), .)
      
    } else {
      sdovars <- c()
      dsdovars <- c()
    }
    
  }
  
  #added by YP-->test if the normalization makes sense-->it makes sense...
  # tt_site<-df_dday[df_dday$sitename=="CA-Man",]
  # tt_site_seb<-tt_site[,c("sitename","date","Year","dday","inbin","gpp_obs","medgpp_obs","dsgpp_obs",
  #                        "temp_day_fluxnet2015","medtemp_day_fluxnet2015","dstemp_day_fluxnet2015")]
  #test plot:
  # library(ggplot2)
  # tt_site_seb$Year<-as.factor(tt_site_seb$Year)
  # ggplot(data = tt_site_seb,aes(x=dday,y=dsgpp_obs,col=Year))+
  #   geom_point()+
  #   geom_line()+
  #   geom_hline(yintercept = 0,lty=2,size=1.1)+
  #   geom_vline(xintercept = 0,lty=2,size=1.1)+
  #   theme_classic()
    
  ######
  out <- list( 
    df_dday = df_dday, 
    # df_noevent_dday=df_noevent_dday,
    bins = bins,
    norm = norm
  )
  return( out )
}
