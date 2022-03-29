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

align_events <- function( df, dovars, leng_threshold, before, after, nbins, do_norm=FALSE ){
  #!comment by YPL:test using some df
  # df<-df_norm_RU_Fyo
  # dovars<-c("gpp")
  # df.sep<-df.sep20_RU_Fyo
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
  
  
  ## Bins for different variables XXX a bit weird with default values
  bins  <- seq( from=-before, to=after, by=(after+before)/nbins )
  
  ## merge df_isevent into df
  df <- df %>% mutate( idx_df = 1:nrow(df) )
  
  ##--------------------------------------------------------
  ## Identify events ()
  ##--------------------------------------------------------
  events <- get_consecutive( 
    df$is_event, 
    leng_threshold = leng_threshold, 
    do_merge       = FALSE
  )
  
  ##--------------------------------------------------------
  ## Re-arrange data, aligning by beginning of events
  ## Creates data frame where not all rows are retained from df
  ## and columns added for 'dday' (number of day relative to onset of event)
  ## and 'iinst' number of event to which row belongs.
  ##--------------------------------------------------------
  if (nrow(events)>1){
    
    df_dday <- c()
    pos_event<-c()   ##added by YP-->the position of "events"
    
    for ( iinst in 1:nrow(events) ){
      #after_inst <- min( after, events$len[iinst] )
      after_inst<-events$len[iinst]  #revised by YP:selected the data until end of events
      dday <- seq( from=-before, to=after_inst, by=1 )
      idxs <- dday + events$idx_start[iinst]
      pos_event<-c(pos_event,idxs)
      drophead <- which( idxs < 1 )  #comment by YP:keep the data when the idxs>1
      if (length(drophead)>0){
        idxs <- idxs[ -drophead ]
        dday <- dday[ -drophead ]
      }
      addrows <- df %>% slice( idxs ) %>% mutate( dday=dday, inst=iinst )
      df_dday <- df_dday %>% bind_rows(addrows)   #added by YP: still not figure out why this is not right   
      # df_dday <- rbind(df_dday,addrows) 
      print(i)
    }
    #added by YP:tidy the non-event data
    pos_nonevent<-setdiff(df$idx_df,pos_event)
    df_noevent_dday<-df[pos_nonevent,]
    
    ##--------------------------------------------------------
    ## Normalise re-arranged data relative to a certain bin's median for each site
    ##--------------------------------------------------------
    #commented by YP:now do not use this part-->so everytime set do_norm==FALSE
    if (do_norm){
      ## add column for bin
      df_dday <- df_dday %>% mutate( inbin  = cut( as.numeric(dday), breaks = bins ) )
      
      ## Normalise by median value in dday-bin before drought onset ("zero-bin")
      ## Get median in zero-bin
      sdovars  <- paste0("s",  dovars)
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
      
      norm <- df_dday %>% group_by( site, inbin ) %>%
        summarise_at( vars(one_of(sdovars)), funs(median( ., na.rm=TRUE )) ) %>%
        filter( !is.na(inbin) ) %>% 
        filter( grepl(",0]", inbin) ) %>% 
        setNames( c( "site", "inbin", paste0("d", sdovars) ) ) %>% 
        select(-inbin)
      
      df_dday <- df_dday %>% 
        left_join(norm, by="site") %>% 
        ungroup()
      
      ## Divide by median in zero-bin
      get_dsdovar <- function(df, sdovar){
        dsdovar <- paste0("d", sdovar)
        df[[dsdovar]] <- df[[sdovar]] / df[[dsdovar]]
        return(select(df, dsdovar))
      }
      df_dday <- purrr::map_dfc(as.list(sdovars), ~get_dsdovar(df_dday, .)) %>% 
        bind_cols( select(df_dday, -one_of(dsdovars)), .)
      
    } else {
      sdovars <- c()
      dsdovars <- c()
    }
    
    ##--------------------------------------------------------
    ## Aggregate accross events, by site
    ##--------------------------------------------------------
    df_dday_agg_inst <- df_dday %>%  
      group_by( sitename, dday ) %>%    #comment by YP:here change site-->sitename
      summarise_at( 
        vars(one_of(dovars, sdovars, dsdovars)), 
        list( ~median( ., na.rm=TRUE), ~q33( ., na.rm=TRUE), ~q66( ., na.rm=TRUE) ) )
    
    ##--------------------------------------------------------
    ## Aggregate accross events and sites
    ##--------------------------------------------------------
    df_dday_agg_inst_site <- df_dday %>%  
      group_by( dday ) %>% 
      summarise_at( 
        vars(one_of(dovars, sdovars, dsdovars)), 
        list( ~median( ., na.rm=TRUE), ~q33( ., na.rm=TRUE), ~q66( ., na.rm=TRUE) ) )
    
  } else {
    
    df_dday_agg_inst      <- NULL
    df_dday_agg_inst_site <- NULL
    
  }
  
  out <- list( 
    df_dday = df_dday, 
    df_noevent_dday=df_noevent_dday,
    df_dday_agg_inst = df_dday_agg_inst, 
    df_dday_agg_inst_site = df_dday_agg_inst_site,
    bins = bins,
    norm = norm
  )
  return( out )
  
}
