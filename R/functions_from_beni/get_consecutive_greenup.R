#---------------------
#code adopted mainly from Beni for the drought-->most focusing on the "GPP overestimation" in my project:
#---------------------
get_consecutive_greenup <- function( dry, leng_threshold=5, anom=NULL, do_merge=FALSE,merge_len){
  ##------------------------------------
  ## Returns a dataframe that contains information about events (starting index and length) 
  ## of consecutive conditions (TRUE) in a boolean vector ('dry' - naming is a legacy).
  ##------------------------------------
  #!comment by YPL:test using some df:
  # dry<-df$greenup
  # leng_threshold=5
  # do_merge=TRUE
  # merge_len=10
    
  ## replace NAs with FALSE (no drought). This is needed because of NAs at head or tail
  dry[ which(is.na(dry)) ] <- FALSE 
  #!comment by YPL: change the value of "yes" and "no" to 1, and NA as I used "no" and "yes" in my data.frame
  dry[dry=="yes"]<-1
  dry[dry=="no"]<-NA
  dry<-as.vector(as.numeric(dry))
  
  #
  ## identifies periods where 'dry' true for consecutive days of length>leng_threshold and 
  ## creates data frame holding each instance's info: start of drought by index in 'dry' and length (number of days thereafter)
  ##!!have a test#
  #using 1:365
  # dry<-dry[1:365]
  
  instances <- data.frame( idx_start=c(), len=c() )
  consecutive_dry <- rep( NA, length( dry ) )
  ndry  <- 0
  ninst <- 0
  for ( idx in 1:length( dry ) ){
    if (!is.na(dry[idx])){ 
      ndry <- ndry + 1 
    } else {
      if (ndry>=leng_threshold) { 
        ## create instance
        ninst <- ninst + 1
        addrow <- data.frame( idx_start=idx-(ndry), len=ndry )
        instances <- rbind( instances, addrow )
      }
      ndry <- 0
    }
    consecutive_dry[idx] <- ndry
  }
  if (ndry>leng_threshold){
    ## create a last instance if the last dry period extends to the end of the time series
    ninst <- ninst + 1
    addrow <- data.frame( idx_start=idx-(ndry), len=ndry )
    instances <- rbind( instances, addrow )
  }
  
  
  if (nrow(instances)>0){

  #   ## Get cumulative deficit per instance (deficit w.r.t. 1, where 'anom' is a vector with values 0-1)
  #   if (!is.null(anom)){
  #     instances$deficit <- rep( NA, nrow(instances) )
  #     for ( idx in 1:nrow(instances) ){
  #       instances$deficit[idx] <- sum( anom[ instances$idx_start[idx]:(instances$idx_start[idx]+instances$len[idx]-1) ] )
  #     }
  #   }
    
    ## merge droughts interrupted by short non-drought periods
    ## if in-between non-drought period is shorter than both of the drought periods
    ## before and after non-drought period
    if (do_merge){
      
      print("dimensions of instances before merging short periods")
      print(dim(instances))
      
      # ninst_save <- nrow( instances ) + 1
      # ninst      <- nrow( instances )
      
      # while (ninst < ninst_save){
        
      # ninst_save <- nrow( instances )
        
        instances_merged <- data.frame( idx_start=c(), len=c() )
        
        idx <- 0
        while (idx<(nrow(instances)-1)){
          idx <- idx + 1
          
          len_betweendrought <- instances$idx_start[idx+1] - (instances$idx_start[idx] + instances$len[idx] + 1)
          
          # if (len_betweendrought<instances$len[idx] && len_betweendrought<instances$len[idx+1]){
          #adjusted by YP-->merge crition: if the len_betweendrought - len[idx]<=10-->small gaps-->then merge them
          if(c(len_betweendrought<=merge_len)){
            addrow <- data.frame( idx_start=instances$idx_start[idx], len=(instances$idx_start[idx+1] + instances$len[idx+1]) - instances$idx_start[idx] )
            instances_merged <- rbind( instances_merged, addrow )
            idx <- idx + 1
          } else {
            instances_merged <- rbind( instances_merged, instances[idx,] )
            if (idx==(nrow(instances)-1)){
              instances_merged <- rbind( instances_merged, instances[idx+1,] )
            }
          }
        }
        
        instances <- instances_merged
        rownames(instances)<-c(1:nrow(instances))
        ninst <- nrow( instances )
        
        print( "dimensions of instances after merging short periods" )
        print( dim( instances ) )
        
      }
      
    }
    
  # }
  
  return( instances )
}