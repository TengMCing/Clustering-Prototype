Define_Interval <- function(time_id, t, ActiveTime){
  left <- max(1, t - ActiveTime)
  right <- t
  indexes <- (time_id >= left) & (time_id <= right)
  if (sum(indexes) == 0) return(NULL)
  which(indexes)
}

Dist_Point_To_Vector <- function(plon, plat, vlon, vlat){
  
  dist_mat <- suppressMessages(geodist::geodist_vec(plon, plat, vlon, vlat))
  as.vector(dist_mat)
}

Nearby_Hotspots <- function(hotspot_list, pointer, lon, lat, AdjDist){

  if (length(lon) == 1) return(NULL)
  
  dist_vector <- Dist_Point_To_Vector(lon[pointer], lat[pointer], lon, lat)
  potential <- which(dist_vector <= AdjDist)
  indexes <- !(potential %in% hotspot_list)
  
  if (sum(indexes) == 0) return(NULL)
  
  potential[!(potential %in% hotspot_list)]
}


Local_Clustering <- function(lon, lat, AdjDist){
  
  if (length(lon) == 1) return(c(1))
  
  hotspots_list <- c(1)
  pointer <- c(1)
  pointer_pos <- 1
  
  memberships <- NULL
  label <- NULL
  
  while (TRUE){
    
    while (TRUE) {
      
      nearby_hospots <- Nearby_Hotspots(hotspots_list, pointer, lon, lat, AdjDist)
      if (!is.null(nearby_hospots)) hotspots_list <- c(hotspots_list, nearby_hospots)
      
      if (pointer_pos < length(hotspots_list)) {
        pointer_pos <- pointer_pos + 1
        pointer <- hotspots_list[pointer_pos]
      } else {
        break
      }
      
    }
    
    if (is.null(memberships)) {
      memberships <- rep(1, length(hotspots_list))
      label <- 1
    } else {
      label <- label + 1
      new_len <- length(hotspots_list) - length(memberships)
      memberships <- c(memberships, rep(label, new_len))
    }
    
    
    indexes <- (!(1:length(lon) %in% hotspots_list))
    
    if (sum(indexes) == 0) break
    
    hotspots_list <- c(hotspots_list, min(which(indexes)))
    pointer_pos <- pointer_pos + 1
    pointer <- hotspots_list[pointer_pos]
    
  }
  
  memberships

}

Update_Memberships <- function(lon, lat, global_memberships, local_memberships, indexes){
  
  if (sum(global_memberships[indexes]) == 0) {
    global_memberships[indexes] <- local_memberships
    return(global_memberships)
  }
  
  
}

Global_Clustering <- function(lon, lat, time_id, ActiveTime, AdjDist){
  
  global_memberships <- rep(0, length(lon))
  
  start_time <- Sys.time()
  
  pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", 
                                   total = max(time_id))
  pb$tick(0)
  
  for (t in 1:max(time_id)){
    
    pb$tick(1)
    
    indexes <- Define_Interval(time_id, t, ActiveTime)
    if (is.null(indexes)) next
    
    local_memberships <- Local_Clustering(lon[indexes], lat[indexes], AdjDist)
    # global_memberships <-  Update_Memberships(global_memberships, local_memberships)
  }
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  print(time_taken)
  message(paste0("Time taken: ",  
                 as.numeric(time_taken, units = "secs") %/% 60, 
                 " mins ", 
                 round(as.numeric(time_taken, units = "secs") %% 60, 0),
                 " secs for ",
                 length(lon),
                 " data points"
                 )
          )
  
  message(paste0("(", 
                 round(as.numeric(time_taken, units = "secs")/length(lon), 3), 
                 " secs per data point)"))
  
  return(123)
  
  # global_memberships
}

########### TEST ###################################################################

hotspots <- readr::read_csv("data/VIC_hotspots_before_clustering.csv")


Local_Clustering(hotspots$lon[13:20], hotspots$lat[13:20], 3000)

Global_Clustering(hotspots$lon, 
                  hotspots$lat, 
                  hotspots$time_id, 
                  ActiveTime =  24,
                  AdjDist = 3000) -> temp

