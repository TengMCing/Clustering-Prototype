################# DEFINE INTERVAL (STEP 1) #####################################

Define_Interval <- function(time_id, t, ActiveTime){
  left <- max(1, t - ActiveTime)
  right <- t
  indexes <- (time_id >= left) & (time_id <= right)
  if (sum(indexes) == 0) return(NULL)
  which(indexes)
}

################################################################################
################# POINT TO VECTOR GEODISC ######################################

Dist_Point_To_Vector <- function(plon, plat, vlon, vlat){
  
  dist_mat <- suppressMessages(geodist::geodist_vec(plon, plat, vlon, vlat))
  as.vector(dist_mat)
}

################################################################################
################# FIND NEARBY HOTSPOTS FOR A HOTSPOT ###########################

Nearby_Hotspots <- function(hotspot_list, pointer, lon, lat, AdjDist){

  if (length(lon) == 1) return(NULL)
  
  dist_vector <- Dist_Point_To_Vector(lon[pointer], lat[pointer], lon, lat)
  potential <- which(dist_vector <= AdjDist)
  indexes <- !(potential %in% hotspot_list)
  
  if (sum(indexes) == 0) return(NULL)
  
  potential[!(potential %in% hotspot_list)]
}

################################################################################
################# LOCAL CLUSTERING (STEP 2) ####################################


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
  
  memberships[order(hotspots_list)]

}

################################################################################
################# ADJUST MEMBERSHIPS ###########################################

Adjust_Memberships <- function(memberships, max_membership){
  
  as.numeric(factor(memberships)) + max_membership
}

################################################################################
################# UPDATE MEMBERSHIPS (STEP 3) ##################################

Update_Memberships <- function(lon, lat, global_memberships, local_memberships, indexes){
  
  if (sum(global_memberships[indexes]) == 0) {
    global_memberships[indexes] <- Adjust_Memberships(local_memberships, max(global_memberships))
    return(global_memberships)
  }
  
  if (all(global_memberships[indexes] != 0)) {
    return(global_memberships)
  }
  
  fin_memberships <- global_memberships[indexes]
  local_lon <- lon[indexes]
  local_lat <- lat[indexes]
    
  new_p <- which(fin_memberships == 0)
  old_p <- which(fin_memberships != 0)
  
  shared_clusteres <- unique(local_memberships[old_p])
  
  type1 <- new_p[local_memberships[new_p] %in% shared_clusteres]
  type2 <- new_p[!local_memberships[new_p] %in% shared_clusteres]
  
  for (i in type1) {
    
    bool <- local_memberships[old_p] == local_memberships[i]
    
    current_old <- old_p[bool]
    
    dist_vector <- Dist_Point_To_Vector(local_lon[i],
                                        local_lat[i],
                                        local_lon[current_old],
                                        local_lat[current_old])
    
    
    target <- current_old[which.min(dist_vector)]
    
    fin_memberships[i] <- fin_memberships[target]
    
  }
  
  if (length(type2) != 0){
    fin_memberships[type2] <- Adjust_Memberships(local_memberships[type2], max(global_memberships))
  }
  
  global_memberships[indexes] <- fin_memberships
  return(global_memberships)

}

################################################################################
################# GLOBAL CLUSTERING (MAIN) #####################################

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
    if (all(global_memberships[indexes] != 0)) next
    
    local_memberships <- Local_Clustering(lon[indexes], lat[indexes], AdjDist)
    global_memberships <-  Update_Memberships(lon, 
                                              lat, 
                                              global_memberships, 
                                              local_memberships, 
                                              indexes)
  }
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
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
  
  return(global_memberships)
}

################################################################################
################# TEST #########################################################

hotspots <- readr::read_csv("data/VIC_hotspots_before_clustering.csv")


Local_Clustering(hotspots$lon[21:30], hotspots$lat[21:30], 3000)

Global_Clustering(hotspots$lon,
                  hotspots$lat,
                  hotspots$time_id,
                  ActiveTime = 24,
                  AdjDist = 3000) -> temp

max(temp)

# Note: the result is slightly different with the Python version due to the use of different approximation method of the geodesic.
# Time taken: 12 mins 21 secs > 6 mins (Python version)
