library(raster)
library(gstat)
library(automap)
library(parallel)
library(doParallel)


setwd("/media/project/abertoncini/01_Precip_Network/08_Chinook_Runs_20230724")

#setwd("C:/Users/alb818/Dropbox/PHD/PRECIP_NETWORK/03_Precip_Network_20230621/00_Lapse_Check")


#Load precipitation data

precip_input <- read.csv("Precip_Undercatch_Corrected_20210727.csv")[,-c(1)]

precip_input <- ifelse(as.matrix(precip_input) > 160, NA, as.matrix(precip_input)) #based on climate normals

#Load station information

station_info <- read.csv("Station_Info_Inside_Domain_20210727.csv")[,-c(1)]


#Load quadrant information

quad_coord <- read.csv("Quadrants_Coordinates_20220416.csv", sep = ";")[,-c(1)]


#Load DEM

srtm <- raster("SRTM_90m_VoidFilled_GEE_v1_Mosaic_Clip.tif")

#srtm <- raster("SRTM_90m_Tile_6.tif")

srtm[srtm == -32768] <- NA
srtm[srtm == 32767] <- NA


#Parallelization by spatial quadrant

qs <- c(1, 24, 47, 70, 93)

qe <- c(23, 46, 69, 92, 115)


for (q_loop in 1:5) { 

no_cores <- (qe[q_loop] - qs[q_loop]) + 1

cl <- makeCluster(no_cores, type = "PSOCK")

registerDoParallel(cl)


foreach(q = qs[q_loop]:qe[q_loop]) %dopar% {
  
  
  library(raster)
  library(gstat)
  library(automap)
  library(parallel)
  library(doParallel)

  
  #Crop DEM for faster processing
  
  crop_extent <- extent(c(quad_coord$long_min[q], quad_coord$long_max[q], quad_coord$lat_min[q], quad_coord$lat_max[q]))
  
  #crop_extent <- extent(c(-115.3, -115.1, 50.85, 50.95))
  
  srtm <- crop(srtm, crop_extent)
  
  #plot(srtm)
  
  
  #Create data.frame with DEM xyz
  
  srtm_points <- as.data.frame(rasterToPoints(srtm))
  
  names(srtm_points) <- c("x", "y", "z")
  
  coordinates(srtm_points) <- ~ x + y
  
  
  #Select day to calculate variogram
  
  #Define available days
  
  days <- seq(as.POSIXct("1990-10-01", format = "%Y-%m-%d", tz = "GMT"),
              as.POSIXct("2020-09-30", format = "%Y-%m-%d", tz = "GMT"), by = 86400)
  
  
  #Define lapse rate and uncertainty
  
  lapse_rate_vector <- vector()
  
  lapse_uncertainty_vector <- vector()
  
  for (o in (which(substr(days, 1, 10) == "1990-10-01")):(which(substr(days, 1, 10) == "2020-09-30"))) {
    
    
    #Retrieve station information
    
    lat <- station_info$lat
    
    long <- station_info$long
    
    elev <- station_info$elev
    
    
    #Calculate the number of NAs for all stations for this particular date
    
    precip_date <- precip_input[o,]
    
    na_vector <- ifelse(is.na(precip_date), 1, 0)
    
    precip_forloop <- precip_date[which(na_vector == 0)]
    
    
    #Retrieve daily lapse rate (following Thornton et al., 1997)
    
    lapse_stations_lower <- c(154, 147, 149, 146, 99, 100, 41, 42, 167, 44, 44, 134, 134, 134, 68, 68, 179, 168, 49, 101, 101, 129, 125, 178, 180, 171, 171, 176, 176, 176, 181, 65, 65, 96, 77, 159, 28, 35, 36, 110, 139, 199, 138, 102, 80, 80, 80, 66, 66, 66, 143, 72, 182, 187, 183, 205, 194)
    
    lapse_stations_higher <- c(155, 144, 145, 145, 74, 74, 161, 161, 161, 179, 83, 179, 83, 177, 177, 83, 83, 169, 169, 129, 8, 8, 135, 172, 173, 173, 175, 181, 174, 107, 174, 111, 62, 111, 106, 15, 15, 34, 34, 104, 88, 204, 114, 103, 142, 90, 143, 142, 143, 90, 90, 70, 202, 202, 203, 201, 195)
    
    
    precip_date_inp <- ifelse(precip_date <= 0, NA, precip_date)
    
    elev_diff_reg <- abs(elev[lapse_stations_higher] - elev[lapse_stations_lower])/1000
    
    elev_diff_reg <- ifelse(elev_diff_reg > 0.2, elev_diff_reg, NA)
    
    
    #reg test
    
    precip_reg <- (precip_date_inp[lapse_stations_higher] - precip_date_inp[lapse_stations_lower])/(precip_date_inp[lapse_stations_higher] + precip_date_inp[lapse_stations_lower])
    
    #precip_reg <- (precip_date_inp[lapse_stations_higher]/precip_date_inp[lapse_stations_lower])
    
    
    if (sum(precip_reg, na.rm = T) == 0 || sum(!is.na(precip_reg*elev_diff_reg)) < 3) {
      
      lapse_rate_vector[o] <- NA
      
      lapse_uncertainty_vector[o] <- NA
      
    } else {
      
      lapse_model <- lm(precip_reg ~ elev_diff_reg)
      
      lapse_model_summary <- summary(lapse_model)
      
      lapse_rate <- round(lapse_model_summary$coefficients[2,1], 2)
      
      lapse_rate_vector[o] <- lapse_rate
      
      lapse_uncertainty <- round(lapse_model_summary$coefficients[2,2], 2)
      
      lapse_uncertainty <- ifelse(lapse_uncertainty < 0, 0, lapse_uncertainty)
      
      lapse_uncertainty_vector[o] <- lapse_uncertainty
      
    }
    
  }  
  
  
  lapse_rate_vector[1] <- 0.25
  lapse_uncertainty_vector[1] <- 0.25
  
  lapse_rate_vector <- ifelse(lapse_rate_vector > 1, NA, lapse_rate_vector)
  lapse_rate_vector <- ifelse(lapse_rate_vector < 0, NA, lapse_rate_vector)
  
  lapse_uncertainty_vector <- ifelse(lapse_uncertainty_vector > 1, NA, lapse_uncertainty_vector)
  lapse_uncertainty_vector <- ifelse(lapse_uncertainty_vector < 0, NA, lapse_uncertainty_vector)
  
  
  for (o in 2:(which(substr(days, 1, 10) == "2020-09-30"))) {
    
    if(is.na(lapse_rate_vector[o])) {
      
      lapse_rate_vector[o] <- lapse_rate_vector[o-1]
      
    } 
    
    
    if(is.na(lapse_uncertainty_vector[o])) {
      
      lapse_uncertainty_vector[o] <- lapse_uncertainty_vector[o-1]
      
    } 
    
    
  }
  
  
  #Choose period
  
  start <- "2019-10-01"
  end <- "2020-09-30"
  
  monthly_precip <- srtm
  monthly_precip[monthly_precip > 0] <- 0
  
  monthly_sd <- srtm
  monthly_sd[monthly_sd > 0] <- 0
  
  pred_rmse_mm <- vector()
  
  
  for (o in (which(substr(days, 1, 10) == start)):(which(substr(days, 1, 10) == end))) {
    
    
    #Retrieve station information
    
    lat <- station_info$lat
    
    long <- station_info$long
    
    elev <- station_info$elev
    
    
    #Calculate the number of NAs for all stations for this particular date
    
    precip_date <- precip_input[o,]
    
    na_vector <- ifelse(is.na(precip_date), 1, 0)
    
    precip_forloop <- precip_date[which(na_vector == 0)]
    
    
    if (sum(precip_forloop > 0) > 1) {
    
  
    #Lapse precipitation to minimum elevation
    
    lapse_rate <- lapse_rate_vector[o]
    lapse_uncertainty <- lapse_uncertainty_vector[o]
    
    lat <- lat[which(na_vector == 0)]
    long <- long[which(na_vector == 0)]
    elev <- elev[which(na_vector == 0)]
    
    
    #Transform precipitation from mm to normal distribution units
    
    order_p <- order(precip_forloop)
    
    precip_vector <- round(precip_forloop[order_p], 2)
    
    #Find duplicated precipitation
    
    duplica_precip <- ifelse(duplicated(precip_vector), 1, 0)
    
    duplica_index <- which(abs(diff(duplica_precip)) == 1)
    
    #Make the transformation to normal distribution units
    
    normal_dist <- qnorm(seq(1/length(precip_forloop), 1, by = 1/length(precip_forloop)), mean = 0, sd = 1)
    
    normal_dist[length(normal_dist)] <- ifelse(is.infinite(tail(normal_dist,1)),
                                               ((normal_dist[length(normal_dist)-1] -
                                                   normal_dist[length(normal_dist)-2]) +
                                                  normal_dist[length(normal_dist)-1]),
                                               tail(normal_dist,1))
    
    precip_normal_dist <- normal_dist
    
    
    #Remove duplicated precipitation
    
    if (length(duplica_index) > 1) {
      
      for (n in 0:((length(duplica_index)/2)-1)) {
        
        precip_normal_dist[duplica_index[n*2 + 1]:duplica_index[n*2 + 2]] <-
          rep(median(normal_dist[duplica_index[n*2 + 1]:duplica_index[n*2 + 2]]),
              length(normal_dist[duplica_index[n*2 + 1]:duplica_index[n*2 + 2]]))
        
      }
      
    }
    
    
    #Check the transformation
    
    nontrans <- ecdf(precip_vector)
    
    trans <- ecdf(precip_normal_dist)
    
    
    #Create normarlly-distributed data frame
    
    lat <- lat[order_p]
    long <- long[order_p]
    elev <- elev[order_p]
    
    normal_df <- data.frame(lat, long, elev, precip_normal_dist)
    
    colnames(normal_df) <- c("lat", "long", "elev", "precip_normal_dist")
    
    coordinates(normal_df) =~long + lat
    
    
    #Remove duplicated stations
    
    points <- SpatialPoints(normal_df, proj4string = CRS("+init=epsg:4326"))
    
    if (length(zerodist(points)[,1]) == 0) {
      
      normal_df <- normal_df
      
    } else {
      
      normal_df <- normal_df[-c(zerodist(points)[,1]),]
      
    }
    
    
    #Create horizontal variogram
    
    options(warn = 2)
    
    try_Gau <- try(autofitVariogram(precip_normal_dist ~ 1, normal_df, model = c("Gau")), silent = T)
    try_Exp <- try(autofitVariogram(precip_normal_dist ~ 1, normal_df, model = c("Exp")), silent = T)
    try_Sph <- try(autofitVariogram(precip_normal_dist ~ 1, normal_df, model = c("Sph")), silent = T)
    try_Pen <- try(autofitVariogram(precip_normal_dist ~ 1, normal_df, model = c("Pen")), silent = T)
    
    options(warn = 1)
    
    try_vec <- c(substr(try_Gau[1], 1, 5), substr(try_Exp[1], 1, 5), substr(try_Sph[1], 1, 5), substr(try_Pen[1], 1, 5))
    
    var_models <- c("Gau", "Exp", "Sph", "Pen")
    
    vgm <- autofitVariogram(precip_normal_dist ~ 1, normal_df, model = var_models[which(try_vec != "Error")])
    
    vgm_model <- vgm$var_model
    
    
    #Predict horizontal precipitation in Z-score units
    
    pred.reg <- krige(precip_normal_dist ~ 1, normal_df, srtm_points, model = vgm_model, maxdist = 2, nmax = 8)
    
    krig_output <- as.data.frame(pred.reg)
    
    
    krig_precip <- data.frame(krig_output$x, krig_output$y, krig_output$var1.pred)
    
    colnames(krig_precip) <- c("x", "y", "z")
    
    precip_raster_hor <- rasterFromXYZ(krig_precip)
    
    
    
    krig_var <- data.frame(krig_output$x, krig_output$y, krig_output$var1.var)
    
    colnames(krig_var) <- c("x", "y", "z")
    
    var_raster_hor <- rasterFromXYZ(krig_var)
    
    
    #Indicator kriging for zero precipitation
    
    precip_zero <- ifelse(precip_vector < 0.2, 0, 1) #0.2 mm is a trace value for ECCC
    
    precip_zero_df <- data.frame(lat, long, elev, precip_zero)
    
    colnames(precip_zero_df) <- c("lat", "long", "elev", "precip_zero")
    
    coordinates(precip_zero_df) =~long + lat
    
    
    if (sum(precip_zero_df$precip_zero) == 0) {
      
      krig_zero_raster_std <- srtm
      krig_zero_raster_std <- setValues(krig_zero_raster_std, 0)
      
    } else if ((sum(precip_zero_df$precip_zero) == length(precip_zero_df$precip_zero))) {
      
      krig_zero_raster_std <- srtm
      krig_zero_raster_std <- setValues(krig_zero_raster_std, 1)
      
    } else {
      
      points <- SpatialPoints(precip_zero_df, proj4string = CRS("+init=epsg:4326"))
      
      if (length(zerodist(points)[,1]) == 0) {
        
        precip_zero_df <- precip_zero_df
        
      } else {
        
        precip_zero_df <- precip_zero_df[-c(zerodist(points)[,1]),]
        
      }
      
      
      options(warn = 2)
      
      try_Gau <- try(autofitVariogram(precip_zero ~ 1, precip_zero_df, model = c("Gau")), silent = T)
      try_Exp <- try(autofitVariogram(precip_zero ~ 1, precip_zero_df, model = c("Exp")), silent = T)
      try_Sph <- try(autofitVariogram(precip_zero ~ 1, precip_zero_df, model = c("Sph")), silent = T)
      try_Pen <- try(autofitVariogram(precip_zero ~ 1, precip_zero_df, model = c("Pen")), silent = T)
      
      options(warn = 1)
      
      try_vec <- c(substr(try_Gau[1], 1, 5), substr(try_Exp[1], 1, 5), substr(try_Sph[1], 1, 5), substr(try_Pen[1], 1, 5))
      
      var_models <- c("Gau", "Exp", "Sph", "Pen")
      
      
      vgm_zero <- autofitVariogram(precip_zero ~ 1, precip_zero_df, model = var_models[which(try_vec != "Error")])
      
      vgm_model_zero <- vgm_zero$var_model
      
      
      pred.reg_zero <- krige(precip_zero ~ 1, precip_zero_df, srtm_points, model = vgm_model_zero, maxdist = 2, nmax = 8)
      
      krig_output_zero <- as.data.frame(pred.reg_zero)
      
      krig_zero <- data.frame(krig_output_zero$x, krig_output_zero$y, krig_output_zero$var1.pred)
      
      colnames(krig_zero) <- c("x", "y", "z")
      
      krig_zero_raster <- rasterFromXYZ(krig_zero)
      
      krig_zero_raster[krig_zero_raster <= 0.45] <- 0
      krig_zero_raster[krig_zero_raster > 0.55] <- 1
      
      krig_zero_raster_std <- (krig_zero_raster - 0.45)/(0.55 - 0.45)
      
      krig_zero_raster_std[krig_zero_raster_std < 0] <- 0
      krig_zero_raster_std[krig_zero_raster_std > 1] <- 1
      
    }
    
    
    #Interpolate station elevations
    
    error_msg <- try(autofitVariogram(precip_zero_df$elev ~ 1, precip_zero_df, model = c("Lin")), silent = T)
    
    if (substr(error_msg[1], 1, 5) == "Error") {
    
    vgm_zero <- autofitVariogram(precip_zero_df$elev ~ 1, precip_zero_df, model = c("Sph"))
    
    vgm_model_zero <- vgm_zero$var_model
    
    interp_elev <- krige(precip_zero_df$elev ~ 1, precip_zero_df, srtm_points, model = vgm_model_zero, maxdist = 2, nmax = 8)
    
    } else {
    
    vgm_zero <- autofitVariogram(precip_zero_df$elev ~ 1, precip_zero_df, model = c("Lin"))
    
    vgm_model_zero <- vgm_zero$var_model
    
    interp_elev <- krige(precip_zero_df$elev ~ 1, precip_zero_df, srtm_points, model = vgm_model_zero, maxdist = 2, nmax = 8)
    
    }
    
    
    interp_elev_df <- as.data.frame(interp_elev)
    
    interp_elev_df_xyz <- data.frame(interp_elev_df$x, interp_elev_df$y, interp_elev_df$var1.pred)
    
    colnames(interp_elev_df_xyz) <- c("x", "y", "z")
    
    interp_elev_final <- rasterFromXYZ(interp_elev_df_xyz)
    
    
    #Calculate lapse rates for elevation model
    
    dif_srtm <- (srtm/1000 - interp_elev_final/1000)
    
    
    f_lapse <- lapse_rate*dif_srtm
    
    f_unc <- lapse_uncertainty*dif_srtm
    
    f_lapse[f_lapse < 0] <- 0
    
    f_lapse[f_lapse > 0.95] <- 0.95
    
    f_unc[f_unc < 0] <- 0
    
    f_unc[f_unc > 0.95] <- 0.95
    
    
    pixel_lapse_multiplier <- (1 + f_lapse)/(1 - f_lapse)
    
    pixel_lapse_unc_multiplier <- (1 + f_unc)/(1 - f_unc)
    
    
    #The multiplier cap of 8 was based on the Fisera Ridge and Hay Meadow lapse rate times 2, because there are up to 2 km of dif_srtm
    
    lapse_reference <- precip_input[,174]/precip_input[,176]
    
    lapse_reference <- ifelse(is.infinite(lapse_reference), NA, lapse_reference)
    
    multiplier_cap <- round(mean(lapse_reference, na.rm = T))*cellStats(dif_srtm, max, na.rm = T)
    
    
    pixel_lapse_multiplier[pixel_lapse_multiplier > multiplier_cap] <- multiplier_cap
    
    #pixel_lapse_unc_multiplier[pixel_lapse_unc_multiplier < 0] <- multiplier_cap
    
    pixel_lapse_unc_multiplier[pixel_lapse_unc_multiplier > multiplier_cap] <- multiplier_cap
    
    #pixel_lapse_unc_multiplier[pixel_lapse_unc_multiplier < 0] <- multiplier_cap
    
    
    #Create reclassification matrix
    
    precip_mm_unique <- unique(precip_vector)
    
    precip_mm_diff <- diff(precip_mm_unique)
    
    precip_dist_unique <- unique(precip_normal_dist)
    
    
    precip_dist_extra <- vector()
    
    for (m in 2:length(precip_mm_unique)) {
      
      precip_dist_extra <- append(precip_dist_extra, seq(precip_dist_unique[m-1], precip_dist_unique[m], abs(precip_dist_unique[m-1] - precip_dist_unique[m])/(precip_mm_diff[m-1]*100)))
      
    }
    
    
    precip_mm_extra <- vector()
    
    for (m in 2:length(precip_mm_unique)) {
      
      precip_mm_extra <- append(precip_mm_extra, seq(precip_mm_unique[m-1], precip_mm_unique[m], abs(precip_mm_unique[m-1] - precip_mm_unique[m])/(precip_mm_diff[m-1]*100)))
      
    }
    
    
    rec_matrix <- cbind(c(precip_dist_extra[1], precip_dist_extra[-c(length(precip_dist_extra))]), precip_dist_extra, precip_mm_extra)
    
    
    #Prediction RMSE
    
    pred_rmse_mm[o] <- quantile(precip_vector, trans(vgm$sserr/length(normal_df$precip_normal_dist)))
    
    
    #Interpolate precipitation
    
    precip_mm <- reclassify(precip_raster_hor, rec_matrix)
    
    precip_mm <- precip_mm*krig_zero_raster_std
    
    precip_lapsed_mm <- precip_mm*pixel_lapse_multiplier
    
    precip_lapsed_mm[precip_mm < 0.2] <- 0 #Based on 0.2 trace value of ECCC
    
    precip_lapsed_mm[precip_lapsed_mm < 0] <- 0
    
    
    uncertainty_mm <- reclassify(sqrt(var_raster_hor), rec_matrix)
    
    uncertainty_mm <- uncertainty_mm*krig_zero_raster_std
    
    uncertainty_lapsed_mm <- uncertainty_mm*pixel_lapse_unc_multiplier
    
    uncertainty_lapsed_mm[precip_mm < 0.2] <- 0
    
    uncertainty_lapsed_mm[uncertainty_lapsed_mm < 0] <- 0
    
    
    monthly_precip <- monthly_precip + precip_lapsed_mm
    
    monthly_sd <- monthly_sd + uncertainty_lapsed_mm
    
    
    } else {
      
      precip_lapsed_mm <- srtm
      precip_lapsed_mm <- setValues(precip_lapsed_mm, 0)
      
      uncertainty_lapsed_mm <- srtm
      uncertainty_lapsed_mm <- setValues(uncertainty_lapsed_mm, 0)
      
      monthly_precip <- monthly_precip + precip_lapsed_mm
      
      monthly_sd <- monthly_sd + uncertainty_lapsed_mm
      
    }
      
    
  }
  
  
  writeRaster(monthly_precip, filename = paste0("/media/project/abertoncini/01_Precip_Network/08_Chinook_Runs_20230724/OUTPUTS/monthly_precip_20230628_", substr(end, 1, 4), "_q", q), format = "GTiff", overwrite = T)
  
  writeRaster(monthly_sd, filename = paste0("/media/project/abertoncini/01_Precip_Network/08_Chinook_Runs_20230724/OUTPUTS/monthly_sd_20230628_", substr(end, 1, 4), "_q", q), format = "GTiff", overwrite = T)
  
  
}

stopCluster(cl)

}

