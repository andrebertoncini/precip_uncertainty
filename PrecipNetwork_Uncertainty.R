library(raster)
library(ggplot2)
library(gstat)
library(automap)
library(parallel)
library(doParallel)
library(RColorBrewer)
library(MASS)


setwd("/home/aberton/projects/rpp-kshook/aberton/precip_network/run_20220413")


#Load precipitation data

precip_input <- read.csv("Precip_Undercatch_Corrected_20210727.csv")[,-c(1)]

precip_input <- ifelse(as.matrix(precip_input) > 160, NA, as.matrix(precip_input)) #based on climate normals

#Load station information

station_info <- read.csv("Station_Info_Inside_Domain_20210727.csv")[,-c(1)]


#Load quadrant information

quad_coord <- read.csv("Quadrants_Coordinates_20220416.csv", sep = ";")[,-c(1)]


#Load DEM

srtm <- raster("SRTM_90m_VoidFilled_GEE_v1_Mosaic_Clip.tif")

srtm[srtm == -32768] <- NA
srtm[srtm == 32767] <- NA


#Parallelization by spatial quadrant

no_cores <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

cl <- makeCluster(no_cores, type = "PSOCK")

registerDoParallel(cl)

foreach(q = 97:191) %dopar% {


library(raster)
library(ggplot2)
library(gstat)
library(automap)
library(parallel)
library(doParallel)
library(RColorBrewer)
library(MASS)


#Crop DEM for faster processing

crop_extent <- extent(c(quad_coord$long_min[q], quad_coord$long_max[q], quad_coord$lat_min[q], quad_coord$lat_max[q]))

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
  
  
  precip_reg <- (precip_date[lapse_stations_higher] - precip_date[lapse_stations_lower])
  
  precip_reg <- ifelse(precip_reg < 0.01, NA, precip_reg)
  
  elev_diff_reg <- abs(elev[lapse_stations_higher] - elev[lapse_stations_lower])/1000
  
  elev_diff_reg <- ifelse(elev_diff_reg > 0.2, elev_diff_reg, NA)
  
  
  if (sum(precip_reg, na.rm = T) == 0 || sum(!is.na(precip_reg)) < 3) {
    
    lapse_rate <- lapse_rate_vector[o-1]
    
    lapse_rate_vector[o] <- lapse_rate
    
    lapse_uncertainty <- lapse_uncertainty_vector[o-1]
    
    lapse_uncertainty_vector[o] <- lapse_uncertainty
    
  } else {
    
    lapse_rate <- mean(elev_diff_reg/precip_reg, na.rm = T)
    
    lapse_rate <- ifelse(lapse_rate < 0, 0, lapse_rate)
    
    #lapse_rate <- ifelse(lapse_rate > 0.4, 0.4, lapse_rate)
    
    lapse_rate_vector[o] <- lapse_rate
    
    lapse_uncertainty <- sd(elev_diff_reg/precip_reg, na.rm = T)
    
    lapse_uncertainty <- ifelse(lapse_uncertainty < 0, 0, lapse_uncertainty)
    
    #lapse_uncertainty <- ifelse(lapse_uncertainty > 0.4, 0.4, lapse_uncertainty)
    
    lapse_uncertainty_vector[o] <- lapse_uncertainty
    
  }
  
}  


#Choose period

start <- "2019-10-01"
end <- "2020-09-30"

pb <- txtProgressBar(min = which(substr(days, 1, 10) == start), max = which(substr(days, 1, 10) == end), style = 3)

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


#Lapse precipitation to minimum elevation

lapse_rate <- lapse_rate_vector[o]
lapse_uncertainty <- lapse_uncertainty_vector[o]

lat <- lat[which(na_vector == 0)]
long <- long[which(na_vector == 0)]
elev <- elev[which(na_vector == 0)]


difference <- ((mean(elev, na.rm = T)/1000) - (elev/1000))

neg_multiplier <- (1.0 - lapse_rate*(abs((mean(elev, na.rm = T)/1000) - (elev/1000))))
neg_multiplier <- ifelse(neg_multiplier < 0, 0.01, neg_multiplier)
neg_multiplier <- ifelse(neg_multiplier > 1, 1, neg_multiplier)

pos_multiplier <- (1.0 + lapse_rate*(abs((mean(elev, na.rm = T)/1000) - (elev/1000))))
pos_multiplier <- ifelse(pos_multiplier < 1, 1, pos_multiplier)

precip_forloop <- ifelse(difference < 0, precip_forloop*neg_multiplier, precip_forloop*pos_multiplier)


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
#plot(nontrans, main = "Nontransformed Daily Precipitation", ylab = "CDF [ ]", xlab = "Daily Precipitaiton [mm]")

trans <- ecdf(precip_normal_dist)
#plot(trans, main = "Transformed Daily Precipitation", ylab = "CDF [ ]", xlab = "Daily Precipitaiton [ ]")


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

vgm <- autofitVariogram(precip_normal_dist ~ 1, normal_df, model = c("Gau", "Exp", "Sph", "Log", "Pen"))

vgm_model <- vgm$var_model


#Predict horizontal precipitation in Z-score units

pred.reg <- krige(precip_normal_dist ~ 1, normal_df, srtm_points, model = vgm_model, maxdist = 2, nmax = 8)

if (length(warnings()) > 0 && substr(warnings()[length(warnings())], 1, 13) == "predict.gstat") {
  
  vgm <- autofitVariogram(precip_normal_dist ~ 1, normal_df, model = c("Mat"))
  
  vgm_model <- vgm$var_model
  
  pred.reg <- krige(precip_normal_dist ~ 1, normal_df, srtm_points, model = vgm_model, maxdist = 2, nmax = 8)
  
}


krig_output <- as.data.frame(pred.reg)


krig_precip <- data.frame(krig_output$x, krig_output$y, krig_output$var1.pred)

colnames(krig_precip) <- c("x", "y", "z")

precip_raster_hor <- rasterFromXYZ(krig_precip)

palette_1 <- brewer.pal(9, "YlGnBu")

#plot(precip_raster_hor, col = palette_1)


krig_var <- data.frame(krig_output$x, krig_output$y, krig_output$var1.var)

colnames(krig_var) <- c("x", "y", "z")

var_raster_hor <- rasterFromXYZ(krig_var)

palette_2 <- brewer.pal(11, "RdYlGn")

#plot(var_raster_hor, col = rev(palette_2))


#Indicator kriging for zero precipitation

precip_zero <- ifelse(precip_vector < 0.1, 0, 1)

precip_zero_df <- data.frame(lat, long, elev, precip_zero)

colnames(precip_zero_df) <- c("lat", "long", "elev", "precip_zero")

if (sum(precip_zero_df$precip_zero) == 0) {
  
  krig_zero_raster_std <- srtm
  krig_zero_raster_std[krig_zero_raster_std > 0] <- 0
  
} else {

coordinates(precip_zero_df) =~long + lat


  if (length(zerodist(points)[,1]) == 0) {

    precip_zero_df <- precip_zero_df

  } else {

    precip_zero_df <- precip_zero_df[-c(zerodist(points)[,1]),]

  }


vgm_zero <- autofitVariogram(precip_zero ~ 1, precip_zero_df, model = paste(vgm_model$model[2]))

vgm_model_zero <- vgm_zero$var_model


pred.reg_zero <- krige(precip_zero ~ 1, precip_zero_df, srtm_points, model = vgm_model_zero, maxdist = 2, nmax = 8)

if (length(warnings()) > 0 && substr(warnings()[length(warnings())], 1, 13) == "predict.gstat") {
  
  vgm_zero <- autofitVariogram(precip_zero ~ 1, precip_zero_df, model = c("Mat"))
  
  vgm_model_zero <- vgm_zero$var_model
  
  pred.reg_zero <- krige(precip_zero ~ 1, precip_zero_df, srtm_points, model = vgm_model_zero, maxdist = 2, nmax = 8)
  
}


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


#Create stations spatial points

station_with_data <- station_info[which(na_vector == 0),]

station_x <- station_with_data$long

station_y <- station_with_data$lat

station_z <- station_with_data$elev


station_xyz <- data.frame(station_x, station_y, station_z)

names(station_xyz) <- c("x", "y", "z")

coordinates(station_xyz) <-  ~ x + y


#Plot to monitor progress

palette_1 <- brewer.pal(9, "YlGnBu")

#plot(monthly_precip, col = palette_1)


#Interpolate precipitation using uncertainty weights

srtm_points <- as.data.frame(rasterToPoints(srtm))

names(srtm_points) <- c("x", "y", "z")

coordinates(srtm_points) <- ~ x + y


pixel_lapse_multiplier <- vector(length = nrow(srtm_points))

pixel_lapse_unc_multiplier <- vector(length = nrow(srtm_points))

  for (l in 1:nrow(srtm_points)) {

    interp_hor <- sqrt((abs(srtm_points@coords[l,2] - lat)*110.574)^2 + (abs(srtm_points@coords[l,1] - long)*(111.320*cos(51.27304*(pi/180))))^2)

    interp_ver <- (abs(srtm_points$z[l] - elev)/1000)

    interp_hor_nona <- interp_hor[which(na_vector == 0)]

    interp_ver_nona <- interp_ver[which(na_vector == 0)]

    lat_nona <- lat[which(na_vector == 0)]
    long_nona <- long[which(na_vector == 0)]
    elev_nona <- elev[which(na_vector == 0)]


    #Bring station list to original order

    order_e <- order(interp_hor_nona)


    #Interpolate elevation using bilinear interpolation of nearest 12 stations

    dist_elev <- interp_hor_nona[order_e][1:8]

    elev4 <- (elev_nona[order_e][1:8])/1000

    interp_elev <- sum((dist_elev*elev4), na.rm = T)/sum(dist_elev, na.rm = T)

    
    #Fixing for instability in the lapse rate equation for high lapse rates and elevation differences
    
    difference_2 <- ((srtm_points$z[l]/1000) - (mean(elev, na.rm = T)/1000))
    
    neg_multiplier_2 <- (1.0 - lapse_rate*(abs(((srtm_points$z[l]/1000) - (mean(elev, na.rm = T)/1000)))))
    neg_multiplier_2 <- ifelse(neg_multiplier_2 < 0, 0.01, neg_multiplier_2)
    neg_multiplier_2 <- ifelse(neg_multiplier_2 > 1, 1, neg_multiplier_2)
    
    pos_multiplier_2 <- (1.0 + lapse_rate*(abs(((srtm_points$z[l]/1000) - (mean(elev, na.rm = T)/1000)))))
    pos_multiplier_2 <- ifelse(pos_multiplier_2 < 1, 1, pos_multiplier_2)
    
    pixel_lapse_multiplier[l] <- ifelse(difference_2 < 0, neg_multiplier_2, pos_multiplier_2)
    

    difference_3 <- ((srtm_points$z[l]/1000) - interp_elev)
    
    neg_multiplier_3 <- (1.0 - lapse_uncertainty*(abs(((srtm_points$z[l]/1000) - interp_elev))))
    neg_multiplier_3 <- ifelse(neg_multiplier_3 < 0, 0.01, neg_multiplier_3)
    neg_multiplier_3 <- ifelse(neg_multiplier_3 > 1, 1, neg_multiplier_3)
    
    pos_multiplier_3 <- (1.0 + lapse_uncertainty*(abs(((srtm_points$z[l]/1000) - interp_elev))))
    pos_multiplier_3 <- ifelse(pos_multiplier_3 < 1, 1, pos_multiplier_3)
    
    pixel_lapse_unc_multiplier[l] <- ifelse(difference_3 < 0, neg_multiplier_3, pos_multiplier_3)
    
    }


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

precip_df <- data.frame(srtm_points@coords, pixel_lapse_multiplier)

colnames(precip_df) <- c("x", "y", "z")

precip_lapse_multi <- rasterFromXYZ(precip_df)

palette_1 <- brewer.pal(9, "YlGnBu")

precip_mm <- reclassify(precip_raster_hor, rec_matrix)

precip_mm <- precip_mm*krig_zero_raster_std

precip_mm[precip_mm < 0] <- 0

precip_lapsed_mm <- precip_mm*precip_lapse_multi

precip_lapsed_mm[precip_lapsed_mm < 0] <- 0

#plot(precip_lapsed_mm, col = palette_1)


precip_df_var <- data.frame(srtm_points@coords, pixel_lapse_unc_multiplier)

colnames(precip_df_var) <- c("x", "y", "z")

precip_lapse_sd_multi <- rasterFromXYZ(precip_df_var)

uncertainty_mm <- reclassify(sqrt(var_raster_hor), rec_matrix)

uncertainty_lapsed_mm <- uncertainty_mm*precip_lapse_sd_multi

uncertainty_lapsed_mm[uncertainty_lapsed_mm < 0] <- 0

#plot(uncertainty_lapsed_mm, col = rev(palette_2))


monthly_precip <- monthly_precip + precip_lapsed_mm

monthly_sd <- monthly_sd + uncertainty_lapsed_mm


#plot(monthly_precip, col = palette_1)


Sys.sleep(0.0001)
setTxtProgressBar(pb, o)

  }


writeRaster(monthly_precip, filename = paste0("monthly_precip_", substr(end, 1, 4), "_q", q), format = "GTiff", overwrite = T)

writeRaster(monthly_sd, filename = paste0("monthly_sd_", substr(end, 1, 4), "_q", q), format = "GTiff", overwrite = T)


}

stopCluster(cl)


#plot(monthly_precip, col = palette_1)
#plot(station_xyz, add = T)

#plot(monthly_sd, col = rev(palette_2))
#plot(station_xyz, add = T)


#Calculate CV

#cv <- monthly_sd/monthly_precip

#plot(cv, col = rev(palette_2))
#plot(station_xyz, add = T)


#Export lapse rate, lapse uncertainty, and RMSE

#lapse_export <- data.frame(lapse_rate_vector[(which(substr(days, 1, 10) == start)):(which(substr(days, 1, 10) == end))],
#                lapse_uncertainty_vector[(which(substr(days, 1, 10) == start)):(which(substr(days, 1, 10) == end))],
#                pred_rmse_mm[(which(substr(days, 1, 10) == start)):(which(substr(days, 1, 10) == end))])

#colnames(lapse_export) <- c("lapse_rate", "lapse_uncertainty", "pred_rmse_mm")

#write.csv(lapse_export, file = paste0("lapse_export_", substr(end, 1, 4), ".csv"))
