library(raster)
library(ggplot2)
library(gstat)
library(automap)
library(parallel)
library(doParallel)
library(RColorBrewer)
library(MASS)


setwd("/media/project/abertoncini/01_Precip_Network/01_3D_Variogram_20220125")

#Load precipitation data

precip_input <- read.csv("Precip_Undercatch_Corrected_20210727.csv")[,-c(1)]

precip_input <- ifelse(as.matrix(precip_input) > 160, NA, as.matrix(precip_input)) #based on climate normals

#Load station information

station_info <- read.csv("Station_Info_Inside_Domain_20210727.csv")[,-c(1)]


#Load DEM

srtm <- raster("SRTM_90m_Tile_6.tif")

srtm[srtm == -32768] <- NA
srtm[srtm == 32767] <- NA

#Crop DEM for faster processing

crop_extent <- extent(c(-115.3, -115.2, 50.83, 50.89))

srtm <- crop(srtm, crop_extent)

plot(srtm)

#Create data.frame with DEM xyz

srtm_points <- as.data.frame(rasterToPoints(srtm))

names(srtm_points) <- c("x", "y", "z")

coordinates(srtm_points) <- ~ x + y


#Select day to calculate variogram

#Define available days

days <- seq(as.POSIXct("1990-10-01", format = "%Y-%m-%d", tz = "GMT"),
            as.POSIXct("2020-09-30", format = "%Y-%m-%d", tz = "GMT"), by = 86400)


#Choose period

pb <- txtProgressBar(min = which(substr(days, 1, 10) == "2000-09-29"), max = which(substr(days, 1, 10) == "2001-09-30"), style = 3)

monthly_precip <- stack()

monthly_sd <- stack()

lapse_rate_vector <- vector()

lapse_uncertainty_vector <- vector()

for (o in (which(substr(days, 1, 10) == "2019-09-30")):(which(substr(days, 1, 10) == "2020-09-30"))) {


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

precip_reg <- ifelse(precip_reg <= 0, NA, precip_reg)

elev_diff_reg <- abs(elev[lapse_stations_higher] - elev[lapse_stations_lower])/1000

elev_diff_reg <- ifelse(elev_diff_reg > 0.2, elev_diff_reg, NA)


if (sum(precip_reg, na.rm = T) == 0 || sum(!is.na(precip_reg)) < 3) {

  lapse_rate <- lapse_rate_vector[o-1]

  lapse_rate_vector[o] <- lapse_rate

  lapse_uncertainty <- lapse_uncertainty_vector[o-1]

  lapse_uncertainty_vector[o] <- lapse_uncertainty

} else {
  
  #lapse_rate <- mean(elev_diff_reg/precip_reg, na.rm = T)
  
  lapse_rate <- mean(elev_diff_reg, na.rm = T)/mean(precip_reg, na.rm = T)

  lapse_rate <- ifelse(lapse_rate < 0, 0, lapse_rate)

  lapse_rate <- ifelse(lapse_rate > 0.4, 0.4, lapse_rate)

  lapse_rate_vector[o] <- lapse_rate

  #lapse_uncertainty <- sd(elev_diff_reg/precip_reg, na.rm = T)
  
  lapse_uncertainty <- mean(elev_diff_reg, na.rm = T)/sd(precip_reg, na.rm = T)

  lapse_uncertainty <- ifelse(lapse_uncertainty < 0, 0, lapse_uncertainty)

  lapse_uncertainty <- ifelse(lapse_uncertainty > 0.4, 0.4, lapse_uncertainty)

  lapse_uncertainty_vector[o] <- lapse_uncertainty

}


#Lapse precipitation to minimum elevation

lat <- lat[which(na_vector == 0)]
long <- long[which(na_vector == 0)]
elev <- elev[which(na_vector == 0)]

precip_forloop <- ((1 + (lapse_rate*((min(elev, na.rm = T)/1000) - (elev/1000))))/(1 - (lapse_rate*((min(elev, na.rm = T)/1000) - (elev/1000)))))*precip_forloop


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

if (length(duplica_index > 0)) {

for (n in 0:((length(duplica_index)/2)-1)) {

  precip_normal_dist[duplica_index[n*2 + 1]:duplica_index[n*2 + 2]] <-
    rep(median(normal_dist[duplica_index[n*2 + 1]:duplica_index[n*2 + 2]]),
        length(normal_dist[duplica_index[n*2 + 1]:duplica_index[n*2 + 2]]))

  }
  
}


#Check the transformation

nontrans <- ecdf(precip_vector)
plot(nontrans, main = "Nontransformed Daily Precipitation", ylab = "CDF [ ]", xlab = "Daily Precipitaiton [mm]")

trans <- ecdf(precip_normal_dist)
plot(trans, main = "Transformed Daily Precipitation", ylab = "CDF [ ]", xlab = "Daily Precipitaiton [ ]")


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

vgm <- autofitVariogram(precip_normal_dist ~ 1, normal_df, model = c("Gau", "Exp", "Sph", "Lin"))

vgm_model <- vgm$var_model

#Predict horizontal precipitation in Z-score units

pred.reg <- krige(precip_normal_dist ~ 1, normal_df, srtm_points, model = vgm_model, maxdist = 2, nmax = 12)

krig_output <- as.data.frame(pred.reg)


krig_precip <- data.frame(krig_output$x, krig_output$y, krig_output$var1.pred)

colnames(krig_precip) <- c("x", "y", "z")

precip_raster_hor <- rasterFromXYZ(krig_precip)

palette_1 <- brewer.pal(9, "YlGnBu")

plot(precip_raster_hor, col = palette_1)


krig_var <- data.frame(krig_output$x, krig_output$y, krig_output$var1.var)

colnames(krig_var) <- c("x", "y", "z")

var_raster_hor <- rasterFromXYZ(krig_var)

palette_2 <- brewer.pal(11, "RdYlGn")

plot(var_raster_hor, col = rev(palette_2))


#Indicator kriging for zero precipitation

precip_zero <- ifelse(precip_vector < 0.2, 0, 1)

precip_zero_df <- data.frame(lat, long, elev, precip_zero)

colnames(precip_zero_df) <- c("lat", "long", "elev", "precip_zero")

coordinates(precip_zero_df) =~long + lat


  if (length(zerodist(points)[,1]) == 0) {

    precip_zero_df <- precip_zero_df

  } else {

    precip_zero_df <- precip_zero_df[-c(zerodist(points)[,1]),]

  }


vgm_zero <- autofitVariogram(precip_zero ~ 1, precip_zero_df, model = paste(vgm_model$model[2]))

vgm_model_zero <- vgm_zero$var_model


pred.reg_zero <- krige(precip_zero ~ 1, precip_zero_df, srtm_points, model = vgm_model_zero, maxdist = 2, nmax = 12)

krig_output_zero <- as.data.frame(pred.reg_zero)

krig_zero <- data.frame(krig_output_zero$x, krig_output_zero$y, krig_output_zero$var1.pred)

colnames(krig_zero) <- c("x", "y", "z")

krig_zero_raster <- rasterFromXYZ(krig_zero)

krig_zero_raster[krig_zero_raster <= 0.45] <- 0
krig_zero_raster[krig_zero_raster > 0.55] <- 1

krig_zero_raster_std <- (krig_zero_raster - 0.45)/(0.55 - 0.45)

krig_zero_raster_std[krig_zero_raster_std < 0] <- 0
krig_zero_raster_std[krig_zero_raster_std > 1] <- 1


#Create stations spatial points

station_with_data <- station_info[which(na_vector == 0),]

station_x <- station_with_data$long

station_y <- station_with_data$lat

station_z <- station_with_data$elev


station_xyz <- data.frame(station_x, station_y, station_z)

names(station_xyz) <- c("x", "y", "z")

coordinates(station_xyz) <-  ~ x + y


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

    dist_elev <- interp_hor_nona[order_e][1:12]

    elev4 <- (elev_nona[order_e][1:12])/1000

    interp_elev <- sum((dist_elev*elev4), na.rm = T)/sum(dist_elev, na.rm = T)

    
    #Fixing for instability in the lapse rate equation for high lapse rates and elevation differences
    
    pixel_lapse_multiplier[l] <- ((1 + (lapse_rate*((srtm_points$z[l]/1000) - (min(elev, na.rm = T)/1000))))/(1 - (lapse_rate*((srtm_points$z[l]/1000) - (min(elev, na.rm = T)/1000)))))
      
    pixel_lapse_unc_multiplier[l] <- ((1 + (lapse_uncertainty*((srtm_points$z[l]/1000) - (min(elev, na.rm = T)/1000))))/(1 - (lapse_uncertainty*((srtm_points$z[l]/1000) - (min(elev, na.rm = T)/1000)))))

    
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


#Interpolate precipitation

precip_df <- data.frame(srtm_points@coords, pixel_lapse_multiplier)

colnames(precip_df) <- c("x", "y", "z")

precip_lapse_multi <- rasterFromXYZ(precip_df)

palette_1 <- brewer.pal(9, "YlGnBu")

precip_mm <- reclassify(precip_raster_hor, rec_matrix)

precip_mm <- precip_mm*krig_zero_raster_std

precip_lapsed_mm <- precip_mm*precip_lapse_multi

precip_lapsed_mm[precip_lapsed_mm < 0] <- 0

plot(precip_lapsed_mm, col = palette_1)


precip_df_var <- data.frame(srtm_points@coords, pixel_lapse_unc_multiplier)

colnames(precip_df_var) <- c("x", "y", "z")

precip_lapse_sd_multi <- rasterFromXYZ(precip_df_var)

uncertainty_mm <- reclassify(sqrt(var_raster_hor), rec_matrix)

uncertainty_lapsed_mm <- uncertainty_mm*precip_lapse_sd_multi

uncertainty_lapsed_mm[uncertainty_lapsed_mm < 0] <- 0

plot(uncertainty_lapsed_mm, col = rev(palette_2))


monthly_precip <- stack(monthly_precip, precip_lapsed_mm)

monthly_sd <- stack(monthly_sd, uncertainty_lapsed_mm)


Sys.sleep(0.0001)
setTxtProgressBar(pb, o)

}


plot(sum(monthly_precip[[-c(1)]]), col = palette_1)

plot(sum(monthly_sd[[-c(1)]]), col = rev(palette_2))


#Calculate CV

cv <- sum(monthly_sd[[-c(1)]])/sum(monthly_precip[[-c(1)]])

plot(cv, col = rev(palette_2))


#Export rasters

yearmonth <- "201910"

writeRaster(sum(monthly_precip), filename = paste0("/media/project/abertoncini/01_Precip_Network/02_Uncertainty_SampleArea_20220224/monthly_precip_", yearmonth), format = "GTiff", overwrite = T)

writeRaster(sum(monthly_sd), filename = paste0("/media/project/abertoncini/01_Precip_Network/02_Uncertainty_SampleArea_20220224/monthly_sd_", yearmonth), format = "GTiff", overwrite = T)

writeRaster(cv, filename = paste0("/media/project/abertoncini/01_Precip_Network/02_Uncertainty_SampleArea_20220224/monthly_cv_", yearmonth), format = "GTiff", overwrite = T)