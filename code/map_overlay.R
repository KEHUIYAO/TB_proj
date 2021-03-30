library(ggplot2)
library(rgdal)
library(raster)
library(ggmap)
library(sf)
library(sp)
library(classInt)
library(tmap)
library(leaflet)
library(RColorBrewer)
library(maptools)
library(readxl)

# this function is to create a spatial dataframe, which contains the location and the log hazard ratio for each location unit
# create_spatial_data_frame <- function(effect,lat,long,lat_band=10,long_band=10){
#   boundary = c(min(lat),max(lat),min(long),max(long))
#   x <- seq(1,lat_band)
#   y <- seq(1,long_band)
#   a = (boundary[2]-boundary[1]) / lat_band
#   b = (boundary[4]-boundary[3]) / long_band
#   data <- expand.grid(X=as.character(x), Y=as.character(y))
#   z = c()
#   id = c()
#   polygon_list = list()
#   for(i in x){
#     for (j in y){
#       id <- c(id,paste(i,",",j,sep = ""))
#       cur_val = c()
#       cur_bound = c(boundary[1]+(i-1)*a, 
#                     boundary[1]+ i*a,
#                     boundary[3]+(j-1)*a,
#                     boundary[3]+ j*a)
#       
#       cur_polygon = Polygon(rbind(c(cur_bound[1],cur_bound[3]),
#                                   c(cur_bound[2],cur_bound[3]),
#                                   c(cur_bound[2],cur_bound[4]),
#                                   c(cur_bound[1],cur_bound[4]),
#                                   c(cur_bound[1],cur_bound[3])))
#       
#       cur_polygon = Polygons(list(cur_polygon), ID = paste(i,",",j,sep = ""))
#       
#       #print(cur_polygon)
#       polygon_list = append(polygon_list,cur_polygon)
#       for (k in 1:length(lat)){
#         if (lat[k]<cur_bound[2] & lat[k]>cur_bound[1] & long[k]<cur_bound[4] & long[k]>cur_bound[3]){
#           cur_val = c(cur_val,effect[k])
#         }
#         
#       }
#       if (length(cur_val) == 0){cur_val = 0}else{cur_val = mean(cur_val)}
#       
#       z = c(z,cur_val)
#     }
#   }
#   
#   spatial_polygons = SpatialPolygons(polygon_list)
#   z = data.frame(ID = id, value = z)
#   
#   spatial_data_frame = SpatialPolygonsDataFrame(spatial_polygons,z,match.ID = "ID")
#   return (spatial_data_frame)
# 
# }


# Based on our analysis, we get a log hazard ratio estimate for each spatial unit, 
# this spatial unit can be identified by a turple (TOWN, RANGE, SECTION). So we need 
# to draw each spatial unit on the actual map and fill in the value. First, we down
# -load a shape file from 
# https://gis-michigan.opendata.arcgis.com/datasets/public-land-survey-sections?geometry=-86.152%2C44.568%2C-82.683%2C45.249,
# then load it into R. This shape file contains the spatial Polygon of each spatial unit.





# read the shape file, convert it into a spatial dataframe, add another column id 
# for the spatial data frame. Then match the spatial data frame with our estimated 
# spatial log hazard by matching id column.
tbdata_adj <- read_excel("../data/tbdata_adj.xlsx")
tbdata_adj$SECTION = as.character(tbdata_adj$SECTION)
tbdata_adj$TOWN = as.character(tbdata_adj$TOWN)
tbdata_adj$RANGE = as.character(tbdata_adj$RANGE)
tbdata_adj$SECTION[tbdata_adj$SECTION=="1"] = "01"
tbdata_adj$SECTION[tbdata_adj$SECTION=="2"] = "02"
tbdata_adj$SECTION[tbdata_adj$SECTION=="3"] = "03"
tbdata_adj$SECTION[tbdata_adj$SECTION=="4"] = "04"
tbdata_adj$SECTION[tbdata_adj$SECTION=="5"] = "05"
tbdata_adj$SECTION[tbdata_adj$SECTION=="6"] = "06"
tbdata_adj$SECTION[tbdata_adj$SECTION=="7"] = "07"
tbdata_adj$SECTION[tbdata_adj$SECTION=="8"] = "08"
tbdata_adj$SECTION[tbdata_adj$SECTION=="9"] = "09"
tbdata_adj = tbdata_adj %>% select(SECTION, RANGE,TOWN, Adj_ID)
data = tbdata_adj
shapefile =readShapeSpatial("../data/Public_Land_Survey_Sections-shp/Public_Land_Survey_Sections.shp",delete_null_obj=TRUE)
shapefile$TOWN = as.character(shapefile$TOWN)
shapefile$RANGE = as.character(shapefile$RANGE)
shapefile$SECTION = as.character(shapefile$SECTION)
data_cleaned = shapefile[shapefile$TOWN %in% unique(data$TOWN),]
data_cleaned = data_cleaned[data_cleaned$RANGE %in% unique(data$RANGE),]
data_cleaned = data_cleaned[data_cleaned$SECTION %in% unique(data$SECTION),]
data_cleaned$id = 0
keep_ind = c()
for (i in 1:nrow(data_cleaned)){
  for (j in 1:nrow(data)){
    if (data$TOWN[j] == data_cleaned$TOWN[i]){
      if (data$RANGE[j] == data_cleaned$RANGE[i]){
        if (data$SECTION[j] == data_cleaned$SECTION[i]){
          keep_ind = c(keep_ind, i)
          data_cleaned$id[i] = data$Adj_ID[j]
        }
      }
    }
  }
}
data_cleaned = data_cleaned[keep_ind,]
#which(duplicated(data_cleaned$id))[-1]
data_cleaned = data_cleaned[-which(duplicated(data_cleaned$id))[-1],]


# using model24 result analysis to get the spatial log hazard ratio mspace_effect
load("model24_for_map.RData")

# male log hazard ratio
z = data.frame(id = data$Adj_ID, value = mspace_effect)
temp = merge(data_cleaned,z)
#temp = create_spatial_data_frame(mspace_effect,location_data[,1],location_data[,2],20,20)
#temp = create_spatial_data_frame(fspace_effect,location_data[,1],location_data[,2],20,20)

# female log hazard ratio
#z = data.frame(id = data$Adj_ID, value = fspace_effect)
#temp = merge(data_cleaned,z)

# sd of male spatial term
#temp = create_spatial_data_frame(mvspace_sd[1,],location_data[,1],location_data[,2],20,20)
# z = data.frame(id = data$Adj_ID, value = mvspace_sd[1,])
# temp = merge(data_cleaned,z)

# sd of female spatial term
# temp = create_spatial_data_frame(mvspace_sd[2,],location_data[,1],location_data[,2],20,20)


# using EDA result analysis to get the spatial positive rate distribution
# load("eda_for_map.RData")
# z = data.frame(id = data$Adj_ID, value = space_positive_rate)
# temp = merge(data_cleaned,z)
#temp = create_spatial_data_frame(space_positive_rate,location_data[,1],location_data[,2],20,20)




# use ggmap to get a base map, then plot the spatial polygons on the basemap and fill
# in the values of either positive rate or spatial hazard ratio for each unit.

# first, to use ggmap function, you need to register a google api
register_google(key="AIzaSyBx4Key-P9lg0FUE5S0SyZwFSeQYQz0Szg")

# get michigan's map by longtitude and latitude
michigan_base_map <- get_map(location=c(lon = -83.8, lat = 44.85), zoom=10, maptype = 'terrain-background', source = 'stamen')


# transform it into sf object
temp = st_as_sf(temp)

# set its crs to michigan georeference system, that's because our original dataset uses this crs.
temp <- st_set_crs(temp,3591)

# convert continuous response variables to discrete intervals
num_of_intervals = 6
breaks_qt = classIntervals(temp$value,n=num_of_intervals,style = "quantile")
br = breaks_qt$brks
br = unique(br)
offs = 0.0000001
br[1] = br[1] - offs
br[length(br)] = br[length(br)] + offs
temp$value_bracket <- cut(temp$value, br)

# project the spatial dataframe on the real map. Notice that the real map uses the WGS84 crs
# which means that you should transform the crs to WGS84.
temp2 <- st_transform(temp,4326)

# main plot function to obtain the static plot
ggmap(michigan_base_map)+
  geom_sf(data = temp, aes(fill=value_bracket),inherit.aes = FALSE,alpha = 0.4) +
  scale_fill_brewer(palette="OrRd")+
  coord_sf(crs = st_crs(4326))+
  labs(fill="standard deviation of female log hazard ratio") # +  labs(fill="positive rate")


# another way to make the interactive plot
# tm_shape(temp2) +
#   tm_polygons("value",style = "quantile")
# 
# tmap_mode("view")
# tmap_last()

# main plot function to obtain the interactive plot using leaflet
pal_fun <- colorQuantile("YlOrRd", NULL, n = num_of_intervals)

p_popup <- paste0("<strong> Log hazard ratio: </strong>", temp2$value)

leaflet(temp2) %>%
  addPolygons(
    stroke = TRUE, # remove polygon borders
    opacity = 0.1,
    fillColor = ~pal_fun(value), # set fill color with function from above and value
    fillOpacity = 0.8, smoothFactor = 0.5,# make it nicer
    popup = p_popup) %>%
    addTiles() %>%
    addLegend("bottomright",
              colors = brewer.pal(num_of_intervals,"YlOrRd"),
              labels = levels(temp2$value_bracket),
              #title = 'prevalence rate')
              title = 'Spatial log hazard ratio for male') # legend title

# publish the leaflet dynamic plot   
# https://rpubs.com/kehui/686144  Spatial log hazard ratio for male
# https://rpubs.com/kehui/686143 Spatial log hazard ratio for female

# ancilary functions 
# boulder_df <- data.frame(lon =data$X,
#                          lat = data$Y
#                          )

# convert to spatial points
# coordinates(boulder_df) <- 1:2
# crs(boulder_df) <- CRS("+init=EPSG:3591")
#crs(boulder_df) <- CRS("+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0
#")
# boulder_df_geog <- spTransform(boulder_df,CRS("+proj=longlat +datum=WGS84"))
# 
# coordinates(boulder_df_geog)

# boulder_df_geog <- as.data.frame(coordinates(boulder_df_geog))
# boulder_df_geog$value = as.factor(sample(c(1,2,3),nrow(data),replace = TRUE))
# map <- get_map("michigan")
# map = get_map()
# MichiganMap <- ggmap(map)
# MichiganMap
# 
# ggplot()+ geom_point(aes(x = lon, y = lat, colour = value, size = value),
#                      data = boulder_df_geog)
# 
# MichiganMap +
#   geom_point(aes(x = lon, y = lat, colour = value, size = value),
#              data = boulder_df_geog)

# us<-getData('GADM', country='USA', level=2)  #Get the County Shapefile for the US
# 
# michigan = subset(us,NAME_1=="Michigan")
# MichiganMap + geom_sf(michigan)
# 
# test = st_as_sf(michigan)
# str(test)
# st_crs(test)
# 
# test$value = rnorm(nrow(test))
# 
# plot(test)
# MichiganMap+geom_sf(test,aes(fill = value))


# 
# ggplot(test)+geom_sf(aes(fill=value))
# 
# ggmap(michigan_base_map)+
#   geom_sf(data = test, aes(fill=value),inherit.aes = FALSE) + 
#   coord_sf(crs = st_crs(4326))
# 
# 
# ggmap(michigan_base_map)+geom_point(data = boulder_df_geog, aes(x=lon,y=lat,fill=value))
# 
  
  
