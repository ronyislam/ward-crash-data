install.packages('pacman')

library(pacman)

p_load('lubridate',
       'RSocrata',
       'dplyr',
       'leaflet',
       'leaflet.extras',
       'sp',
       'sf',
       'geosphere',
       'dismo',
       'rgeos',
       'mapview',
       'rgdal',
       'pbapply',
       'formattable',
       'Hmisc',
       'forecast')


chicagoCrashPeople <- read.socrata("https://data.cityofchicago.org/resource/u6pd-qa9d.json")
chicagoCrashCrashes <- read.socrata("https://data.cityofchicago.org/resource/85ca-t3if.json")

dateQueryString <- paste0("crash_date between '",ymd('2018-01-01'), "' and '", ymd('2023-11-13'),"'")

chicagoCrashPeople <- read.socrata(paste0("https://data.cityofchicago.org/resource/u6pd-qa9d.json?$where=",dateQueryString))
chicagoCrashData <- read.socrata(paste0("https://data.cityofchicago.org/resource/85ca-t3if.json?$where=",dateQueryString)) %>%
  mutate(latitude = as.numeric(latitude),
         longitude = as.numeric(longitude),
         injuries_total = as.integer(injuries_total),
         injuries_fatal = as.integer(injuries_fatal)
  )

wardMap <- read_sf("https://data.cityofchicago.org/resource/p293-wvbd.geojson")


##==============================================================================
## DOWNLOAD DATA
##==============================================================================

## Contained code to download Chicago's wards
shpCityWards <- local({
  cur <- getwd()
  on.exit(setwd(cur))
  
  tmp <- tempfile(fileext = ".zip")
  setwd(dirname(tmp))
  # url <- "https://data.cityofchicago.org/api/geospatial/sp34-6z76?method=export&format=Shapefile"
  url <- "https://data.cityofchicago.org/api/geospatial/sp34-6z76?method=export&format=GeoJSON"
  download.file(url, destfile = tmp)
  shp <- rgdal::readOGR(basename(tmp), stringsAsFactors = FALSE)
  shp
})

## Generate city outline
shpCityOutline <- rgeos::gUnaryUnion(as(shpCityWards, "SpatialPolygons"))



##==============================================================================
## LABEL WITH centroid FUNCTION
##==============================================================================
## Source for functions:
## https://gis.stackexchange.com/a/265475/78424

#' find the center of mass / furthest away from any boundary
#' 
#' Takes as input a spatial polygon
#' @param pol One or more polygons as input
#' @param ultimate optional Boolean, TRUE = find polygon furthest away from centroid. False = ordinary centroid

require(rgeos)
require(sp)

centroid <- function(pol,ultimate=TRUE,iterations=5,initial_width_step=10){
  if (ultimate){
    new_pol <- pol
    # For every polygon do this:
    for (i in 1:length(pol)){
      width <- -initial_width_step
      area <- gArea(pol[i,])
      centr <- pol[i,]
      wasNull <- FALSE
      for (j in 1:iterations){
        if (!wasNull){ # stop when buffer polygon was alread too small
          centr_new <- gBuffer(centr,width=width)
          # if the buffer has a negative size:
          substract_width <- width/20
          while (is.null(centr_new)){ #gradually decrease the buffer size until it has positive area
            width <- width-substract_width
            centr_new <- gBuffer(centr,width=width)
            wasNull <- TRUE
          }
          # if (!(is.null(centr_new))){
          #   plot(centr_new,add=T)
          # }
          new_area <- gArea(centr_new)
          #linear regression:
          slope <- (new_area-area)/width
          #aiming at quarter of the area for the new polygon
          width <- (area/4-area)/slope
          #preparing for next step:
          area <- new_area
          centr<- centr_new
        }
      }
      #take the biggest polygon in case of multiple polygons:
      d <- disaggregate(centr)
      if (length(d)>1){
        biggest_area <- gArea(d[1,])
        which_pol <- 1                             
        for (k in 2:length(d)){
          if (gArea(d[k,]) > biggest_area){
            biggest_area <- gArea(d[k,])
            which_pol <- k
          }
        }
        centr <- d[which_pol,]
      }
      #add to class polygons:
      new_pol@polygons[[i]] <- remove.holes(new_pol@polygons[[i]])
      new_pol@polygons[[i]]@Polygons[[1]]@coords <- centr@polygons[[1]]@Polygons[[1]]@coords
    }
    centroids <- gCentroid(new_pol,byid=TRUE)
  }else{
    centroids <- gCentroid(pol,byid=TRUE)  
  }  
  return(centroids)
}

#Given an object of class Polygons, returns
#a similar object with no holes

remove.holes <- function(Poly){
  # remove holes
  is.hole <- lapply(Poly@Polygons,function(P)P@hole)
  is.hole <- unlist(is.hole)
  polys <- Poly@Polygons[!is.hole]
  Poly <- Polygons(polys,ID=Poly@ID)
  # remove 'islands'
  max_area <- largest_area(Poly)
  is.sub <- lapply(Poly@Polygons,function(P)P@area<max_area)  
  is.sub <- unlist(is.sub)
  polys <- Poly@Polygons[!is.sub]
  Poly <- Polygons(polys,ID=Poly@ID)
  Poly
}
largest_area <- function(Poly){
  total_polygons <- length(Poly@Polygons)
  max_area <- 0
  for (i in 1:total_polygons){
    max_area <- max(max_area,Poly@Polygons[[i]]@area)
  }
  max_area
}


labs <- centroid(pol = shpCityWards, 
                 ultimate = TRUE,
                 iterations = 150,
                 initial_width_step = .01)
labs$ward <- shpCityWards$ward


###


alderList <- read.socrata("https://data.cityofchicago.org/resource/htai-wnw4.json")

write.csv2(chicagoCrashPeople, file = "data/chicago_crash_people.csv")
write.csv2(chicagoCrashCrashes, file = "data/chicago_crash_crashes.csv")

chicagoCrashPeople <- read.csv2(file = "data/chicago_crash_people.csv")
chicagoCrashCrashes <- read.csv2(file = "data/all_chicago_crashes.csv")

injuriesPerCrash <- chicagoCrashPeople %>%
  group_by(crash_record_id) %>%
  summarise(
    None = sum(injury_classification == 'NO INDICATION OF INJURY' | is.na(injury_classification)),
    Reported = sum(injury_classification == 'REPORTED, NOT EVIDENT'),
    Minor = sum(injury_classification == 'NONINCAPACITATING INJURY'),
    Severe = sum(injury_classification == 'INCAPACITATING INJURY'),
    Fatal = sum(injury_classification == 'FATAL'),
    Cost = None*7429.2+Reported*27190.8+Minor*44634+Severe*172022.4+Fatal*1927973
  ) %>%
  mutate(CostPretty = formattable::currency(Cost)) 

chicagoCrashCrashes <- chicagoCrashData %>%
  left_join(injuriesPerCrash)


chicagoCrashCrashes["Cost"][is.na(chicagoCrashCrashes["Cost"])] <- 0
chicagoCrashCrashes["injuries_total"][is.na(chicagoCrashCrashes["injuries_total"])] <- 0
chicagoCrashCrashes["injuries_fatal"][is.na(chicagoCrashCrashes["injuries_fatal"])] <- 0

chicagoCrash <- chicagoCrashPeople %>%
  left_join(chicagoCrashCrashes, by = c("crash_record_id", "crash_date"))





crashes <- chicagoCrash %>%
  filter(crash_date >= '2018-01-01', longitude != 0)
  #filter(injury_classification == "NONINCAPACITATING INJURY" | injury_classification ==  "REPORTED, NOT EVIDENT"
  #       | injury_classification == "INCAPACITATING INJURY"  | injury_classification == "FATAL" ,
  #       crash_date >= '2017-01-01', longitude != 0)

chicagoCrash %>%
  filter(injury_classification=='FATAL') %>%
  group_by(year(crash_date)) %>%
  summarise(n())

crashesLongLats <- crashes %>%
  dplyr::select(longitude, latitude, crash_record_id, person_id)

# create a points collection
pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(crashesLongLats), 
                                     function(i) {st_point(as.numeric(crashesLongLats[i, 1:2]))}), list("crs" = 4326))) 

pnts_trans <- st_transform(pnts_sf, 2163) # apply transformation to pnts sf
tt1_trans <- st_transform(wardMap, 2163)      # apply transformation to polygons sf

# intersect and extract state name
crashesLongLats$ward <- as.integer(paste(apply(st_intersects(tt1_trans, pnts_trans, sparse = FALSE), 2, 
                     function(col) { 
                       tt1_trans[which(col), ]$ward
                     })))


crashesWard <- crashes %>%
  inner_join(crashesLongLats %>% filter(!is.na(ward))) %>%
  mutate(ward = as.factor(ward))


View(crashesWard %>% group_by(person_type,is.na(age)) %>% summarise(count=n()))

write.csv(crashesWard%>%st_drop_geometry()%>%dplyr::select(-location.coordinates), file="peopleCrashesWard.csv")

crashesWardSum <- crashesWard %>%
#  filter(person_type != 'DRIVER' & person_type != 'PASSENGER') %>%
  count(ward) %>%
  rename(totalCrashes = n) %>%
  mutate(ward = as.integer(paste(ward)))


crashesInjuriesWard <- crashesWard %>%
  filter(injury_classification == "NONINCAPACITATING INJURY" | injury_classification ==  "REPORTED, NOT EVIDENT"
               | injury_classification == "INCAPACITATING INJURY"  | injury_classification == "FATAL" ,
               crash_date >= '2018-01-01', longitude != 0)

crashesByMonth <- crashesInjuriesWard %>%
#  filter(person_type != 'DRIVER' & person_type != 'PASSENGER') %>%
  mutate(year_month = as_date(paste(year(crash_date),month(crash_date),'01', sep ="-"))) %>%
  mutate(month = as.integer(month(crash_date))) %>%
  mutate(year = as.factor(year(crash_date))) %>%
  filter(year_month < rollback(today()) + 1) %>%
  count(year_month, year, month)

plot<-ggplot(crashesByMonth, aes(x = year_month, y =n, color = year)) +
  geom_line() +
  geom_point()

plot

crashesByMonthWard <- crashesInjuriesWard %>%
#  filter(person_type != 'DRIVER' & person_type != 'PASSENGER') %>%
  mutate(year_month = as_date(paste(year(crash_date),month(crash_date),'01', sep ="-"))) %>%
  mutate(month = as.integer(month(crash_date))) %>%
  mutate(year = as.factor(year(crash_date))) %>%
  count(year_month, year, month, ward)

crashesByMonthInAWard <- crashesByMonthWard %>% 
  mutate(ward = as.factor(ward)) %>%
  filter(ward == 44)

ts_data <- ts(crashesByMonthInAWard$n, start = c(2018,11), 
              end = c(2023,11),frequency = 12) 

decomposedRes <- decompose(ts_data, type="mult") # use type = "additive" for additive components
plot(decomposedRes)
stlRes <- stl(ts_data, s.window = "periodic")

postPandemic <- window(stlRes$time.series,start = c(2020,4), 
                       end = c(2023,11),frequency = 12)

plot(postPandemic)

for(i in 1:50){
  crashesByMonthInAWard <- crashesByMonthWard %>% 
    filter(ward == i)
  
  ts_data <- ts(crashesByMonthInAWard$n, start = c(2018,1), 
                end = c(2022,11),frequency = 12) 
  
  decomposedRes <- decompose(ts_data, type="mult") # use type = "additive" for additive components
  stlRes <- stl(ts_data, s.window = "periodic")
  
  postPandemic <- window(stlRes$time.series,start = c(2020,4), 
                         end = c(2022,5),frequency = 12)
  
  
  toDate <- function(tt) {
    year <- as.integer(time(tt))
    month <- as.integer(cycle(tt))  # first week of year is 1, etc.
    as_date(paste(year,month,'01', sep ="-"))
  }
  
  ts_df <- data.frame(dates = toDate(postPandemic), series = c(postPandemic[,2]))
  model <- lm(series~dates, ts_df)
  wardPedCyclistInjuryFatalities[i,5] <- model$coefficients[2]
}
names(wardPedCyclistInjuryFatalities)[5] = "postPandemicTrend"

crashesByYearWard <- crashesInjuriesWard %>%
#  filter(person_type != 'DRIVER' & person_type != 'PASSENGER') %>%
  mutate(year_month = as_date(paste(year(crash_date),month(crash_date),'01', sep ="-"))) %>%
  mutate(month = as.integer(month(crash_date))) %>%
  mutate(year = as.integer(year(crash_date))) %>%
  count(year, ward) %>%
  rename("reported_injuries" = n)

ggplot(crashesByYearWard %>%   mutate(ward = as.factor(ward)) %>% filter(ward == 44), aes(x = year, y = reported_injuries, color = ward)) +
  geom_line() +
  geom_point() 

crashesByYearWardYTD <- crashesInjuriesWard %>%
#  filter(person_type != 'DRIVER' & person_type != 'PASSENGER') %>%
  mutate(year_month = as_date(paste(year(crash_date),month(crash_date),'01', sep ="-"))) %>%
  mutate(month = as.integer(month(crash_date))) %>%
  mutate(year = as.integer(year(crash_date))) %>%
  filter(month < month(today())) %>%
  count(year, ward) %>%
  rename("reported_injuries_ytd" = n)

ggplot(crashesByYearWardYTD %>%   mutate(ward = as.factor(ward)) %>% filter(ward == 44), aes(x = year, y = reported_injuries_ytd, color = ward)) +
  geom_line() +
  geom_point() 

crashesByYearWardYTD %>% filter(year==2019)%>%summarise(sum(reported_injuries_ytd))
#+ 
#  geom_smooth(method = "lm", col = "red")


crashesWard$person_type_factor <- factor(crashesWard$person_type)

factpal <- colorFactor(palette='Dark2', crashesWard$person_type_factor)


map<-leaflet(crashesWard %>%filter(ward==44)) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircleMarkers( color = ~factpal(person_type_factor), stroke=FALSE, fillOpacity =1, radius=2) %>%
  addLegend("topright", pal = factpal, values = ~person_type_factor,
                                                                 title = "Person Type",
                                                                 opacity = 1
  )

map

map<-leaflet(crashesWard %>%filter(ward==44, person_type=='BICYCLE')) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircleMarkers( color = ~factpal(person_type_factor), stroke=FALSE, fillOpacity =1, radius=2) %>%
  addLegend("topright", pal = factpal, values = ~person_type_factor,
            title = "Person Type",
            opacity = 1
  )

map

map<-leaflet(crashesWard %>%filter(ward==44, person_type=='PEDESTRIAN')) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircleMarkers( color = ~factpal(person_type_factor), stroke=FALSE, fillOpacity =1, radius=2) %>%
  addLegend("topright", pal = factpal, values = ~person_type_factor,
            title = "Person Type",
            opacity = 1
  )

map


crashesWard$injury_classification_factor <- factor(crashesWard$injury_classification)

injuryPal <- colorFactor(palette='Dark2', crashesWard$injury_classification_factor)

map<-leaflet(crashesWard %>%filter(ward==44, injury_classification!='NO INDICATION OF INJURY')) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircleMarkers( color = ~injuryPal(injury_classification_factor), stroke=FALSE, fillOpacity =1, radius=2) %>%
  addLegend("topright", pal = injuryPal, values = ~injury_classification_factor,
            title = "Injury Type",
            opacity = 1
  )

map

map<-leaflet(crashesWard %>%filter(ward==44, injury_classification=='FATAL')) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircleMarkers( color = ~injuryPal(injury_classification_factor), stroke=FALSE, fillOpacity =1, radius=2) %>%
  addLegend("topright", pal = injuryPal, values = ~injury_classification_factor,
            title = "Injury Type",
            opacity = 1
  )

map

map<-leaflet(crashesWard %>%filter(ward==44, injury_classification=='INCAPACITATING INJURY')) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircleMarkers( color = ~injuryPal(injury_classification_factor), stroke=FALSE, fillOpacity =1, radius=2) %>%
  addLegend("topright", pal = injuryPal, values = ~injury_classification_factor,
            title = "Injury Type",
            opacity = 1
  )

map

##ALL CRASHES CLUSTER
selectWardCrashes <- crashesWard %>% filter(ward %in% c(44,46,47,43))

de = MASS::kde2d(selectWardCrashes$longitude,selectWardCrashes$latitude, n=75)
image(de)
cl = contourLines(de, nlevels=10)
cllinesPolyLines = do.call(rbind,
                           Map(function(x){
                             st_as_sf(
                               data.frame(
                                 Z=x$level,
                                 geometry=st_sfc(
                                   st_linestring(cbind(x$x, x$y))
                                 )
                               )
                             )},cl
                           ))

cllines = do.call(rbind,
                  Map(function(x){
                    st_as_sf(
                      data.frame(
                        Z=x$level,
                        geometry=st_sfc(
                          st_linestring(cbind(x$x, x$y))
                        )
                      )
                    )},cl
                  ))

plot(cllinesPolyLines)

pal2 <- colorNumeric(
  palette = "Dark2",
  cllines$Z)

leaflet(cllinesPolyLines) %>%
  addCircleMarkers(data=selectWardCrashes, radius=.5) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolylines(weight=2, color = ~pal2(Z))

map<-leaflet(crashesWard %>%filter(ward==44)) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircleMarkers( color = ~factpal(person_type_factor), stroke=FALSE, fillOpacity =.4, radius=2) %>%
  addLegend("topright", pal = factpal, values = ~person_type_factor,
            title = "Person Type",
            opacity = 1
  ) %>%
  addPolylines(data=cllinesPolyLines, weight=2, color = ~pal2(Z)) %>%
  addPolygons(data = shpCityWards, fillOpacity = 0, 
              weight = 3, label = ~paste("Ward:", ward) )
map

###BIKE PED CRASH CLUSTER
selectWardCrashes <- crashesWard %>% filter(ward %in% c(44,46,47,43)) %>% filter(person_type %in% c('BICYCLE', 'PEDESTRIAN'))

de = MASS::kde2d(selectWardCrashes$longitude,selectWardCrashes$latitude, n=50)
image(de)
cl = contourLines(de, nlevels=32)
cllinesPolyLines = do.call(rbind,
                           Map(function(x){
                             st_as_sf(
                               data.frame(
                                 Z=x$level,
                                 geometry=st_sfc(
                                   st_linestring(cbind(x$x, x$y))
                                 )
                               )
                             )},cl
                           ))

cllines = do.call(rbind,
                  Map(function(x){
                    st_as_sf(
                      data.frame(
                        Z=x$level,
                        geometry=st_sfc(
                          st_linestring(cbind(x$x, x$y))
                        )
                      )
                    )},cl
                  ))

plot(cllinesPolyLines)

pal2 <- colorNumeric(
  palette = "Dark2",
  cllines$Z)

map<-leaflet(crashesWard %>%filter(ward==44) %>% filter(person_type %in% c('BICYCLE', 'PEDESTRIAN'))) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircleMarkers( color = ~factpal(person_type_factor), stroke=FALSE, fillOpacity =.4, radius=2) %>%
  addLegend("topright", pal = factpal, values = ~person_type_factor,
            title = "Person Type",
            opacity = 1
  ) %>%
  addPolylines(data=cllinesPolyLines, weight=2, color = ~pal2(Z)) %>%
  addPolygons(data = shpCityWards, fillOpacity = 0, 
              weight = 3, label = ~paste("Ward:", ward) )
map

###Injuries
selectWardCrashes <- crashesWard %>% filter(ward %in% c(44,46,47,43)) %>% filter(injury_classification!='NO INDICATION OF INJURY')

de = MASS::kde2d(selectWardCrashes$longitude,selectWardCrashes$latitude, n=50)
image(de)
cl = contourLines(de, nlevels=32)
cllinesPolyLines = do.call(rbind,
                           Map(function(x){
                             st_as_sf(
                               data.frame(
                                 Z=x$level,
                                 geometry=st_sfc(
                                   st_linestring(cbind(x$x, x$y))
                                 )
                               )
                             )},cl
                           ))

cllines = do.call(rbind,
                  Map(function(x){
                    st_as_sf(
                      data.frame(
                        Z=x$level,
                        geometry=st_sfc(
                          st_linestring(cbind(x$x, x$y))
                        )
                      )
                    )},cl
                  ))

plot(cllinesPolyLines)

pal2 <- colorNumeric(
  palette = "Dark2",
  cllines$Z)

map<-leaflet(crashesWard %>%filter(ward==44, injury_classification!='NO INDICATION OF INJURY')) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircleMarkers( color = ~injuryPal(injury_classification_factor), stroke=FALSE, fillOpacity =.4, radius=2) %>%
  addLegend("topright", pal = injuryPal, values = ~injury_classification_factor,
            title = "Injury Type",
            opacity = 1
  )%>%
  addPolylines(data=cllinesPolyLines, weight=2, color = ~pal2(Z)) %>%
  addPolygons(data = shpCityWards, fillOpacity = 0, 
              weight = 3, label = ~paste("Ward:", ward) )
map


###all intersections

intersections <- read.csv("../quick-build-map/all_intersections.csv") %>% rename(count=longitude,  longitude=latitude, latitude=X)



intersectionLongLats <- intersections %>%
  dplyr::select(longitude, latitude)

# create a points collection
pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(intersectionLongLats), 
                                     function(i) {st_point(as.numeric(intersectionLongLats[i, 1:2]))}), list("crs" = 4326))) 

pnts_trans <- st_transform(pnts_sf, 2163) # apply transformation to pnts sf
tt1_trans <- st_transform(wardMap, 2163)      # apply transformation to polygons sf

# intersect and extract state name
intersectionLongLats$ward <- as.integer(paste(apply(st_intersects(tt1_trans, pnts_trans, sparse = FALSE), 2, 
                                               function(col) { 
                                                 tt1_trans[which(col), ]$ward
                                               })))




# Create buffer around each intersection

intersectionLongLats <- intersectionLongLats%>%
  filter(ward==44) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove=FALSE)

buffer_size <- 34 # meters
intersection_buffer <- as(st_buffer(intersectionLongLats, buffer_size), "sf")


chicagoCrashCrashes.sf <- chicagoCrashData %>%
  filter(year(crash_date)>=2018, latitude!=0) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove=FALSE)

chicagoCrashPeople.sf <- crashesWard %>%
  filter(ward %in% c(44,46,47,43)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove=FALSE)

intersectionLongLats$crash_count <- lengths(st_intersects(intersection_buffer, chicagoCrashCrashes.sf))

intersectionLongLats$injury_count <- lengths(st_intersects(intersection_buffer, chicagoCrashPeople.sf %>% filter(injury_classification!='NO INDICATION OF INJURY')))

intersectionLongLats$crash_bike_ped_count <- lengths(st_intersects(intersection_buffer, chicagoCrashPeople.sf %>% filter(person_type=='CYCLIST' | person_type=='PEDESTRIAN')))

intersectionCrashPal<- colorBin(palette='Dark2', domain=intersectionLongLats$crash_count, bins=10, pretty=TRUE)
intersectionInjuryPal<- colorBin(palette='Dark2', domain=intersectionLongLats$injury_count, bins=5, pretty=TRUE)
intersectionBikePedCrashPal<- colorBin(palette='Dark2', domain=intersectionLongLats$crash_bike_ped_count, bins=3, pretty=TRUE)

leaflet(intersectionLongLats) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircles(color = ~intersectionCrashPal(crash_count), stroke=FALSE, fillOpacity =.4, radius=34) %>%
  addLegend("topright", pal = intersectionCrashPal, values = ~crash_count,
            title = "Crash Count",
            opacity = 1
  )

leaflet(intersectionLongLats) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircles(color = ~intersectionInjuryPal(injury_count), stroke=FALSE, fillOpacity =.4, radius=34) %>%
  addLegend("topright", pal = intersectionInjuryPal, values = ~injury_count,
            title = "Injury Count",
            opacity = 1
  )

leaflet(intersectionLongLats) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addCircles(color = ~intersectionBikePedCrashPal(crash_bike_ped_count), stroke=FALSE, fillOpacity =.4, radius=34) %>%
  addLegend("topright", pal = intersectionBikePedCrashPal, values = ~crash_bike_ped_count,
            title = "Bike/Ped Involved Crash Count",
            opacity = 1
  )


crashInjuriesPedCyclist <- crashesInjuriesWard %>%
  filter(injury_classification != "FATAL") %>% 
  filter(person_type != 'DRIVER' & person_type != 'PASSENGER') %>%
  mutate(year_month = as_date(paste(year(crash_date),month(crash_date),'01', sep ="-"))) %>%
  mutate(month = as.integer(month(crash_date))) %>%
  mutate(year = as.factor(year(crash_date)))

crashFatalitiesPedCyclist <- crashesInjuriesWard %>%
  filter(injury_classification == "FATAL") %>% 
  filter(person_type != 'DRIVER' & person_type != 'PASSENGER') %>%
  mutate(year_month = as_date(paste(year(crash_date),month(crash_date),'01', sep ="-"))) %>%
  mutate(month = as.integer(month(crash_date))) %>%
  mutate(year = as.factor(year(crash_date)))

wardPedCyclistInjuryFatalities <- crashInjuriesPedCyclist %>%
  count(ward) %>%
  rename(injuries = n) %>%
  left_join(crashFatalitiesPedCyclist %>%
              count(ward) %>%
              rename(fatalities = n)) %>%
  arrange(ward)

wardPedCyclistInjuryFatalities[is.na(wardPedCyclistInjuryFatalities)] <- 0

wardPedCyclistInjuryFatalities <- wardPedCyclistInjuryFatalities %>% mutate(ward=as.integer(ward))

wardPedCyclistInjuryFatalities <- wardPedCyclistInjuryFatalities %>%
  left_join(crashesWardSum) %>%
  mutate(crashes = totalCrashes-injuries-fatalities)

wardPedCyclistInjuryFatalities$injuryRate <- wardPedCyclistInjuryFatalities$injuries/wardPedCyclistInjuryFatalities$totalCrashes
wardPedCyclistInjuryFatalities$fatalityRate <- wardPedCyclistInjuryFatalities$fatalities/wardPedCyclistInjuryFatalities$totalCrashes

meanInjury <- mean(wardPedCyclistInjuryFatalities$injuries)
sdInjury <- sd(wardPedCyclistInjuryFatalities$injuries)

meanFatality <- mean(wardPedCyclistInjuryFatalities$fatalities)
sdFatality <- sd(wardPedCyclistInjuryFatalities$fatalities)

wardPedCyclistInjuryFatalities$injuryDeviations <- wardPedCyclistInjuryFatalities$injuries/sdInjury
wardPedCyclistInjuryFatalities$fatalityDeviations <- wardPedCyclistInjuryFatalities$fatalities/sdFatality
wardPedCyclistInjuryFatalities$injuryMeanDistance <- wardPedCyclistInjuryFatalities$injuries/meanInjury
wardPedCyclistInjuryFatalities$fatalityMeanDistance <- wardPedCyclistInjuryFatalities$fatalities/meanFatality

wardPedCyclistInjuryFatalities$crashRank[order(wardPedCyclistInjuryFatalities$totalCrashes)] <- 51-1:nrow(wardPedCyclistInjuryFatalities)
wardPedCyclistInjuryFatalities$injuryRank[order(wardPedCyclistInjuryFatalities$injuries)] <- 51-1:nrow(wardPedCyclistInjuryFatalities)
wardPedCyclistInjuryFatalities$fatalityRank[order(wardPedCyclistInjuryFatalities$fatalities)] <- 51-1:nrow(wardPedCyclistInjuryFatalities)


injuryRanking <- ggplot(data=wardPedCyclistInjuryFatalities, aes(x=injuryRank, y=injuries)) +
  geom_bar(stat="identity", width=.9, position = position_dodge(width=0.2), fill = "grey69") +
  geom_text(aes(label=ward), vjust=2, size=3.0, color = "black") +
  scale_x_reverse()
injuryRanking

fatalityRanking <- ggplot(data=wardPedCyclistInjuryFatalities, aes(x=fatalityRank, y=fatalities)) +
  geom_bar(stat="identity", width=.9, position = position_dodge(width=0.2), fill = "grey69") +
  geom_text(aes(label=ward), vjust=2, size=3.0, color = "black") +
  scale_x_reverse()
fatalityRanking


stackedBarData <- union_all(wardPedCyclistInjuryFatalities %>% 
  dplyr::select(ward, injuries),
  wardPedCyclistInjuryFatalities %>% 
    dplyr::select(ward, crashes)) %>%
  union_all(wardPedCyclistInjuryFatalities %>% 
              dplyr::select(ward, fatalities))%>%
  mutate(category = ifelse(!is.na(injuries),"injury",ifelse(!is.na(fatalities),"fatality","no injury"))) %>%
  left_join(wardPedCyclistInjuryFatalities, by = c("ward"))


ggplot(stackedBarData, aes(fill=category, y=crashes.y, x=fatalityRank)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label=ward), vjust=-0.3, size=3.5) +
  scale_x_reverse()


##MENU MONEY
menuMoney <- list()
for(i in 1:50){
  menuMoney[[i]]<- as_tibble(read_xlsx("data/2019-2020.xlsx", range = "A2:D600",sheet =i, skip = 1))  %>%
    filter(if_any(everything(), ~ !is.na(.))) %>%
    replace(3:3, as.numeric(.[[3]])) %>%
    rename(category = 4, cost = 3) %>%
    mutate(category = if_else(category == 'DS', 'Driver Support', if_else(category == 'PS','Ped Support', if_else(category == 'CE', 'Community Enhance', "??")))) %>%
    mutate(category = as_factor(category))
    
    
    menuMoney[[i]] <- menuMoney[[i]] %>%
      mutate(year = if_else(
      row_number()<which(menuMoney[[i]] == 'WARD 2019 BALANCE', arr.ind=TRUE)[1], 
      '2019', 
      '2020')
      ) %>%
    filter(!is.na(category))
}

ggplot(menuMoney[[1]], aes(fill=category, y=cost, x=year)) + 
  geom_bar(position="stack", stat="identity") + 
  ggtitle('Ward 1 Menu Spending')


injuryCountsByYear <- crashesWard %>%
  filter(month(crash_date)<=6) %>%
  count(year(crash_date), injury_classification) %>% 
  mutate(injury_class = if_else(injury_classification == 'FATAL', 4, if_else(injury_classification == 'INCAPACITATING INJURY', 3,if_else(injury_classification == 'NONINCAPACITATING INJURY', 2,if_else(injury_classification == 'REPORTED, NOT EVIDENT', 1,0))))) %>%
  mutate(injury_classification = if_else(injury_classification == 'FATAL', "Fatal", if_else(injury_classification == 'INCAPACITATING INJURY', "Severe Injury",if_else(injury_classification == 'NONINCAPACITATING INJURY', "Minor Injury",if_else(injury_classification == 'REPORTED, NOT EVIDENT', "Reported Injury","No Injury"))))) %>%
  rename(year = 'year(crash_date)', victims = n)


crashesWard %>% 
  count(year(crash_date))

fatalityVisionZero<-ggplot(injuryCountsByYear %>% filter(injury_class==4), aes(x = year, y =victims, color = injury_classification)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("black")) +
  geom_hline(aes(yintercept=130/2, linetype = "Zero Vision"),  colour = "#e41e2d") + 
  geom_segment(aes(x = 2018, xend = 2026, y = 130/2 , yend = 0, linetype = "Vision Zero"), colour = "#578bc6") +
  scale_linetype_manual(name = "Projections", values = c(2,2), 
                        guide = guide_legend(override.aes = list(color = c("#578bc6", "#e41e2d"))))+
  ylim(0, 170/2) + 
  xlim(2018, 2026) + 
  ggtitle("Vision Zero Chicago First Half 2023 Report \n Overall Fatalities") +
  labs(caption = "ronythebikeczar.substack.com") +
  xlab("Year") + 
  ylab("Victims") +
  labs(color = "Injury Type") +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold", family = "Cubano"),
        axis.title.x = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.title.y = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.text = element_text(color="black", family = "Cubano"),
        legend.title = element_text(color="black", family = "Cubano"),
        legend.text = element_text(color="black", family = "Cubano"))

fatalityVisionZero

severeInjuryVisionZero<-ggplot(injuryCountsByYear %>% filter(injury_class==3), aes(x = year, y =victims, color = injury_classification), color=black) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("black")) +
  geom_hline(aes(yintercept=2481/2, linetype = "Zero Vision"),  colour = "#e41e2d") + 
  geom_segment(aes(x = 2018, xend = 2026, y = 2481/2 , yend = 0, linetype = "Vision Zero"), colour = "#578bc6") +
  scale_linetype_manual(name = "Projections", values = c(2,2), 
                        guide = guide_legend(override.aes = list(color = c("#578bc6", "#e41e2d"))))+
  ylim(0, 2500/2) + 
  xlim(2018, 2026) + 
  ggtitle("Vision Zero Chicago First Half 2023 Report \n Overall Severe Injuries") +
  labs(caption = "ronythebikeczar.substack.com") +
  xlab("Year") + 
  ylab("Victims") +
  labs(color = "Injury Type") +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold", family = "Cubano"),
        axis.title.x = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.title.y = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.text = element_text(color="black", family = "Cubano"),
        legend.title = element_text(color="black", family = "Cubano"),
        legend.text = element_text(color="black", family = "Cubano"))
  
severeInjuryVisionZero


##VISION ZERO FOR VULNERABLE ROAD USERS
injuryCountsByYearPedCyclist <- crashesWard %>%
  filter(month(crash_date)<=6) %>%
  filter(person_type=="PEDESTRIAN") %>%
  count(year(crash_date), injury_classification) %>% 
  mutate(injury_class = if_else(injury_classification == 'FATAL', 4, if_else(injury_classification == 'INCAPACITATING INJURY', 3,if_else(injury_classification == 'NONINCAPACITATING INJURY', 2,if_else(injury_classification == 'REPORTED, NOT EVIDENT', 1,0))))) %>%
  mutate(injury_classification = if_else(injury_classification == 'FATAL', "Fatal", if_else(injury_classification == 'INCAPACITATING INJURY', "Severe Injury",if_else(injury_classification == 'NONINCAPACITATING INJURY', "Minor Injury",if_else(injury_classification == 'REPORTED, NOT EVIDENT', "Reported Injury","No Injury"))))) %>%
  rename(year = 'year(crash_date)', victims = n)

crashesWard %>% 
  count(year(crash_date))

fatalityVisionZeroPedCyclist<-ggplot(injuryCountsByYearPedCyclist %>% filter(injury_class==4), aes(x = year, y =victims, color = injury_classification)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("black")) +
  geom_hline(aes(yintercept=39/2, linetype = "Zero Vision"),  colour = "#e41e2d") + 
  geom_segment(aes(x = 2018, xend = 2026, y = 39/2 , yend = 0, linetype = "Vision Zero"), colour = "#578bc6") +
  scale_linetype_manual(name = "Projections", values = c(2,2), 
                        guide = guide_legend(override.aes = list(color = c("#578bc6", "#e41e2d"))))+
  ylim(0, 60/2) + 
  xlim(2018, 2026) + 
  ggtitle("Vision Zero Chicago First Half 2023 Report \n Ped and Cyclist Fatalities") +
  labs(caption = "ronythebikeczar.substack.com") +
  xlab("Year") + 
  ylab("Victims") +
  labs(color = "Injury Type") +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold", family = "Cubano"),
        axis.title.x = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.title.y = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.text = element_text(color="black", family = "Cubano"),
        legend.title = element_text(color="black", family = "Cubano"),
        legend.text = element_text(color="black", family = "Cubano"))

fatalityVisionZeroPedCyclist

severeInjuryVisionZeroPedCyclist<-ggplot(injuryCountsByYearPedCyclist %>% filter(injury_class==3), aes(x = year, y =victims, color = injury_classification)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("black")) +
  geom_hline(aes(yintercept=573/2, linetype = "Zero Vision"),  colour = "#e41e2d") + 
  geom_segment(aes(x = 2018, xend = 2026, y = 573/2 , yend = 0, linetype = "Vision Zero"), colour = "#578bc6") +
  scale_linetype_manual(name = "Projections", values = c(2,2), 
                        guide = guide_legend(override.aes = list(color = c("#578bc6", "#e41e2d"))))+
  ylim(0, 750) + 
  xlim(2018, 2026) + 
  ggtitle("Vision Zero Chicago First Half 2023 Report \n Ped and Cyclist Severe Injuries") +
  labs(caption = "ronythebikeczar.substack.com") +
  xlab("Year") + 
  ylab("Victims") +
  labs(color = "Injury Type") +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold", family = "Cubano"),
        axis.title.x = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.title.y = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.text = element_text(color="black", family = "Cubano"),
        legend.title = element_text(color="black", family = "Cubano"),
        legend.text = element_text(color="black", family = "Cubano"))

severeInjuryVisionZeroPedCyclist


crashesByYear<-ggplot(chicagoCrashCrashes %>% 
                                           filter(year(crash_date)>=2018, month(crash_date)<=6) %>% 
                                           count(year(crash_date)) %>%
                                           rename(year = 'year(crash_date)', crashes = n), aes(x = year, y =crashes)) +
  geom_line() +
  geom_point() +
  ggtitle("Vision Zero Chicago First Half 2023 Report \n Overall Traffic Crashes") +
  labs(caption = "ronythebikeczar.substack.com") +
  xlab("Year") + 
  ylab("Crashes") +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold", family = "Cubano"),
        axis.title.x = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.title.y = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.text = element_text(color="black", family = "Cubano"),
        legend.title = element_text(color="black", family = "Cubano"),
        legend.text = element_text(color="black", family = "Cubano"))

crashesByYear

crashesWithInjuryByYear<-ggplot(chicagoCrashCrashes %>% 
                        filter(year(crash_date)>=2018, month(crash_date)<=6, most_severe_injury != 'NO INDICATION OF INJURY') %>% 
                        count(year(crash_date)) %>%
                        rename(year = 'year(crash_date)', crashes = n), aes(x = year, y =crashes)) +
  geom_line() +
  geom_point() +
  ggtitle("Vision Zero Chicago First Half 2023 Report \n Traffic Crashes w/ Injury Involved") +
  labs(caption = "ronythebikeczar.substack.com") +
  xlab("Year") + 
  ylab("Crashes") +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold", family = "Cubano"),
        axis.title.x = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.title.y = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.text = element_text(color="black", family = "Cubano"),
        legend.title = element_text(color="black", family = "Cubano"),
        legend.text = element_text(color="black", family = "Cubano"))

crashesWithInjuryByYear

hitAndRunsFatalByYear<-ggplot(chicagoCrashCrashes %>% 
                        filter(year(crash_date)>=2018, month(crash_date)<=6, hit_and_run_i == 'Y', most_severe_injury=='FATAL') %>% 
                        count(year(crash_date)) %>%
                        rename(year = 'year(crash_date)', crashes = n), aes(x = year, y =crashes)) +
  geom_line() +
  geom_point() +
  ggtitle("Vision Zero Chicago First Half 2023 Report \n Fatal Hit and Runs") +
  labs(caption = "ronythebikeczar.substack.com") +
  xlab("Year") + 
  ylab("Crashes") +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold", family = "Cubano"),
        axis.title.x = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.title.y = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.text = element_text(color="black", family = "Cubano"),
        legend.title = element_text(color="black", family = "Cubano"),
        legend.text = element_text(color="black", family = "Cubano"))

hitAndRunsFatalByYear

hitAndRunsByYear<-ggplot(chicagoCrashCrashes %>% 
                           filter(year(crash_date)>=2018, month(crash_date)<=6, hit_and_run_i == 'Y') %>% 
                           count(year(crash_date)) %>%
                           rename(year = 'year(crash_date)', crashes = n), aes(x = year, y =crashes)) +
  geom_line() +
  geom_point() +
  ggtitle("Vision Zero Chicago First Half 2023 Report \n Overall Hit and Runs") +
  labs(caption = "ronythebikeczar.substack.com") +
  xlab("Year") + 
  ylab("Crashes") +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold", family = "Cubano"),
        axis.title.x = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.title.y = element_text(color="black", size=12, face="bold", family = "Cubano"),
        axis.text = element_text(color="black", family = "Cubano"),
        legend.title = element_text(color="black", family = "Cubano"),
        legend.text = element_text(color="black", family = "Cubano"))

hitAndRunsByYear

crashesWard %>%
  filter(ward==1,
         injury_classification=='FATAL') %>%
  count()

wardMap$ward = as.integer(wardMap$ward)

plot()
