library(RSocrata)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(sf)
library(Hmisc)
library(forecast)
library(readxl)
library(leaflet)
library(rgeos)
library(rgdal)
library(viridis)

chicagoCrashPeople <- read.socrata("https://data.cityofchicago.org/resource/u6pd-qa9d.json")
chicagoCrashCrashes <- read.socrata("https://data.cityofchicago.org/resource/85ca-t3if.json")
wardMap <- read_sf("https://data.cityofchicago.org/resource/p293-wvbd.geojson")
alderList <- read.socrata("https://data.cityofchicago.org/resource/htai-wnw4.json")

write.csv2(chicagoCrashPeople, file = "data/chicago_crash_people.csv")
write.csv2(chicagoCrashCrashes, file = "data/chicago_crash_crashes.csv")

chicagoCrashPeople <- read.csv2(file = "data/chicago_crash_people.csv")
chicagoCrashCrashes <- read.csv2(file = "data/all_chicago_crashes.csv")


chicagoCrash <- chicagoCrashPeople %>%
  left_join(chicagoCrashCrashes, by = c("crash_record_id", "crash_date"))

crashes <- chicagoCrash %>%
  filter(crash_date >= '2018-01-01', longitude != 0)
#filter(injury_classification == "NONINCAPACITATING INJURY" | injury_classification ==  "REPORTED, NOT EVIDENT"
#       | injury_classification == "INCAPACITATING INJURY"  | injury_classification == "FATAL" ,
#       crash_date >= '2017-01-01', longitude != 0)

chicagoCrash %>%
  filter(month(crash_date)>6, injury_classification=='FATAL') %>%
  group_by(year(crash_date)) %>%
  summarise(n())

crashesLongLats <- crashes %>%
  select(longitude, latitude, crash_record_id, person_id)

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


#########
#WARD LABELS
#########


## Contained code to download Chicago's wards
shpCityWards <- local({
  cur <- getwd()
  on.exit(setwd(cur))
  
  tmp <- tempfile(fileext = ".zip")
  setwd(dirname(tmp))
  url <- "https://data.cityofchicago.org/api/geospatial/p293-wvbd?method=export&format=GeoJSON"
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


###WARD LABELS END

crashesWard <- crashes %>%
  inner_join(crashesLongLats %>% filter(!is.na(ward))) %>%
  mutate(ward = as.factor(ward))

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
  filter(year_month < '2020-12-01') %>%
  count(year_month, year, month, ward)

crashesByMonthInAWard <- crashesByMonthWard %>% 
  mutate(ward = as.factor(ward)) %>%
  filter(ward == 4)

ts_data <- ts(crashesByMonthInAWard$n, start = c(2018,1), 
              end = c(2022,12),frequency = 12) 

decomposedRes <- decompose(ts_data, type="mult") # use type = "additive" for additive components
plot(decomposedRes)
stlRes <- stl(ts_data, s.window = "periodic")

postPandemic <- window(stlRes$time.series,start = c(2020,4), 
                       end = c(2022,11),frequency = 12)

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
  filter(year_month < '2023-01-01') %>%
  count(year, ward) %>%
  rename("reported_injuries" = n)

ggplot(crashesByYearWard %>%   mutate(ward = as.factor(ward)) %>% filter(ward == 4), aes(x = year, y = reported_injuries, color = ward)) +
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

ggplot(crashesByYearWardYTD %>%   mutate(ward = as.factor(ward)) %>% filter(ward == 4), aes(x = year, y = reported_injuries_ytd, color = ward)) +
  geom_line() +
  geom_point() 

crashesByYearWardYTD %>% filter(year==2019)%>%summarise(sum(reported_injuries_ytd))
#+ 
#  geom_smooth(method = "lm", col = "red")



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
  geom_bar(stat="identity", width=.9, position = position_dodge(width=0.2), fill = "black69") +
  geom_text(aes(label=ward), vjust=2, size=3.0, color = "black") +
  scale_x_reverse()
injuryRanking

fatalityRanking <- ggplot(data=wardPedCyclistInjuryFatalities, aes(x=fatalityRank, y=fatalities)) +
  geom_bar(stat="identity", width=.9, position = position_dodge(width=0.2), fill = "black69") +
  geom_text(aes(label=ward), vjust=2, size=3.0, color = "black") +
  scale_x_reverse()
fatalityRanking


stackedBarData <- union_all(wardPedCyclistInjuryFatalities %>% 
                              select(ward, injuries),
                            wardPedCyclistInjuryFatalities %>% 
                              select(ward, crashes)) %>%
  union_all(wardPedCyclistInjuryFatalities %>% 
              select(ward, fatalities))%>%
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

fatalitiesByWard <- data.frame(
ward=c(1:50), fatalities=NA) %>%
  left_join(crashesWard %>%
              mutate(ward = as.integer(ward) ,
                     injury_class = if_else(injury_classification == 'FATAL', 4, if_else(injury_classification == 'INCAPACITATING INJURY', 3,if_else(injury_classification == 'NONINCAPACITATING INJURY', 2,if_else(injury_classification == 'REPORTED, NOT EVIDENT', 1,0))))) %>%
              filter(month(crash_date)<=6,
                     year(crash_date)==2023,
                     injury_classification=='FATAL') %>%
              group_by(ward) %>%
              count() %>%
              rename("fatalities"="n"),
             by="ward"
            ) %>%
  mutate(fatalities=coalesce(fatalities.x, fatalities.y)) %>%
  select(ward,fatalities)


wardMap$fatalities <- fatalitiesByWard$fatalities

injuriesByWard <- data.frame(
  ward=c(1:50), injuries=NA) %>%
  left_join(crashesWard %>%
              mutate(ward = as.integer(ward) ,
                     injury_class = if_else(injury_classification == 'FATAL', 4, if_else(injury_classification == 'INCAPACITATING INJURY', 3,if_else(injury_classification == 'NONINCAPACITATING INJURY', 2,if_else(injury_classification == 'REPORTED, NOT EVIDENT', 1,0))))) %>%
              filter(month(crash_date)<=6,
                     year(crash_date)==2023,
                     injury_class>0) %>%
              group_by(ward) %>%
              count() %>%
              rename("injuries"="n"),
            by="ward"
  ) %>%
  mutate(injuries=coalesce(injuries.x, injuries.y)) %>%
  select(ward,injuries)


wardMap$injuries <- injuriesByWard$injuries

severeInjuriesByWard <- data.frame(
  ward=c(1:50), severeInjuries=NA) %>%
  left_join(crashesWard %>%
              mutate(ward = as.integer(ward) ,
                     injury_class = if_else(injury_classification == 'FATAL', 4, if_else(injury_classification == 'INCAPACITATING INJURY', 3,if_else(injury_classification == 'NONINCAPACITATING INJURY', 2,if_else(injury_classification == 'REPORTED, NOT EVIDENT', 1,0))))) %>%
              filter(month(crash_date)<=6,
                     year(crash_date)==2023,
                     injury_class==3) %>%
              group_by(ward) %>%
              count() %>%
              rename("severeInjuries"="n"),
            by="ward"
  ) %>%
  mutate(severeInjuries=coalesce(severeInjuries.x, severeInjuries.y)) %>%
  select(ward,severeInjuries)


wardMap$severeInjuries <- severeInjuriesByWard$severeInjuries

crashesByWard <- data.frame(
  ward=c(1:50), crashes=NA) %>%
  left_join(crashesWard %>%
              select(ward,crash_record_id,crash_date) %>%
              distinct() %>%
              mutate(ward = as.integer(ward)) %>%
              filter(month(crash_date)<=6,
                     year(crash_date)==2023) %>%
              group_by(ward) %>%
              count() %>%
              rename("crashes"="n"),
            by="ward"
  ) %>%
  mutate(crashes=coalesce(crashes.x, crashes.y)) %>%
  select(ward,crashes)


wardMap$crashes <- crashesByWard$crashes

wardMap$injuriesPerCrash <- round(wardMap$injuries/wardMap$crashes,2)
wardMap$severeInjuriesPerCrash <- round(wardMap$severeInjuries/wardMap$crashes,3)
wardMap$fatalitiesPerCrash <- round(wardMap$fatalities/wardMap$crashes,5)

wardMap$crashesPerDay <- round(wardMap$crashes/180,2)
wardMap$injuriesPerDay <- round(wardMap$injuries/180,2)
wardMap$severeInjuriesPerDay <- round(wardMap$severeInjuries/180,3)
wardMap$fatalitiesPerDay <- round(wardMap$fatalities/180,5)

wardMap$percentCrashes <- wardMap$crashes/sum(wardMap$crashes)*100
wardMap$percentInjuries <- wardMap$injuries/sum(wardMap$injuries)*100
wardMap$percentSevereInjuries <- wardMap$severeInjuries/sum(wardMap$severeInjuries)*100
wardMap$percentFatalities <- wardMap$fatalities/sum(wardMap$fatalities, na.rm=TRUE)*100

wardMap$percentInjuriesOverCrashes <- wardMap$injuries/sum(wardMap$injuries)/wardMap$percentCrashes
wardMap$percentSevereInjuriesOverCrashes <- wardMap$severeInjuries/sum(wardMap$severeInjuries)/wardMap$percentCrashes
wardMap$percentFatalitiesOverCrashes <- wardMap$fatalities/sum(wardMap$fatalities, na.rm=TRUE)/wardMap$percentCrashes

map_attr <- "© <a href='http://ronythebikeczar.substack.com'>ronythebikeczar.substack.com - Rony I.</a>"

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$severeInjuries
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(severeInjuries)) %>%
  addLegend(pal = pal, values = ~severeInjuries, title="# of Severe Injuries") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$injuries
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(injuries)) %>%
  addLegend(pal = pal, values = ~injuries, title="# of Injuries") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$fatalities
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(fatalities)) %>%
  addLegend(pal = pal, values = ~fatalities, title="# of Fatalities") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$severeInjuries/180
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(severeInjuries/180)) %>%
  addLegend(pal = pal, values = ~severeInjuries/180, title="Severe Injuries Per Day") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$injuries/180
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(injuries/180)) %>%
  addLegend(pal = pal, values = ~injuries/180, title="Injuries Per Day") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$fatalities/180
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(fatalities/180)) %>%
  addLegend(pal = pal, values = ~fatalities/180, title="Fatalities Per Day") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$crashesByWard
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(crashes)) %>%
  addLegend(pal = pal, values = ~crashes, title="# of Crashes") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$injuriesPerCrash
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(injuriesPerCrash)) %>%
  addLegend(pal = pal, values = ~injuriesPerCrash, title="Injuries Per Crash") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$severeInjuriesPerCrash
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(severeInjuriesPerCrash)) %>%
  addLegend(pal = pal, values = ~severeInjuriesPerCrash, title="Severe Injuries Per Crash") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$fatalitiesPerCrash
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(fatalitiesPerCrash)) %>%
  addLegend(pal = pal, values = ~fatalitiesPerCrash, title="Fatalities Per Crash") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$percentInjuries
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(percentInjuries)) %>%
  addLegend(pal = pal, values = ~percentInjuries, title="% Injuries") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$percentSevereInjuries
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(percentSevereInjuries)) %>%
  addLegend(pal = pal, values = ~percentSevereInjuries, title="% Severe Injuries") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$percentCrashes
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(percentCrashes)) %>%
  addLegend(pal = pal, values = ~percentCrashes, title="% Crashes") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$percentInjuries/wardMap$percentCrashes
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(percentInjuries/percentCrashes)) %>%
  addLegend(pal = pal, values = ~percentInjuries/percentCrashes, title="% Injuries v. % Crashes") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$percentSevereInjuries/wardMap$percentCrashes
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(percentSevereInjuries/percentCrashes)) %>%
  addLegend(pal = pal, values = ~percentSevereInjuries/percentCrashes, title="% Severe Injuries v. % Crashes") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 

pal <- colorNumeric(
  palette = "viridis",
  domain = wardMap$percentFatalities/wardMap$percentCrashes
)

leaflet(data=wardMap,options = leafletOptions(zoomControl = TRUE,
                                              zoomSnap = 0,
                                              zoomDelta = 0.1)) %>%
  addTiles(attribution = map_attr) %>%
  addProviderTiles("CartoDB.Voyager") %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(percentFatalities/percentCrashes)) %>%
  addLegend(pal = pal, values = ~percentFatalities/percentCrashes, title="% Fatalities v. % Crashes") %>%
  addLabelOnlyMarkers(data = labs, ~labs$x, ~labs$y, label = ~as.character(ward),
                      labelOptions = labelOptions(noHide = TRUE,
                                                  direction = "center",
                                                  offset = c(0, 0), opacity = 1, 
                                                  textsize = "12px", textOnly = TRUE, 
                                                  style = list("font-style" = "bold",
                                                               color = "black"))) 


View(as.data.frame(wardMap))
  
