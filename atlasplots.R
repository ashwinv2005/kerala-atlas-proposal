library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(sp)
library(sf)
library(rgeos)





################################# create dataset

## read eBird data, better to use more current info

rawpath = "ebd_IN_relSep-2019.txt"
sensitivepath = "Sensitive_India_may 2019.csv"

## file you sent

klpath = "Atlas.csv"

atlas = read.csv(klpath, header = F, sep = " ")
atlas = as.vector(atlas[,3])
atlas = atlas[atlas != "SNA"]

## read eBird data

preimp = c("CATEGORY","COMMON.NAME","SCIENTIFIC.NAME","OBSERVATION.COUNT",
           "LOCALITY.ID","LOCALITY.TYPE","REVIEWED","APPROVED","STATE","COUNTY",
           "LATITUDE","LONGITUDE","OBSERVATION.DATE","OBSERVER.ID",
           "DURATION.MINUTES","EFFORT.DISTANCE.KM",
           "NUMBER.OBSERVERS","GROUP.IDENTIFIER","SAMPLING.EVENT.IDENTIFIER")

nms = read.delim(rawpath, nrows = 1, sep = "\t", header = T, quote = "", stringsAsFactors = F, 
                 na.strings = c(""," ",NA))
nms = names(nms)
nms[!(nms %in% preimp)] = "NULL"
nms[nms %in% preimp] = NA

# read data from certain columns only

data = read.delim(rawpath, colClasses = nms, sep = "\t", header = T, quote = "", 
                  stringsAsFactors = F, na.strings = c(""," ",NA))

# read sensitive species data

nms = nms[-47]
sesp = read.csv(sensitivepath, colClasses = nms, stringsAsFactors = F)
stdformat = data.frame(date = as.character(sesp$OBSERVATION.DATE))
stdformat = stdformat %>%
  separate(date, c("month","day","year"), "/")
stdformat$year = as.numeric(stdformat$year)
sesp$OBSERVATION.DATE = paste(stdformat$year,"-",stdformat$month,"-",stdformat$day, sep = "")
sesp = sesp %>% mutate(GROUP.IDENTIFIER = ifelse(GROUP.IDENTIFIER == "", NA, GROUP.IDENTIFIER))

## merge both data frames

data = rbind(data,sesp)

## subset atlas lists

data = data %>% filter(SAMPLING.EVENT.IDENTIFIER %in% atlas, STATE == "Kerala")

## no of days in every month, and cumulative number

days = c(31,28,31,30,31,30,31,31,30,31,30,31)
cdays = c(0,31,59,90,120,151,181,212,243,273,304,334)

## create a column "group.id" which can help remove duplicate checklists

data = data %>%
  mutate(group.id = ifelse(is.na(GROUP.IDENTIFIER), SAMPLING.EVENT.IDENTIFIER, GROUP.IDENTIFIER))

## modify data

data = data %>%
  group_by(group.id,COMMON.NAME) %>% slice(1) %>% ungroup %>%
  mutate(OBSERVATION.DATE = as.Date(OBSERVATION.DATE), 
         month = month(OBSERVATION.DATE),
         day = day(OBSERVATION.DATE) + cdays[month], 
         #week = week(OBSERVATION.DATE),
         #fort = ceiling(day/14),
         year = year(OBSERVATION.DATE)) %>%
  dplyr::select(-c("OBSERVATION.DATE")) %>%
  ungroup
data  = data %>% filter(month %in% c(1:3,7:9)) %>% mutate(season = ifelse(month %in% 1:3, "D", "W"))

klbase = raster("KeralaDEMBounadary1.tif")

bb = bbox(klbase)
cs = c(6*1000/111111,6*1000/111111)  # cell size 9 km
cc = bb[, 1] + (cs/2)  # cell offset
cd = ceiling(diff(t(bb))/cs)  # number of cells per direction
grd = GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd) # create required grids
sp_grd = SpatialGridDataFrame(grd, data=data.frame(id=1:prod(cd))) # create spatial grid data frame
sp_grd_poly = as(sp_grd, "SpatialPolygonsDataFrame") # SGDF to SPDF

temp = data %>% group_by(group.id) %>% slice(1)

rownames(temp) = temp$group.id
coordinates(temp) = ~LONGITUDE + LATITUDE
temp = over(temp,sp_grd_poly)
temp = data.frame(temp)
temp$group.id = rownames(temp)
data = left_join(temp,data)
names(data)[1] = "gridkm"

gdat = data %>% distinct(COMMON.NAME,gridkm,season)
gdat = gdat %>% group_by(COMMON.NAME,gridkm) %>% summarize(nseasons = n())
data = left_join(data,gdat)
data  = data %>% mutate(season = ifelse(nseasons == 1, season, "R"))
data = data %>% filter(REVIEWED == 0 | APPROVED == 1) %>% 
  dplyr::select(-LOCALITY.ID,-LOCALITY.TYPE,-OBSERVER.ID,-SAMPLING.EVENT.IDENTIFIER,-DURATION.MINUTES,
         -EFFORT.DISTANCE.KM,-NUMBER.OBSERVERS,-GROUP.IDENTIFIER,-APPROVED,
         -REVIEWED,-nseasons)



gridmap = sp_grd_poly

#klsp = klbase > -Inf
#klsp = rasterToPolygons(klsp, dissolve=TRUE)

#bb = bbox(klsp)
#cs = c((25*1000)/111111,(25*1000)/111111)  # cell size grid*grid
#cc = bb[, 1] + (cs/2)  # cell offset
#cd = ceiling(diff(t(bb))/cs)  # number of cells per direction
#grd = GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
#sp_grd = SpatialGridDataFrame(grd, data=data.frame(id=1:prod(cd)))
#sp_grd_poly = as(sp_grd, "SpatialPolygonsDataFrame")
#mask = sp_grd_poly - klsp # SPDF that excludes KL map
#mask = crop(mask,extent(klbase))

save(data,gridmap,klbase,klsp,file = "datatoplot.RData")








################################### plotting

library(tidyverse)
library(ggfortify)
library(raster)
library(rgdal)
library(extrafont)
library(ggthemes)
library(rgeos)
library(scales)
library(grid)
library(gtable)
library(ggsci)
library(ggnewscale)
theme_set(theme_tufte())

load("datatoplot.RData")

species = "Green Warbler"

klbase_df = as.data.frame(klbase, xy = TRUE)
#klspdf = fortify(klsp)

plotbase = ggplot() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_raster(data = klbase_df , aes(x = x, y = y, fill = KeralaDEMBounadary1)) +
  scale_fill_material("brown", na.value = "white", reverse = F) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 12)) +
  theme(text=element_text(family="Gill Sans MT")) +
  theme(legend.text = element_text(size = 12)) +
  coord_quickmap()

data  = data %>%
  filter(!is.na(gridkm))

data$gridkm = as.factor(data$gridkm)

temp = data %>% 
  group_by(gridkm) %>%
  mutate(lists = n_distinct(group.id)) %>% ungroup() %>%
  filter(COMMON.NAME == species) %>%
  group_by(gridkm) %>%
  summarize(freq = n_distinct(group.id)/max(lists))

seasons = data %>% filter(COMMON.NAME == species) %>% distinct(gridkm,season)


fortified = fortify(gridmap, region = c("id"))
fortified$id = as.factor(fortified$id)
plotdf = na.omit(left_join(fortified,temp, by = c('id' = "gridkm"))) # SPDF to plot
plotdf = na.omit(left_join(plotdf,seasons, by = c('id' = "gridkm"))) # SPDF to plot

plotdfdry = plotdf %>% filter(season == "D")
plotdfwet = plotdf %>% filter(season == "W")
plotdfboth = plotdf %>% filter(season == "R")

## scale_fill_material colour schemes
# https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html

plotfin = plotbase +
  new_scale("fill") +
  theme(text=element_text(family="Gill Sans MT")) +
  geom_polygon(data = plotdfdry, aes(x = long, y = lat, group = group, fill = freq), alpha = 0.8) +
  scale_fill_material("indigo", na.value = "white") +
  new_scale("fill") +
  theme(text=element_text(family="Gill Sans MT")) +
  geom_polygon(data = plotdfwet, aes(x = long, y = lat, group = group, fill = freq), alpha = 0.8) +
  scale_fill_material("red", na.value = "white") +
  new_scale("fill") +
  theme(text=element_text(family="Gill Sans MT")) +
  geom_polygon(data = plotdfboth, aes(x = long, y = lat, group = group, fill = freq), alpha = 0.8) +
  scale_fill_material("green", na.value = "white") +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.position = "none")
  

# plot
n2 = paste(species,".png",sep="")

png(n2, units="in", width=10, height=7, res=1000)
grid::grid.draw(plotfin)
dev.off()





########################### high resolution occupancy based plots

library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(sp)
library(sf)
library(rgeos)
require(reshape2)
require(data.table)
require(unmarked)

rawpath = "ebd_IN_relSep-2019.txt"
sensitivepath = "Sensitive_India_may 2019.csv"

preimp = c("CATEGORY","COMMON.NAME","SCIENTIFIC.NAME","OBSERVATION.COUNT",
           "LOCALITY.ID","LOCALITY.TYPE","REVIEWED","APPROVED","STATE","COUNTY",
           "LATITUDE","LONGITUDE","OBSERVATION.DATE","OBSERVER.ID",
           "DURATION.MINUTES","EFFORT.DISTANCE.KM","PROTOCOL.TYPE",
           "NUMBER.OBSERVERS","GROUP.IDENTIFIER","SAMPLING.EVENT.IDENTIFIER")

nms = read.delim(rawpath, nrows = 1, sep = "\t", header = T, quote = "", stringsAsFactors = F, 
                 na.strings = c(""," ",NA))
nms = names(nms)
nms[!(nms %in% preimp)] = "NULL"
nms[nms %in% preimp] = NA

# read data from certain columns only

data = read.delim(rawpath, colClasses = nms, sep = "\t", header = T, quote = "", 
                  stringsAsFactors = F, na.strings = c(""," ",NA))

# read sensitive species data

nms = nms[-47]
sesp = read.csv(sensitivepath, colClasses = nms, stringsAsFactors = F)
stdformat = data.frame(date = as.character(sesp$OBSERVATION.DATE))
stdformat = stdformat %>%
  separate(date, c("month","day","year"), "/")
stdformat$year = as.numeric(stdformat$year)
sesp$OBSERVATION.DATE = paste(stdformat$year,"-",stdformat$month,"-",stdformat$day, sep = "")
sesp = sesp %>% mutate(GROUP.IDENTIFIER = ifelse(GROUP.IDENTIFIER == "", NA, GROUP.IDENTIFIER))

## merge both data frames

data = rbind(data,sesp)

## subset KL lists

data = data %>% filter(STATE == "Kerala", PROTOCOL.TYPE == "Stationary" | 
                         (PROTOCOL.TYPE == "Traveling" & EFFORT.DISTANCE.KM < 0.3))


# no of days in every month, and cumulative number
days = c(31,28,31,30,31,30,31,31,30,31,30,31)
cdays = c(0,31,59,90,120,151,181,212,243,273,304,334)

# create a column "group.id" which can help remove duplicate checklists
data = data %>%
  mutate(group.id = ifelse(is.na(GROUP.IDENTIFIER), SAMPLING.EVENT.IDENTIFIER, GROUP.IDENTIFIER))

imp = c("COMMON.NAME","OBSERVATION.COUNT",
        "COUNTY","LOCALITY.ID",
        "LATITUDE","LONGITUDE","OBSERVATION.DATE",
        "group.id")

data = data %>%
  group_by(group.id,COMMON.NAME) %>% slice(1) %>% ungroup %>%
  dplyr::select(imp) %>%
  mutate(OBSERVATION.DATE = as.Date(OBSERVATION.DATE), 
         month = month(OBSERVATION.DATE),
         day = day(OBSERVATION.DATE) + cdays[month], 
         #week = week(OBSERVATION.DATE),
         #fort = ceiling(day/14),
         cyear = year(OBSERVATION.DATE)) %>%
  dplyr::select(-c("OBSERVATION.DATE")) %>%
  mutate(year = ifelse(day <= 151, cyear-1, cyear)) %>%
  group_by(group.id) %>% mutate(no.sp = n_distinct(COMMON.NAME)) %>%
  ungroup

data = data %>%
  filter(year > 2013)

data = data %>%
  mutate(OBSERVATION.COUNT = replace(OBSERVATION.COUNT, !is.na(OBSERVATION.COUNT), "1"))

data$OBSERVATION.COUNT = as.numeric(data$OBSERVATION.COUNT)

klbase = raster("KeralaDEMBounadary1.tif")
klbase_df = as.data.frame(klbase, xy = TRUE)


klextract = extract(klbase, SpatialPoints(klbase), sp = T)

#klsp = klbase > -Inf
#klsp = rasterToPolygons(klsp, dissolve=FALSE)


temp = data %>% group_by(LOCALITY.ID) %>% slice(1)

rownames(temp) = temp$LOCALITY.ID
coordinates(temp) = ~LONGITUDE + LATITUDE
temp = over(temp,klextract)
temp = data.frame(temp)
temp$LOCALITY.ID = rownames(temp)
data1 = left_join(temp,data)
names(data1)[1] = "elevation"







###################### for Praveen

plotkl = function(path = "klmap.csv", a = 6, b = 5, c = -1, d = 2.1, exp1 = -1.5)
{
  require(tidyverse)
  require(ggfortify)
  require(raster)
  require(rgdal)
  require(viridis)
  require(extrafont)
  require(ggthemes)
  require(rgeos)
  require(scales)
  require(grid)
  require(gtable)
  require(ggsci)
  require(ggnewscale)
  theme_set(theme_tufte())
  
  klbase_df = read.csv(path)
  
  klbase_df$long = scale(klbase_df$y, center = F)
  klbase_df$alt = scale(log(klbase_df$elevation+80), center = F)
  
  ## here I have created a frequency column using a simple function including altitude and longitude
  ## simulate any other data as you see fit
  
  klbase_df = klbase_df %>%
    mutate(freq = (1/(1+exp(c + a*alt^exp1 + b*(long^2-d*long)))))
  klbase_df$freq = round(klbase_df$freq,1)
  
  klbase = rasterFromXYZ(klbase_df[,c(1,2,5)])
  klbase = focal(klbase, w=matrix(1,nrow=3,ncol=3), fun=mean)
  klbase_df = as.data.frame(klbase, xy = TRUE)
  names(klbase_df)[3] = "freq"
  
  # I have used this viridis scale but you can use scale_fill_material instead or any others
  
  ## scale_fill_material colour schemes
  # https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
  
  klbase_df$freq = klbase_df$freq - 0.3
  klbase_df$freq[klbase_df$freq < -0.1] = NA
  klbase_df$freq[klbase_df$freq < 0] = 0

  plotbase = ggplot() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    geom_raster(data = klbase_df , aes(x = x, y = y, fill = freq)) +
    #scale_fill_material("brown", na.value = "white", reverse = F) +
    scale_fill_viridis(na.value = "grey") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 12)) +
    theme(text=element_text(family="Gill Sans MT")) +
    theme(legend.text = element_text(size = 12)) +
    #theme(legend.position = "none") +
    coord_quickmap()
  
  ########## to plot a new scale, create another column freq2 and you can overlay the colour scale
  
  #plotbase = plotbase +
    #new_scale("fill") +
    #theme(text=element_text(family="Gill Sans MT")) +
    #geom_raster(data = klbase_df , aes(x = x, y = y, fill = freq2)) +
    #scale_fill_material("indigo", na.value = "white")
  
  return(plotbase)
}

plotkl(a = 8, b = 8, c = -1, d = 1.9, exp1 = -2)







