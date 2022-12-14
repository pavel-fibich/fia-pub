rm(list=ls())

# read FIA and envi data, format the data

library(rFIA)
getFIA("FL", dir = "./states", common = F, tables = c("COND") )
for (op in foplc[1,2]) getFIA(op, dir = "./states", common = F, tables = c("PLOT","TREE") )
for (op in foplc[,2]) getFIA(op, dir = "./states", common = F, tables = c("COND") )
db <- readFIA("./states", tables = c("PLOT","TREE") )


foplots=c("Florida","Texas","Louisiana","Arkansas","Mississippi","Alabama","Georgia","South Carolina",
          "North Carolina","Tennessee","Kentucky","Virginia","West Virginia","Delaware","Maryland","New Jersey",
          "Pennsylvania","Connecticut","New York","Massachusetts","New Hampshire","Rhode Island",
          "Vermont","Maine")
survs<-read.csv("states_all/SURVEY.csv");names(survs)                                          
foplc<-merge(data.frame(STATENM=foplots),survs[,c(5:6,4)],all.x=T); foplc<-foplc[!duplicated(foplc),]; foplc
rm(survs,foplots)

allp <- NULL
dist52<-FALSE

for (pl in foplc[,2]) {
  dbpl<-read.csv(paste("states/",pl,"_PLOT.csv",sep=""));dim(dbpl)
  dbco<-read.csv(paste("states/",pl,"_COND.csv",sep=""));dim(dbco)
  names(dbco)
  head(names(dbpl));head(names(dbco))
  names(dbco)[1:2]=c("CNx","CN")
  dbpl<-merge(dbpl,dbco[,c(2,grep("DSTRBCD1$",names(dbco)),grep("TRTCD1$",names(dbco)),
                           grep("FOREST_COMMUNITY_PNWRS",names(dbco)),grep("BALIVE",names(dbco)),
                           grep("DOMINANT_SPECIES1_PNWRS",names(dbco)),
                           grep("SLOPE",names(dbco)), grep("ASPECT",names(dbco)), grep("PHYSCLCD",names(dbco)),
                           grep("ELEV$",names(dbco)), grep("LIVE_CANOPY_CVR_PCT",names(dbco)) ) ], all.x=T)
  rm(dbco)
  #table(dbpl$DSTRBCD1);table(dbpl$TRTCD1);table(dbpl$FOREST_COMMUNITY_PNWRS)
  if (dist52) {
    dbpl <- dbpl[ ( dbpl$DSTRBCD1 %in% c(52) )  # disturbances
                  & (is.na(dbpl$TRTCD1) | (dbpl$TRTCD1 == 0)) # no forestry treatment
                  & (is.na(dbpl$FOREST_COMMUNITY_PNWRS) | (!dbpl$FOREST_COMMUNITY_PNWRS %in% c(8,9))),] # not agri and plantations  
  } else {
    dbpl <- dbpl[ (is.na(dbpl$DSTRBCD1) | (!dbpl$DSTRBCD1 %in% c(80, 30, 31, 32, 46, 52)))  # disturbances
                  & (is.na(dbpl$TRTCD1) | (dbpl$TRTCD1 == 0)) # no forestry treatment
                  & (is.na(dbpl$FOREST_COMMUNITY_PNWRS) | (!dbpl$FOREST_COMMUNITY_PNWRS %in% c(8,9))),] # not agri and plantations
  }
  print(pl)
  if ( nrow(dbpl)< 5) next;
  print(table(dbpl$DSTRBCD1))
  #}
  dbpl$latlon<-paste0(dbpl$LAT,"x",dbpl$LON)
  reme<-table(dbpl$latlon); (remet<-table(reme)); sum(remet)
  table(dbpl$MEASYEAR)
  
  # trees
  dbtr<-read.csv(paste("states/",pl,"_TREE.csv",sep=""));dim(dbtr)
  dbtr<-dbtr[ dbtr$PLT_CN %in% dbpl$CN,]
  dbtr$DIAHT<-dbtr$DIA/dbtr$HT
  # add weighted mean by BA
  dbtr$BA <- pi*(dbtr$DIA/2)^2
  
  # aggregate trees chars
  
  #x<-dbtr[dbtr$PLT_CN == dbtr[500,]$PLT_CN,]
  
  maxaht<-aggregate(ACTUALHT ~ PLT_CN,data=dbtr,FUN=max,na.rm=T)
  names(maxaht) <- c("CN","maxAHT")
  meaaht<-aggregate(dbtr$ACTUALHT ~ dbtr$PLT_CN,FUN=mean,na.rm=T)
  names(meaaht) <- c("CN","meanAHT")
  maxaht9<-aggregate(ACTUALHT ~ PLT_CN,data=dbtr,FUN=quantile,probs=c(0.99),na.rm=T)
  names(maxaht9) <- c("CN","maxAHT9")
  maxaht5<-aggregate(ACTUALHT ~ PLT_CN,data=dbtr,FUN=quantile,probs=c(0.95),na.rm=T)
  names(maxaht5) <- c("CN","maxAHT5")
  maxahtsk<-aggregate(ACTUALHT ~ PLT_CN,data=dbtr,FUN=skewness,na.rm=T)
  names(maxahtsk) <- c("CN","maxAHTsk")
  maxahtskk<-aggregate(ACTUALHT ~ PLT_CN,data=dbtr,FUN=kurtosis,na.rm=T)
  names(maxahtskk) <- c("CN","maxAHTskk")
  
  maxht<-aggregate(HT ~ PLT_CN,data=dbtr,FUN=max,na.rm=T)
  names(maxht) <- c("CN","maxHT")
  meaht<-aggregate(dbtr$HT ~ dbtr$PLT_CN,FUN=mean,na.rm=T)
  names(meaht) <- c("CN","meanHT")
  maxht9<-aggregate(HT ~ PLT_CN,data=dbtr,FUN=quantile,probs=c(0.99),na.rm=T)
  names(maxht9) <- c("CN","maxHT9")
  maxht5<-aggregate(HT ~ PLT_CN,data=dbtr,FUN=quantile,probs=c(0.95),na.rm=T)
  names(maxht5) <- c("CN","maxHT5")
  maxhtsk<-aggregate(HT ~ PLT_CN,data=dbtr,FUN=skewness,na.rm=T)
  names(maxhtsk) <- c("CN","maxHTsk")
  maxhtskk<-aggregate(HT ~ PLT_CN,data=dbtr,FUN=kurtosis,na.rm=T)
  names(maxhtskk) <- c("CN","maxHTskk")
  
  maxdia<-aggregate(dbtr$DIA ~ dbtr$PLT_CN,FUN=max,na.rm=T)
  names(maxdia) <- c("CN","maxDIA")
  meadia<-aggregate(dbtr$DIA ~ dbtr$PLT_CN,FUN=mean,na.rm=T)
  names(meadia) <- c("CN","meanDIA")
  maxdia9<-aggregate(DIA ~ PLT_CN,data=dbtr,FUN=quantile,probs=c(0.99),na.rm=T)
  names(maxdia9) <- c("CN","maxDIA9")
  maxdia5<-aggregate(DIA ~ PLT_CN,data=dbtr,FUN=quantile,probs=c(0.95),na.rm=T)
  names(maxdia5) <- c("CN","maxDIA5")
  maxdiask<-aggregate(DIA ~ PLT_CN,data=dbtr,FUN=skewness,na.rm=T)
  names(maxdiask) <- c("CN","maxDIAsk")
  maxdiaskk<-aggregate(DIA ~ PLT_CN,data=dbtr,FUN=kurtosis,na.rm=T)
  names(maxdiaskk) <- c("CN","maxDIAskk")
  
  maxdiaht<-aggregate(dbtr$DIAHT ~ dbtr$PLT_CN,FUN=max,na.rm=T)
  names(maxdiaht) <- c("CN","maxDIAHT")
  meadiaht<-aggregate(dbtr$DIAHT ~ dbtr$PLT_CN,FUN=mean,na.rm=T)
  names(meadiaht) <- c("CN","meanDIAHT")
  
  sumba<-aggregate(dbtr$BA ~ dbtr$PLT_CN,FUN=sum,na.rm=T)
  names(sumba) <- c("CN","sumBA")
  
  sumno<-aggregate(dbtr$CN ~ dbtr$PLT_CN,FUN=function(x) {length(unique(x,na.rm=T))})
  names(sumno) <- c("CN","sumNO")
  
  nospe<-aggregate(dbtr$SPCD ~ dbtr$PLT_CN,FUN=function(x){ return(length(unique(x))) })
  names(nospe) <- c("CN","noSPE")
  nospes<-aggregate(dbtr$SPCD ~ dbtr$PLT_CN,FUN=function(x){ return( diversity(table(x)) )  })
  names(nospes) <- c("CN","noSPEsha")
  nospee<-aggregate(dbtr$SPCD ~ dbtr$PLT_CN,FUN=function(x){ return( diversity(table(x))/ log(length(unique(x))) )  })
  names(nospee) <- c("CN","noSPEeve")
  
  dbpl <- merge(dbpl,maxht, all.x=T)
  dbpl <- merge(dbpl,maxht9, all.x=T)
  dbpl <- merge(dbpl,maxht5, all.x=T)
  dbpl <- merge(dbpl,maxhtsk, all.x=T)
  dbpl <- merge(dbpl,maxhtskk, all.x=T)
  dbpl <- merge(dbpl,meaht, all.x=T)
  dbpl <- merge(dbpl,maxaht, all.x=T)
  dbpl <- merge(dbpl,maxaht9, all.x=T)
  dbpl <- merge(dbpl,maxaht5, all.x=T)
  dbpl <- merge(dbpl,maxahtsk, all.x=T)
  dbpl <- merge(dbpl,maxahtskk, all.x=T)
  dbpl <- merge(dbpl,meaaht, all.x=T)
  dbpl <- merge(dbpl,maxdia, all.x=T)
  dbpl <- merge(dbpl,maxdia9, all.x=T)
  dbpl <- merge(dbpl,maxdia5, all.x=T)
  dbpl <- merge(dbpl,maxdiask, all.x=T)
  dbpl <- merge(dbpl,maxdiaskk, all.x=T)
  dbpl <- merge(dbpl,meadia, all.x=T)
  dbpl <- merge(dbpl,maxdiaht, all.x=T)
  dbpl <- merge(dbpl,meadiaht, all.x=T)
  dbpl <- merge(dbpl,sumba, all.x=T)
  dbpl <- merge(dbpl,sumno, all.x=T)
  dbpl <- merge(dbpl,nospe, all.x=T)
  dbpl <- merge(dbpl,nospes, all.x=T)
  dbpl <- merge(dbpl,nospee, all.x=T)
  
  # seedlings
  dbse<-read.csv(paste("states/",pl,"_SEEDLING.csv",sep=""));dim(dbse)
  dbse<-dbse[ dbse$PLT_CN %in% dbpl$CN,]
  
  sumse<-aggregate(dbse$TREECOUNT ~ dbse$PLT_CN,FUN=sum,na.rm=T)
  names(sumse) <- c("CN","sumSEEtc")
  sumse2<-aggregate(dbse$CN ~ dbse$PLT_CN,FUN=function(x) { return( length(unique(x)))})
  names(sumse2) <- c("CN","sumSEEun")
  sumse3<-aggregate(dbse$SPCD ~ dbse$PLT_CN,FUN=function(x) { return( length(unique(x)))})
  names(sumse3) <- c("CN","sumSEEspe")
  sumse4<-aggregate(dbse$SPCD ~ dbse$PLT_CN,FUN=function(x) { return( diversity(table(x)) )})
  names(sumse4) <- c("CN","sumSEEsha")
  sumse5<-aggregate(dbse$SPCD ~ dbse$PLT_CN,FUN=function(x) { return( diversity(table(x))/ log(length(unique(x))) )})
  names(sumse5) <- c("CN","sumSEEeve")
  
  dbpl <- merge(dbpl,sumse, all.x=T)
  dbpl <- merge(dbpl,sumse2, all.x=T)
  dbpl <- merge(dbpl,sumse3, all.x=T)
  dbpl <- merge(dbpl,sumse4, all.x=T)
  dbpl <- merge(dbpl,sumse5, all.x=T)
  if (dist52) {  
    write.csv(dbpl,paste("states_sel/",pl,"_PLOTselDIS52.csv",sep=""),row.names=F)
  } else{
    write.csv(dbpl,paste("states_sel/",pl,"_PLOTselno52.csv",sep=""),row.names=F)
  }
  rm(dbtr,maxdia,meadia,maxht,meaht,maxdiaht,meadiaht, sumba,sumno,nospe, nospes,nospee, dbse,sumse,sumse2,sumse3,sumse4,sumse5)
  if ( is.null(allp) ){
    allp <-dbpl
    #allt<-dbtr
  } else{
    allp <- rbind(allp,dbpl)
    #allt <- rbind(allt,dbtr)
  }
}


# worldclim
library(sp)
library(raster)
wc<-getData("worldclim",var="bio",res=10) # 10 minutes resolution
#human impact
hfp <- raster("HFP2009.tif")
#tracks
#https://www.ncdc.noaa.gov/ibtracs/index.php?name=ib-v4-access
library(geosphere)
library(ncdf4)
ct<-nc_open("cyclones/IBTrACS.NA.v04r00.nc")
ct_lon <-ncvar_get(ct,"lon") 
ct_lat <-ncvar_get(ct,"lat")
ct_wind <-ncvar_get(ct,"usa_wind")
ct_sshs <-ncvar_get(ct,"usa_sshs")
ct_year <-ncvar_get(ct,"season")
last40<-(ct_year <=2019) & (ct_year>=1980)
windmax<-apply(ct_sshs,2,function(x) { any(x>=1,na.rm=T) }); table(windmax)
table(last40&windmax)


xy <- SpatialPoints(dbpl[ cooplot,c("LON","LAT")]) 
crs(xy) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")   #priradi geograficky system (WGS84)
xy; #plot(xy)


wc_2 <- raster::extract(wc,y=xy,method="bilinear")
wc_h <- raster::extract(hfp,y=xy,method="bilinear")
dbpl[cooplot, grep("^bio",names(dbpl))] = round(wc_2,2)
dbpl[cooplot,]$HFP = round(wc_h,2)

dbplx<-dbpl[(cooplot) & (!duplicated(dbpl$latlon)),]

# add ibtracs
onc<-unique(oncm$id)
i=onc[155]
for( i in onc[1:length(onc)] ){
  print(paste0(i," ",round( which(i==onc) / length(onc) * 100), '% completed'))
  oncmi <- oncm[ oncm$id == i,]
  oncmi <- oncmi[(!is.na(oncmi$x)) & (!is.na(oncmi$y)), ]
  if (nrow(oncmi) == 1)
    oncmi <- rbind(oncmi,oncmi)
  fceck <- distm ( y=dbplx[,c("LON","LAT")], x=oncmi[,1:2], fun = distHaversine)
  if ( (any(fceck<(updi2+updi1))) ){
    cloplo<-which(colSums(fceck<(updi2+updi1),na.rm=T )>0)
    ddic1 <-dist2Line(p = dbplx[cloplo,c("LON","LAT")],line=oncmi[,c("x","y")], distfun=distHaversine )
    if ( any(ddic1[,1] < updi2) ) { # up to 200km
      ddic2 <- which( ddic1[,1] < updi2)#858000)
      if (length(ddic2)>=1) {
        # increment freq - by default is 0
        dbplx[cloplo[ddic2],]$ITf200 <- dbplx[ cloplo[ddic2],]$ITf200 + 1
        cloplo2<-colSums(fceck<(updi2),na.rm=T )
        dbplx$ITfp200 <- dbplx$ITfp200 + cloplo2
        # distance to closest point on cyclone
        clo2 <- distm ( y=ddic1[ddic2,2:3], x=oncmi[,1:2], fun = distHaversine ) 
        clo2min <- apply(clo2,2,which.min)
        clo2wind <- oncmi[clo2min,]$z
        # add max wind speed - default is NA
        dbplx[cloplo[ddic2],]$ITwm200 <- pmax(dbplx[cloplo[ddic2],]$ITwm200,clo2wind,na.rm=T)
        # add max wind speed - default is NA
        dbplx[cloplo[ddic2],]$ITws200 <- rowSums(cbind(dbplx[cloplo[ddic2],]$ITws200,clo2wind),na.rm=T)
        dbplx[cloplo[ddic2],]$ITfs200 <- dbplx[cloplo[ddic2],]$ITfs200 + (!is.na(clo2wind))
        
        if ( any(ddic1[,1] < updi1) ) { # up to 100km
          ddic2 <- which( ddic1[,1] < updi1)#878000)
          if (length(ddic2)>=1) {
            # increment freq - by default is 0
            dbplx[cloplo[ddic2],]$ITf100 <- dbplx[cloplo[ddic2],]$ITf100 + 1
            cloplo1<-colSums(fceck<(updi1),na.rm=T )
            dbplx$ITfp100 <- dbplx$ITfp100 + cloplo1
            # distance to closest point on cyclone
            clo2 <- distm ( y=ddic1[ddic2,2:3], x=oncmi[,1:2], fun = distGeo)
            clo2min <- apply(clo2,2,which.min)
            clo2wind <- oncmi[clo2min,]$z
            # add max wind speed - default is NA
            dbplx[cloplo[ddic2],]$ITwm100 <- pmax(dbplx[cloplo[ddic2],]$ITwm100,clo2wind,na.rm=T)
            # add mean wind speed - default is NA
            dbplx[cloplo[ddic2],]$ITws100 <- rowSums(cbind(dbplx[cloplo[ddic2],]$ITws100,clo2wind),na.rm=T)
            dbplx[cloplo[ddic2],]$ITfs100 <- dbplx[cloplo[ddic2],]$ITfs100 + (!is.na(clo2wind))
          }
        }
      }
    }
  }
}


## distance to coastline
library(rnaturalearth)
library(sf)
library(raster)
library(rgdal)
library(sp)
library(tidyverse)
library(geosphere)
coastline <- ne_coastline(scale = 110,returnclass="sf")
coastline <- st_transform(coastline, 4326)
coastline <- st_cast(coastline, "MULTILINESTRING")
# crop just interesting part
box = c(xmin = -110, ymin = 20, xmax = -60, ymax = 50)
coastline1<-st_crop(coastline, box)
coastline1 <- st_transform(coastline1, 4326)
coastline1 <- st_cast(coastline1, "MULTILINESTRING")
d1_sf <- dbplx %>% st_as_sf(coords = c("LON","LAT")) %>% st_set_crs(crs(coastline))

# get distance 
#ddic<-geosphere::dist2Line(p = st_coordinates(d1_sf),line=as(coastline,'Spatial') )
ddic1<-geosphere::dist2Line(p = st_coordinates(d1_sf),line=as(coastline1,'Spatial') )
dbplx$CoastDist <- ddic1[,1]

table(dbplx$ITf200>=10)
dbplx<-dbplx[ dbplx$ITf200>=10, ]



# aggregation

round2<-function(x,decs=0){
  return ( round(x/0.5,decs)*0.5 )
}


dbplxT<-dbplx

sround<-0
dbplx$LAT1 <- round(dbplx$LAT,sround)
dbplx$LON1 <- round(dbplx$LON,sround)
# expla
(onag<-paste0("cbind(",paste0(oncy,collapse = ","),") ~ LON1+LAT1"))
dbpla<-aggregate( as.formula(onag),data=dbplx,FUN=mean,na.rm=T,na.action = na.omit)
names(dbpla)
#x<-dbplx[,oncy]
#colSums(is.na(x))
# max
#onagm<-"cbind(maxHT,maxAHT,maxDIA,maxHT9,maxAHT9,maxDIA9,maxHT5,maxAHT5,maxDIA5,maxHTsk,maxAHTsk,maxDIAsk,maxDIAHT,one) ~ LON1+LAT1"
fpa<-grep("noSPEeve",yrep)
(onagm <- paste0("cbind(",paste0(yrep[1:fpa],collapse=","),",one) ~ LON1+LAT1"))
dbpla1<-aggregate( as.formula(onagm), data=dbplx, FUN=mean,na.rm=T)
names(dbpla1)
#x<-dbplx[,yrep]
#colSums(!is.na(x))
#mean
#onagme<-"cbind(SLOPE,ASPECT,ELEV, BALIVE,meanHT,meanAHT,meanDIA,meanDIAHT,sumBA,sumNO,noSPE,noSPEsha,noSPEeve,sumSEEun,sumSEEspe,sumSEEsha,sumSEEeve,two) ~ LON1+LAT1"
(onagme<-paste0("cbind(",paste0(yrep[(fpa+1):length(yrep)],collapse=","),",SLOPE,ASPECT,ELEV,CoastDist,two) ~ LON1+LAT1"))
dbpla2<-aggregate( as.formula(onagme), data=dbplx, FUN=mean,na.rm=T)
names(dbpla2)
dbpla = merge(dbpla,dbpla1,all.x=T)
dbpla = merge(dbpla,dbpla2,all.x=T)

sround<-0
dbplx$LAT1 <- round2(dbplx$LAT,sround)
dbplx$LON1 <- round2(dbplx$LON,sround)
# expla
# dbplaa<-aggregate( as.formula("cbind(ITf200,ITf100)~LON1+LAT1"), data=dbplx,FUN=mean,na.rm=T)
# write.csv(dbplaa,"rawITdata05aggALL.csv",row.names=F)
# dbplaa<-aggregate( as.formula(onag), data=dbplx,FUN=sum,na.rm=T)
dbplaa<-aggregate( as.formula(onag), data=dbplx,FUN=mean,na.rm=T)
names(dbplaa)
# max
dbplaa1<-aggregate( as.formula(onagm), data=dbplx, FUN=mean,na.rm=T)
names(dbplaa1)
#mean
dbplaa2<-aggregate( as.formula(onagme),data=dbplx, FUN=mean,na.rm=T, drop=F)
names(dbplaa2)
dbplaa = merge(dbplaa,dbplaa1,all.x=T)
dbplaa = merge(dbplaa,dbplaa2,all.x=T)

if (onpl != "ALL") {
  if (is.null(alldbpla) ){
    dbpla$state<-onpl
    alldbpla<-dbpla
    dbplaa$state<-onpl
    alldbpla2<-dbplaa
  } else {
    dbpla$state<-onpl
    alldbpla<-rbind(alldbpla,dbpla)
    dbplaa$state<-onpl
    alldbpla2<-rbind(alldbpla2,dbplaa)
  }
}

save(dbpl,file=paste0("aggres",fsu,"/dbplPc",onpl,".RData"))
save(dbpla,file=paste0("aggres",fsu,"/dbplPc",onpl,"_agg.RData"))
save(dbplaa,file=paste0("aggres",fsu,"/dbplPc",onpl,"_agg2.RData"))

#type 1
#type 1
#type 1
ontype<-"1"; dbplx<-dbplxT[(!duplicated(dbplxT$latlon)) & (dbplxT$SiteType ==ontype),]

sround<-0
dbplx$LAT1 <- round(dbplx$LAT,sround)
dbplx$LON1 <- round(dbplx$LON,sround)
# expla
#onag<-paste0("cbind(",paste0(oncy,collapse = ","),") ~ LON1+LAT1")

dbpla<-aggregate( as.formula(onag),data=dbplx[rowSums(!is.na(dbplx[,oncy]))>2,],FUN=mean,na.rm=T,na.action = na.omit) 
names(dbpla)
x<-dbplx[rowSums(!is.na(dbplx[,oncy]))>2,]
x<-dbplx[,oncy]
colSums(!is.na(x))

# max
#onagm<-"cbind(maxHT,maxAHT,maxDIA,maxDIAHT) ~ LON1+LAT1"
dbpla1<-aggregate( as.formula(onagm), data=dbplx, FUN=mean,na.rm=T)
names(dbpla1)
#mean
#onagme<-"cbind(SLOPE,ASPECT,ELEV, BALIVE,LiveCover,meanHT,meanAHT,meanDIA,meanDIAHT,sumBA,sumNO,noSPE,noSPEsha,noSPEeve,sumSEEtc,sumSEEun,sumSEEspe,sumSEEsha,sumSEEeve) ~ LON1+LAT1"
dbpla2<-aggregate( as.formula(onagme), data=dbplx, FUN=mean,na.rm=T, drop=F)
names(dbpla2)
dbpla = merge(dbpla,dbpla1,all.x=T)
dbpla = merge(dbpla,dbpla2,all.x=T)

sround<-0
dbplx$LAT1 <- round2(dbplx$LAT,sround)
dbplx$LON1 <- round2(dbplx$LON,sround)
# expla
dbplaa<-aggregate( as.formula(onag), data=dbplx,FUN=mean,na.rm=T)
names(dbplaa)
# max
dbplaa1<-aggregate( as.formula(onagm), data=dbplx, FUN=mean,na.rm=T)
names(dbplaa1)
#mean
dbplaa2<-aggregate( as.formula(onagme),data=dbplx, FUN=mean,na.rm=T, drop=F)
names(dbplaa2)
dbplaa = merge(dbplaa,dbplaa1,all.x=T)
dbplaa = merge(dbplaa,dbplaa2,all.x=T)

if (onpl != "ALL") {
  if (is.null(alldbpla) ){
    dbpla$state<-onpl
    alldbpla_1<-dbpla
    dbplaa$state<-onpl
    alldbpla2_1<-dbplaa
  } else {
    dbpla$state<-onpl
    alldbpla_1<-rbind(alldbpla_1,dbpla)
    dbplaa$state<-onpl
    alldbpla2_1<-rbind(alldbpla2_1,dbplaa)
  }
}
save(dbpla,file=paste0("aggres",fsu,"/dbplPc",onpl,"_agg_",ontype,".RData"))
save(dbplaa,file=paste0("aggres",fsu,"/dbplPc",onpl,"_agg2_",ontype,".RData"))


#type 2
#type 2
#type 2
ontype<-"2"; dbplx<-dbplxT[(!duplicated(dbplxT$latlon)) & (dbplxT$SiteType ==ontype),]

sround<-0
dbplx$LAT1 <- round(dbplx$LAT,sround)
dbplx$LON1 <- round(dbplx$LON,sround)
# expla
#onag<-paste0("cbind(",paste0(oncy,collapse = ","),") ~ LON1+LAT1")
dbpla<-aggregate( as.formula(onag),data=dbplx,FUN=mean,na.rm=T,na.action = na.omit)
names(dbpla)
# max
#onagm<-"cbind(maxHT,maxAHT,maxDIA,maxDIAHT) ~ LON1+LAT1"
dbpla1<-aggregate( as.formula(onagm), data=dbplx, FUN=mean,na.rm=T)
names(dbpla1)
#mean
#onagme<-"cbind(SLOPE,ASPECT,ELEV, BALIVE,LiveCover,meanHT,meanAHT,meanDIA,meanDIAHT,sumBA,sumNO,noSPE,noSPEsha,noSPEeve,sumSEEtc,sumSEEun,sumSEEspe,sumSEEsha,sumSEEeve) ~ LON1+LAT1"
dbpla2<-aggregate( as.formula(onagme), data=dbplx, FUN=mean,na.rm=T, drop=F)
names(dbpla2)
dbpla = merge(dbpla,dbpla1,all.x=T)
dbpla = merge(dbpla,dbpla2,all.x=T)

sround<-0
dbplx$LAT1 <- round2(dbplx$LAT,sround)
dbplx$LON1 <- round2(dbplx$LON,sround)
# expla
dbplaa<-aggregate( as.formula(onag), data=dbplx,FUN=mean,na.rm=T)
names(dbplaa)
# max
dbplaa1<-aggregate( as.formula(onagm), data=dbplx, FUN=mean,na.rm=T)
names(dbplaa1)
#mean
dbplaa2<-aggregate( as.formula(onagme),data=dbplx, FUN=mean,na.rm=T, drop=F)
names(dbplaa2)
dbplaa = merge(dbplaa,dbplaa1,all.x=T)
dbplaa = merge(dbplaa,dbplaa2,all.x=T)

if (onpl != "ALL") {
  if (is.null(alldbpla) ){
    dbpla$state<-onpl
    alldbpla_2<-dbpla
    dbplaa$state<-onpl
    alldbpla2_2<-dbplaa
  } else {
    dbpla$state<-onpl
    alldbpla_2<-rbind(alldbpla_1,dbpla)
    dbplaa$state<-onpl
    alldbpla2_2<-rbind(alldbpla2_1,dbplaa)
  }
}
save(dbpla,file=paste0("aggres",fsu,"/dbplPc",onpl,"_agg_",ontype,".RData"))
save(dbplaa,file=paste0("aggres",fsu,"/dbplPc",onpl,"_agg2_",ontype,".RData"))

#type 3
#type 3
#type 3
ontype<-"3"; dbplx<-dbplxT[(!duplicated(dbplxT$latlon)) & (dbplxT$SiteType ==ontype),]

sround<-0
dbplx$LAT1 <- round(dbplx$LAT,sround)
dbplx$LON1 <- round(dbplx$LON,sround)
# expla
#onag<-paste0("cbind(",paste0(oncy,collapse = ","),") ~ LON1+LAT1")
dbpla<-aggregate( as.formula(onag),data=dbplx,FUN=mean,na.rm=T,na.action = na.omit)
names(dbpla)
# max
#onagm<-"cbind(maxHT,maxAHT,maxDIA,maxDIAHT) ~ LON1+LAT1"
dbpla1<-aggregate( as.formula(onagm), data=dbplx, FUN=mean,na.rm=T)
names(dbpla1)
#mean
#onagme<-"cbind(SLOPE,ASPECT,ELEV, BALIVE,LiveCover,meanHT,meanAHT,meanDIA,meanDIAHT,sumBA,sumNO,noSPE,noSPEsha,noSPEeve,sumSEEtc,sumSEEun,sumSEEspe,sumSEEsha,sumSEEeve) ~ LON1+LAT1"
dbpla2<-aggregate( as.formula(onagme), data=dbplx, FUN=mean,na.rm=T, drop=F)
names(dbpla2)
dbpla = merge(dbpla,dbpla1,all.x=T)
dbpla = merge(dbpla,dbpla2,all.x=T)

sround<-0
dbplx$LAT1 <- round2(dbplx$LAT,sround)
dbplx$LON1 <- round2(dbplx$LON,sround)
# expla
dbplaa<-aggregate( as.formula(onag), data=dbplx,FUN=mean,na.rm=T)
names(dbplaa)
# max
dbplaa1<-aggregate( as.formula(onagm), data=dbplx, FUN=mean,na.rm=T)
names(dbplaa1)
#mean
dbplaa2<-aggregate( as.formula(onagme),data=dbplx, FUN=mean,na.rm=T, drop=F)
names(dbplaa2)
dbplaa = merge(dbplaa,dbplaa1,all.x=T)
dbplaa = merge(dbplaa,dbplaa2,all.x=T)

if (onpl != "ALL") {
  if (is.null(alldbpla) ){
    dbpla$state<-onpl
    alldbpla_3<-dbpla
    dbplaa$state<-onpl
    alldbpla2_3<-dbplaa
  } else {
    dbpla$state<-onpl
    alldbpla_3<-rbind(alldbpla_3,dbpla)
    dbplaa$state<-onpl
    alldbpla2_3<-rbind(alldbpla2_3,dbplaa)
  }
}
save(dbpla,file=paste0("aggres",fsu,"/dbplPc",onpl,"_agg_",ontype,".RData"))
save(dbplaa,file=paste0("aggres",fsu,"/dbplPc",onpl,"_agg2_",ontype,".RData"))



save.image(paste0("afterTracs_",fsu,"no52_2.RData"))