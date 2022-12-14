rm(list=ls())

# SEM analyses

oncyi<-c(paste0("bio",c(1,12)), "CoastDist","ELEV", "IT200","HFP")
yrepi<-c("maxHT", "maxHTsk","maxDIA","maxDIAsk","sumBA","sumNO","noSPE","sumSEEun")
fpref<-"J1_fws_log"
(tajs<-oncyi[grep("IT",oncyi)])

# all forests
load("all_res/dbplPcALL_agg.RData") # 1deg
load("all_res/dbplPcALL_agg2.RData") # 0.5deg
dbpla_0<-dbpla
dbplaa_0<-dbplaa
# forest types 1 xeric, 2 mesic, 3 hydric
load("all_res/dbplPcALL_agg_1.RData") # 1deg 1type
load("all_res/dbplPcALL_agg2_1.RData") # 0.5deg 1type
dbpla_1 <-dbpla
dbplaa_1<-dbplaa
load("all_res/dbplPcALL_agg_2.RData") # 1deg 2type
load("all_res/dbplPcALL_agg2_2.RData") # 0.5deg 2type
dbpla_2 <-dbpla
dbplaa_2<-dbplaa
load("all_res/dbplPcALL_agg_3.RData") # 1deg 3type
load("all_res/dbplPcALL_agg2_3.RData") # 0.5deg 3type
dbpla_3 <-dbpla
dbplaa_3<-dbplaa

# back to all
dbpla<-dbpla_0
dbplaa<-dbplaa_0

if (zeroone){
  dbpla <- dbplaa
  dbpla_0 <- dbplaa_0  
  dbpla_1 <- dbplaa_1  
  dbpla_2 <- dbplaa_2  
  dbpla_3 <- dbplaa_3  
}

dbpla$ITf200 <- dbpla$ITf200/40
dbpla_1$ITf200 <- dbpla_1$ITf200/40
dbpla_2$ITf200 <- dbpla_2$ITf200/40
dbpla_3$ITf200 <- dbpla_3$ITf200/40
dbpla$ITf100 <- dbpla$ITf200/40
dbpla_1$ITf100 <- dbpla_1$ITf100/40
dbpla_2$ITf100 <- dbpla_2$ITf100/40
dbpla_3$ITf100 <- dbpla_3$ITf100/40

dbpla$ITws200 <- ifelse( (!is.na(dbpla$ITws200)) & (!is.na(dbpla$ITws200)) & (dbpla$ITws200>0),
                         dbpla$ITws200/dbpla$ITfs200,dbpla$ITws200)
dbpla_1$ITws200 <- ifelse( (!is.na(dbpla_1$ITws200)) & (!is.na(dbpla_1$ITws200)) & (dbpla_1$ITws200>0),
                           dbpla_1$ITws200/dbpla_1$ITfs200,dbpla_1$ITws200)
dbpla_2$ITws200 <- ifelse( (!is.na(dbpla_2$ITws200)) & (!is.na(dbpla_2$ITws200)) & (dbpla_2$ITws200>0),
                           dbpla_2$ITws200/dbpla_2$ITfs200,dbpla_2$ITws200)
dbpla_3$ITws200 <- ifelse( (!is.na(dbpla_3$ITws200)) & (!is.na(dbpla_3$ITws200)) & (dbpla_3$ITws200>0),
                           dbpla_3$ITws200/dbpla_3$ITfs200,dbpla_3$ITws200)

te<-prcomp(dbpla[,c("ITf200","ITws200")],scale=T,center = T) #;fviz_pca_biplot(te)
dbpla[["IT200"]] <-te$x[,1]
te<-prcomp(dbpla_1[,c("ITf200","ITws200")],scale=T,center = T) #;fviz_pca_biplot(te)
dbpla_1[["IT200"]] <-te$x[,1]
te<-prcomp(dbpla_2[,c("ITf200","ITws200")],scale=T,center = T) #;fviz_pca_biplot(te)
dbpla_2[["IT200"]] <-te$x[,1]
te<-prcomp(dbpla_3[,c("ITf200","ITws200")],scale=T,center = T) #;fviz_pca_biplot(te)
dbpla_3[["IT200"]] <-te$x[,1]


respf<-c("maxHT", "maxHTsk","maxDIA","maxDIAsk","meanDIAHT","sumBA","sumNO","noSPE","sumSEEun")
respt<-c("maxHeight", "heightSkew", "maxDBH","DBHSkew", "DBHtoHeight","basalArea","NOtrees","NOspecies","NOrecruits")
expf<-c("bio1", "bio12", "ELEV", "HFP","IT200","ITf200","ITws200")
expt<-c("MAT", "MAP", "Elev", "HumanImp","Cyclone","CycloneFreq","CycloneWindSum")

onall<-dbpla[,oncyi]

onall<-scale(log(2+onall))
onallrun<-cbind(onall, log(2+dbpla[,yrepi]), dbpla[,c("LON1","LAT1")] )
#onallrun<-cbind(onall, log(2+dbpla[dbpla$ELEV<1000,yrepi]))
head(onallrun)
#onallCOO<-dbpla[,c("LON1","LAT1",yrepi)]
#type 1
onall_1<-dbpla_1[,oncyi]
onall_1<-scale(log(2+onall_1))
onallrun_1<-cbind(onall_1,log(2+dbpla_1[,yrepi]), dbpla_1[,c("LON1","LAT1")])
head(onallrun_1)
#onallCOO_1<-dbpla_1[,c("LON1","LAT1",yrepi)]
#type 2
onall_2<-dbpla_2[,oncyi]
onall_2<-scale(log(2+onall_2))
onallrun_2<-cbind(onall_2,log(2+dbpla_2[,yrepi]), dbpla_2[,c("LON1","LAT1")])
#onallCOO_2<-dbpla_2[,c("LON1","LAT1",yrepi)]
#type 3
onall_3<-dbpla_3[,oncyi]
onall_3<-scale(log(2+onall_3))
onallrun_3<-cbind(onall_3,log(2+dbpla_3[,yrepi]), dbpla_3[,c("LON1","LAT1")])


modsel<-data.frame(resp=NA, data=NA, n=NA, modno=NA,FisherC=NA,FisherCp=NA,AIC=NA,AICc=NA,W=NA,
                   R2resp=NA,R2Cyclone=NA, R2MAP=NA)
for(i in sort(oncyi)) modsel[[i]]<-NA
for(i in sort(oncyi[c(1,3)]) ) modsel[[paste0("Cyclone<-",i)]]<-NA
for(i in sort(oncyi[c(1,3,4)]) ) modsel[[paste0("MAP<-",i)]]<-NA

library(piecewiseSEM)
library(qgraph)
library(MASS)
library(semEff)
library(car)
#newpairs( onall,c(oncyi[order(oncyi)]),yrepi,tab=F)
(a<-yrepi[7])
# a ~
mmodels<-mmodels2<-mmodels3<-list()
mmodels[[1]] <- paste0(c("bio1","bio12","IT200","HFP"),collapse="+")
mmodels[[2]] <- paste0(c("bio1","bio12","IT200"),collapse="+")
mmodels[[3]] <- paste0(c("bio12","IT200","HFP"),collapse="+")
mmodels[[4]] <- paste0(c("bio12","IT200"),collapse="+")
mmodels[[5]] <- paste0(c("bio1","IT200"),collapse="+")
# bio12 ~
mmodels2[[1]] <- paste0(c("CoastDist","IT200"),collapse="+")
mmodels2[[2]] <- paste0(c("CoastDist"),collapse="+")
mmodels2[[3]] <- paste0(c("CoastDist","IT200","bio1"),collapse="+")
mmodels2[[4]] <- paste0(c("IT200","bio1"),collapse="+")
# IT200 ~
mmodels3[[1]] <- paste0(c("CoastDist"),collapse="+")
mmodels3[[2]] <- paste0(c("CoastDist","bio1"),collapse="+")

modco<-expand.grid(1:length(mmodels),1:length(mmodels2),1:length(mmodels3))
imods<-data.frame(res=NA,on=NA,modno=NA,fisherC=NA,P=NA,AIC=NA,W=NA,itsig=NA)

expnames<-names(modsel)[c(( grep("^IT200",names(modsel)) -4):grep("^IT200",names(modsel)) )]
expNA<-rep(NA,length(expnames))
names(expNA)<-expnames

expnames1<-names(modsel)[ grep("Cyclone<-",names(modsel)) ]
expNA1<-rep(NA,length(expnames1))
names(expNA1)<- gsub("(.*)<-(.*)","\\2",expnames1)

expnames2<-names(modsel)[ grep("MAP<-",names(modsel)) ]
expNA2<-rep(NA,length(expnames2))
names(expNA2)<- gsub("(.*)<-(.*)","\\2",expnames2)

library(spdep)
library(spatialreg)
library(devtools)
library(piecewiseSEM)
#piecewiseSEM:::stdCoefs
# loop for responses
(a<-yrepi[1])

for(doTypes in c("noTypes","justTypes","all")){
  if (doTypes == "noTypes"){
    pdf(paste0(fpref,"semboxesNOe_resp_",tajs,"_",doTypes,".pdf"),8,15)
    par(mfrow=c(4,2),mar=c(1,1,3,1))
    
  }else if (doTypes == "justTypes"){
    pdf(paste0(fpref,"semboxesNOe_resp_",tajs,"_",doTypes,".pdf"),15,30)
    par(mfrow=c(8,3),mar=c(1,1,3,1))
  }else {
    pdf(paste0(fpref,"semboxesNOe_resp_",tajs,"_",doTypes,".pdf"),20,30)
    par(mfrow=c(8,4),mar=c(1,1,3,1))
  }
  
  #for(a in yrepi[c(1,3)] ) {
  for(a in yrepi[c(1,3,7,6,8 ,5,2,4)] ) {
    onall <- onallrun[ (!is.na(onallrun[,a])) & (!is.infinite(onallrun[,a])),]# & (!is.na(onallrun$ELEV)), ]
    onallcoo <- onall #COO[ (!is.na(onallrun[,a])) & (!is.infinite(onallrun[,a])), ]
    
    # autoregresive spatial SAR
    coordinates(onallcoo)<- ~ LON1 + LAT1
    Wv=autoremv[1]; W="x"
    if (autorem) {
      W<-tryCatch(
        nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
        , error=function(e) return("empty"))#neighbours within 150 km
      if( class(W)[1]!="listw" ) {
        Wv=autoremv[2]
        W<-tryCatch(nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
                    , error=function(e) return("empty"))
      }
      if( class(W)[1]!="listw" ) {
        Wv=autoremv[3]
        W<-tryCatch(nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
                    , error=function(e) return("empty"))          
      }
      Wv0<-Wv # start for loop of sems
      W0<- W
    }
    #plot(W, coordinates(onallcoo)) #visualize network of neighbors
    #moran.test( onall[,a],listw=W)
    # length(onall[,a]); rownames(W$neighbours)
    
    fia<-c(); fiaw<-c(); mle=1
    print(paste("......Start of",a,"with Wv =",Wv))
    for(mle in 1:nrow(modco) ) {
      print(paste("......model no =",mle," on",a))
      (afo<-paste0(a,"~",mmodels[[ modco[mle,1] ]] ) )
      (afo2<-paste0("bio12~",mmodels2[[ modco[mle,2] ]] ) )
      (afo3<-paste0("IT200~",mmodels3[[ modco[mle,3] ]] ) )
      if (! autorem) { #lm
        modell<-list(lm( as.formula(afo), onall), lm(as.formula(afo2),onall), lm(as.formula(afo3),onall))
      } else { #errorsarlm
        Wv<-Wv0
        W<-nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
        modell<-list(errorsarlm( as.formula(afo), data=onall, listw = W), 
                     errorsarlm( as.formula(afo2),data=onall, listw = W), 
                     errorsarlm( as.formula(afo3),data=onall, listw = W))
      }
      model <- psem( modell[[1]], modell[[2]], modell[[3]])
      (moaic<-unlist(summary(model)$IC[1]))
      print(moaic)
      # test of distances for spatial autocorrelation - loop
      if (autorem) {
        for(wdi in autoremv[(grep(Wv0,autoremv)+1):length(autoremv)] ) if (!is.na(wdi)){
          Wvw<-wdi
          Ww<-nb2listw(dnearneigh(coordinates(onallcoo), 0, Wvw, longlat = TRUE), style="W")
          modellw<-list(errorsarlm( as.formula(afo), data=onall, listw = Ww), 
                        errorsarlm( as.formula(afo2),data=onall, listw = Ww), 
                        errorsarlm( as.formula(afo3),data=onall, listw = Ww))
          modelw <- psem( modellw[[1]], modellw[[2]], modellw[[3]])
          (moaicw<-unlist(summary(modelw)$IC[1]))
          print(paste(wdi, moaicw))
          if ( moaicw < moaic ) { # new W is better
            print(paste(moaicw,"<",moaic," for wdi",wdi))
            Wv<-wdi
            W<-nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
            modell<-list(errorsarlm( as.formula(afo), data=onall, listw = W), 
                         errorsarlm( as.formula(afo2),data=onall, listw = W), 
                         errorsarlm( as.formula(afo3),data=onall, listw = W))
            model <- psem( modell[[1]], modell[[2]], modell[[3]])
            moaic<-moaicw
          }
        } 
      }
      #summary(model)
      mlefish <-fisherC(model)
      print( mlefish  )
      (x0<-coefs(model)); x1<-x0[ (x0$Response ==a ),]
      (itsig=x1[x1$Predictor == "IT200",]$P.Value)
      imods<-rbind(imods, c(a,"all",mle,unname(mlefish[1]),unname(mlefish[3]),unname(moaic),Wv,itsig) )
      (fia<-c(fia, mlefish$P.Value ))
      #(fia<-c(fia, unname(moaic) ))
      (fiaw<-c(fiaw, Wv ))
    }
    i22<-imods[ !is.na(imods$res) & (imods$res == a) & (imods$on == "all") & (imods$itsig < 0.05), ]
    bestaicmodel <- i22[which.max(i22$P),]$modno
    #bestaicmodel <- i22[which.min(i22$AIC),]$modno
    if ( length(bestaicmodel) == 0 ) {
      bestaicmodel<-which.max(fia)
    } else {
      if ( i22[ i22$modno== bestaicmodel,]$P ==0 ) {
        bestaicmodel<-which.max(fia)
      }
    }
    #bestaicmodel<-which.max(fia)
    (afo<-paste0(a,"~",mmodels[[ modco[bestaicmodel,1] ]] ) )
    (afo2<-paste0("bio12~",mmodels2[[ modco[bestaicmodel,2] ]] ) )
    (afo3<-paste0("IT200~",mmodels3[[ modco[bestaicmodel,3] ]] ) )
    
    if (! autorem) {
      modell<-list(lm( as.formula(afo), onall), lm(as.formula(afo2),onall), lm(as.formula(afo3),onall))
    } else {
      Wv<-fiaw[bestaicmodel];
      W<-nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
      print (paste("..Best Wv for",a,"is", Wv))
      modell<-list(errorsarlm( as.formula(afo), data=onall, listw = W), 
                   errorsarlm( as.formula(afo2),data=onall, listw = W), 
                   errorsarlm( as.formula(afo3),data=onall, listw = W))
    }
    
    model <- psem( modell[[1]], modell[[2]], modell[[3]])
    #fisherC(model) 
    #summary(model);dSep(model)
    #print(lapply(modell[[1]], vif))
    (msu<-summary(model)); 
    (x0<-coefs(model))
    x0[,ncol(x0)] <- ifelse (x0[,ncol(x0)] == "", ifelse(x0[,ncol(x0)-2]<=0.1,"'",""), x0[,ncol(x0)] )
    
    #resp
    x1<-x0[ (x0$Response ==a ),]
    x1vals<-unname(apply( x1[order(x1$Predictor),ncol(x1)-c(1,0)],1,paste0,collapse=""))
    names(x1vals)<-c(x1$Predictor[order(x1$Predictor)])
    expval<-expNA
    expval[ names(x1vals)]<-x1vals
    #cyclone
    x1<-x0[ (x0$Response =="IT200" ),]
    x1vals<-unname(apply( x1[order(x1$Predictor),ncol(x1)-c(1,0)],1,paste0,collapse=""))
    names(x1vals)<-c(x1$Predictor[order(x1$Predictor)])
    expval1<-expNA1
    expval1[ names(x1vals)]<-x1vals
    #MAP
    x1<-x0[ (x0$Response =="bio12" ),]
    x1vals<-unname(apply( x1[order(x1$Predictor),ncol(x1)-c(1,0)],1,paste0,collapse=""))
    names(x1vals)<-c(x1$Predictor[order(x1$Predictor)])
    expval2<-expNA2
    expval2[ names(x1vals)]<-x1vals
    
    modsel<-rbind(modsel, c(nname(a),"all", nrow(onall), 
                            unname(bestaicmodel),
                            unname(msu$Cstat[c(1,3)]),
                            unname(msu$IC[1:2]),Wv,
                            msu$R2[msu$R2$Response ==a,"R.squared"],
                            msu$R2[msu$R2$Response =="IT200","R.squared"],
                            msu$R2[msu$R2$Response =="bio12","R.squared"],
                            unname(expval),unname(expval1),unname(expval2)
    ) )
    
    y<-as.matrix(x0[,c(2:1,ncol(x0)-c(1) )])
    x0$Std.Estimate <- round( x0$Std.Estimate,3)
    save.image(paste0(fprefo,a,"_type_all.RData"))
    if (doTypes %in% c("noTypes","all") ){
      nicesem(y,x0,ns=ifelse(msu$Cstat[c(3)] <0.05,TRUE, FALSE)); 
      mtext(paste0("FisherC=",msu$Cstat[c(1)]," p=",msu$Cstat[c(3)]," "),3,1,adj=1)
      mtext( nname(a),3,1,adj=0,cex=2)
      if (doTypes %in% c("all") ) mtext( "all",3,-1,cex=2)
      R2val<-msu$R2[msu$R2$Response ==a,"R.squared"]
      #text( 1,0.3, bquote(R^2 == .(R2val)),cex=1.2)
      an<-nname(a)
      ap<-nname("bio12")
      ac<-nname("IT200")
      mtext(bquote(.(an) ~ R^2 == .(R2val)),3,-0.5,adj=1)
      R2val<-msu$R2[msu$R2$Response =="bio12","R.squared"]
      mtext(bquote(.(ap)~R^2 == .(R2val)),3,-2,adj=1)
      R2val<-msu$R2[msu$R2$Response =="IT200","R.squared"]
      mtext(bquote(.(ac)~R^2 == .(R2val)),3,-3.5,adj=1)
    }
    # forest types
    ty=1; for(ty in 1:3){
      if (ty == 1) {
        onall<-onallrun_1[ !is.infinite(onallrun_1[,a]) & !is.na(onallrun_1[,a]),]# & (!is.na(onallrun_1$ELEV)) , ]
      } else if (ty == 2) {
        onall<-onallrun_2[ !is.infinite(onallrun_2[,a]) & !is.na(onallrun_2[,a]),]# & (!is.na(onallrun_2$ELEV)), ]
      } else {
        onall<-onallrun_3[ !is.infinite(onallrun_3[,a]) & !is.na(onallrun_3[,a]),]# & (!is.na(onallrun_3$ELEV)), ]
      }
      
      onallcoo <- onall
      coordinates(onallcoo)<- ~ LON1 + LAT1
      Wv=autoremv[1]; W="x"
      if (autorem) {
        W<-tryCatch(
          nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
          , error=function(e) return("empty"))#neighbours within 150 km
        if( class(W)[1]!="listw" ) {
          Wv=autoremv[2]
          W<-tryCatch(
            nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
            , error=function(e) return("empty"))
        }
        if( class(W)[1]!="listw" ) {
          Wv=autoremv[3]
          W<-tryCatch(nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
                      , error=function(e) return("empty"))
        }
        Wv0<-Wv # start for loop of sems
        W0<- W  
      }
      
      fia<-c();fiaw<-c()
      print(paste("......Start of",a,"with Wv=",Wv,"type=",ty))
      for(mle in 1:nrow(modco) ) {
        print(paste("......model no =",mle," on",a))
        (afo<-paste0(a,"~",mmodels[[ modco[mle,1] ]] ) )
        (afo2<-paste0("bio12~",mmodels2[[ modco[mle,2] ]] ) )
        (afo3<-paste0("IT200~",mmodels3[[ modco[mle,3] ]] ) )
        
        if (!autorem){
          modell<-list(lm( as.formula(afo), onall), lm(as.formula(afo2),onall), lm(as.formula(afo3),onall))
        } else {
          Wv <- Wv0
          W <- nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
          modell<-list(errorsarlm( as.formula(afo), data=onall, listw = W), 
                       errorsarlm( as.formula(afo2),data=onall, listw = W), 
                       errorsarlm( as.formula(afo3),data=onall, listw = W))
          
        }
        model <- psem( modell[[1]], modell[[2]], modell[[3]])
        (moaic<-summary(model)$IC[1])
        if (autorem){ for(wdi in autoremv[grep(Wv0,autoremv)+1:length(autoremv)]) if (!is.na(wdi)) {
          Wvw<-wdi
          Ww<-nb2listw(dnearneigh(coordinates(onallcoo), 0, Wvw, longlat = TRUE), style="W")
          modellw<-list(errorsarlm( as.formula(afo), data=onall, listw = Ww), 
                        errorsarlm( as.formula(afo2),data=onall, listw = Ww), 
                        errorsarlm( as.formula(afo3),data=onall, listw = Ww))
          modelw <- psem( modellw[[1]], modellw[[2]], modellw[[3]])
          (moaicw<-summary(modelw)$IC[1])
          print(paste(wdi, moaicw))
          if ( moaicw < moaic) { # new W is better
            Wv<-wdi
            W<-nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
            modell<-list(errorsarlm( as.formula(afo), data=onall, listw = W), 
                         errorsarlm( as.formula(afo2),data=onall, listw = W), 
                         errorsarlm( as.formula(afo3),data=onall, listw = W))
            model <- psem( modell[[1]], modell[[2]], modell[[3]])
            moaic<-moaicw
            #model<-modelw
          }
        } }
        mlefish <-fisherC(model)
        print( mlefish  )
        (x0<-coefs(model)); x1<-x0[ (x0$Response ==a ),]
        (itsig=x1[x1$Predictor == "IT200",]$P.Value)
        imods<-rbind(imods, c(a,ty,mle,unname(mlefish[1]),unname(mlefish[3]),unname(moaic),Wv,itsig) )
        (fia<-c(fia, mlefish$P.Value ))
        #(fia<-c(fia, unname(moaic)))
        (fiaw<-c(fiaw, Wv ))
        #(fia<-c(fia, fisherC(model)$P.Value ))
      }
      i22<-imods[ !is.na(imods$res) & (imods$res == a) & (imods$on == ty) & (imods$itsig < 0.05), ]
      bestaicmodel <- i22[which.max(i22$P),]$modno
      #bestaicmodel <- i22[which.min(i22$AIC),]$modno
      if ( length(bestaicmodel) ==0 ) {
        bestaicmodel<-which.max(fia)
      } else {
        if ( i22[ i22$modno== bestaicmodel,]$P ==0 ) {
          bestaicmodel<-which.max(fia)
        }
      }  
      (afo<-paste0(a,"~",mmodels[[ modco[bestaicmodel,1] ]] ) )
      (afo2<-paste0("bio12~",mmodels2[[ modco[bestaicmodel,2] ]] ) )
      (afo3<-paste0("IT200~",mmodels3[[ modco[bestaicmodel,3] ]] ) )
      
      if (!autorem){
        modell<-list(lm( as.formula(afo), onall), lm(as.formula(afo2),onall), lm(as.formula(afo3),onall))
      } else {
        Wv<-fiaw[bestaicmodel];
        W<-nb2listw(dnearneigh(coordinates(onallcoo), 0, Wv, longlat = TRUE), style="W")
        print (paste("..Best Wv for",a,"is", Wv))
        modell<-list(errorsarlm( as.formula(afo), data=onall, listw = W), 
                     errorsarlm( as.formula(afo2),data=onall, listw = W), 
                     errorsarlm( as.formula(afo3),data=onall, listw = W))
      }
      model <- psem( modell[[1]], modell[[2]], modell[[3]])
      
      #print(lapply(modell, vif))
      #dSep(model)
      (msu<-summary(model)); (x0<-coefs(model))
      x0[,ncol(x0)] <- ifelse (x0[,ncol(x0)] == "", ifelse(x0[,ncol(x0)-2]<=0.1,"'",""), x0[,ncol(x0)] )
      
      #resp
      x1<-x0[ (x0$Response ==a ),]
      x1vals<-unname(apply( x1[order(x1$Predictor),ncol(x1)-c(1,0)],1,paste0,collapse=""))
      names(x1vals)<-c(x1$Predictor[order(x1$Predictor)])
      expval<-expNA
      expval[ names(x1vals)]<-x1vals
      #cyclone
      x1<-x0[ (x0$Response =="IT200" ),]
      x1vals<-unname(apply( x1[order(x1$Predictor),ncol(x1)-c(1,0)],1,paste0,collapse=""))
      names(x1vals)<-c(x1$Predictor[order(x1$Predictor)])
      expval1<-expNA1
      expval1[ names(x1vals)]<-x1vals
      #MAP
      x1<-x0[ (x0$Response =="bio12" ),]
      x1vals<-unname(apply( x1[order(x1$Predictor),ncol(x1)-c(1,0)],1,paste0,collapse=""))
      names(x1vals)<-c(x1$Predictor[order(x1$Predictor)])
      expval2<-expNA2
      expval2[ names(x1vals)]<-x1vals
      
      modsel<-rbind(modsel, c(nname(a), c("xeric","mesic","hydric")[ty], nrow(onall), 
                              unname(bestaicmodel),
                              unname(msu$Cstat[c(1,3)]),
                              unname(msu$IC[1:2]),Wv,
                              msu$R2[msu$R2$Response ==a,"R.squared"],
                              msu$R2[msu$R2$Response =="IT200","R.squared"],msu$R2[msu$R2$Response =="bio12","R.squared"],
                              unname(expval),unname(expval1),unname(expval2)
      ) )
      save.image(paste0(fprefo,a,"_type_",ty,".RData"))
      if (doTypes %in% c("justTypes","all") ){
        y<-as.matrix(x0[,c(2:1,ncol(x0)-c(1) )])
        x0$Std.Estimate <- round( x0$Std.Estimate,3)
        nicesem(y,x0,ns=ifelse(msu$Cstat[c(3)] <0.05,TRUE, FALSE));
        mtext(paste0("FisherC=",msu$Cstat[c(1)]," p=",msu$Cstat[c(3)]," "),3,1,adj=1)
        mtext(c("xeric","mesic","hydric")[ty],3,1,cex=1.5)
        R2val<-msu$R2[msu$R2$Response ==a,"R.squared"]
        mtext(bquote(.(an) ~ R^2 == .(R2val)),3,-0.5,adj=1)
        R2val<-msu$R2[msu$R2$Response =="bio12","R.squared"]
        mtext(bquote(.(ap) ~ R^2 == .(R2val)),3,-2,adj=1)
        R2val<-msu$R2[msu$R2$Response =="IT200","R.squared"]
        mtext(bquote(.(ac) ~ R^2 == .(R2val)),3,-3.5,adj=1)
      }
    }
    
    foo <- nname( names(modsel))
    names(modsel)<-foo
    foo2 <-nname(modsel$on)
    modsel$on <- foo2
    write.csv(modsel[-c(1),-c(grep("^CoastDist",names(modsel)))],paste0(fpref,"semboxesNOe_resp_",tajs,".csv"),row.names = F)
    
  }  
  dev.off()
  write.csv(imods[-c(1),],paste0(fpref,"_indivmods.csv"),row.names=F)
}
