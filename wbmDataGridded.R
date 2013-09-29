#SHAPEFILE#################################################################################################
pj4s<-CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")
eaf<- readShapeSpatial("GIS/Hydrology/Study_limits.shp",proj4string=pj4s)
ind<- readShapeSpatial("GIS/EastAfrica/indejeWGS.shp",proj4string=pj4s)
wat<- readShapeSpatial("GIS/EastAfrica/Water_Bodeis.shp",proj4string=pj4s)
ctry<-readShapeSpatial("GIS/EastAfrica/East_Africa_ctry.shp",proj4string=pj4s)
getinfo.shape("GIS/EastAfrica/eafcsel.shp")
#RASTERS#####################################################################################################
dem<-raster("GIS\\Hydrology\\DEM\\ea_dem30");save(dem,file="dem")
fac<-raster("GIS\\Hydrology\\DEM\\ea_acc30");save(fac,file="fac")
dir<-raster("GIS\\Hydrology\\DEM\\ea_dir30");save(dir,file="dir")
fx.mpz(dem,main="East Africa elevation",xlim=c(29,42),ylim=c(-15,6.5)) #plot all maps
#get the catchment  bounds
############################################################################################################
(box<-basin@bbox)
#Area of basin
(area<-fx.areasqm(basin)) #
(XY<-coordinates(basin)) #centroid
#RAIN DATA x,y, range #####################################################################################
#latitudes & longitudes CRU
#y<-seq(-89.75,90,Dxy); x<-seq(-179.75,180,Dxy) #africa
y<-seq(-12,2,Dxy); x<-seq(25,45,Dxy) #east Africa
(xmn<-x[which(abs(x-box["x","min"])==min(abs(x-box["x","min"])))]-diff(x)[1]) #minus resolution
(xmx<-x[which(abs(x-box["x","max"])==min(abs(x-box["x","max"])))]+diff(x)[1]) #add   resolution
(ymn<-y[which(abs(y-box["y","min"])==min(abs(y-box["y","min"])))]-diff(x)[1])
(ymx<-y[which(abs(y-box["y","max"])==min(abs(y-box["y","max"])))]+diff(x)[1])
#COLOR################################################################################################
#(basin<-unionSpatialPolygons(eaf,eaf@data$Id) );cname="east_Africa" #East Africa
kagera<-TRUE
#StnID|Sname|District|Northing|Easting|Elevation
#COLORS
kol1<-rgb( colorRamp(c("azure", "blue"))(seq(0, 1, length = 100)), max = 255)
kol2<-rgb( colorRamp(c("yellow","red"))(seq(0, 1, length = 100)), max = 255)
acol=kol1#colorRampPalette(c("orange","blue"))(3000)
acol<-color.scale(0:3000,col=acol)[-1];abrk<-seq(-1,2999,1);
lcol<-((c("white","blue")));lbrk<-c(-1,0.5,1)
zcol<-color.scale(0:6000,col=terrain.colors(6000))[-1];zbrk<-seq(-1,5999,1)
tcol<-color.scale(0:30,colorRampPalette(c("blue","cadetblue","yellow","red"))(60))[-1];tbrk<-seq(-1,29,1)
fmt<-"CRU"
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#mean elevation in Basin
library(rgeos)
basin <- gUnaryUnion(basin) #union the polygons in basin
try(z<- (raster::extract(dem,basin))[[1]],silent=TRUE)
Z<- mean(z,na.rm=TRUE)#(lapply(elevs,FUN=mean)))
#pdf(paste(cname,"hypso.pdf"))
hz<-fx.hypso(z)
##dev.off()
#Plotting Limits###########################################################################################
xL<-c(-18,50); yL<-c(-40,40) #Africa Plotting limits at CRU data extraction
xL<-c(25,45); yL<-c(-12,2) #Africa Plotting limits at CRU data extraction
    
#Use Basin Bounds to select from CRU data####################################################################################  
#SELECTION#africa#lonx<-(-18.25:60);laty<-(-40.25:40)
(lonx<-seq(xmn,xmx,diff(x)[1])); 
(laty<-seq(ymn,ymx,diff(x)[1]))
#WORGING DIRECTORY
#############################################################################################################
(outDir<-paste("C:/DATA/STATIONS/wbm_",cname,"/Calib",sep=""))
dir.create(outDir)

#DATA #cru
#CRU retreieve global & extract ctahcment based on the shapefile (TRY GHCN)
##PRECIPITATION###########################################################################################  
library(audio)#wait
# if(file.exists(paste(cname,"prec.Rdata",sep="."))){
#   Rain<-get(load(paste(cname,"prec.Rdata",sep=".")))
#   } else {
  fil.prc<-"CRU/badc/cru.prec.Rdata"
  wait(system.time(Rain<-cru.retrieve(fil=fil.prc,lonx=lonx,laty=laty,xL=xL,yL=yL,MULTI=0.1,col=kol1)),60)
  save(Rain,file=paste(outDir,"/",cname,".prec.Rdata",sep=""))
# }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#pdf(paste(cname,"Prec.PE.elev.lakes.pdf",sep="."))
#par(mfrow=c(2,2),mar=c(4,4,4,4))
for(i in 1){
  #image.plot(lonx,laty,Rain[,,i],main=paste("Precipitation",month.abb[i]),asp=1,col=acol)
  rain<-apply(Rain,c(1,2),mean)*12
  image.plot(lonx,laty,rain,main=paste("Mean Annual Precipitation 61-90"),asp=1,col=acol)
  plot(basin,add=TRUE)
  plot(wat,add=TRUE,col="transparent",border="blue")
  };
#title("PREC",outer=TRUE);box();#par(mfrow=c(1,1))
#TEMPERATURE##################################################################################################  
# if(file.exists(paste(cname,"tmp.Rdata",sep="."))){
#   Temp<-get(load(paste(cname,"tmp.Rdata",sep=".")))
#   } else{
  fil.tmp<-"CRU/badc/cru.tmp.Rdata"
  wait(system.time(Temp<-cru.retrieve(fil=fil.tmp,lonx=lonx,laty=laty,xL=xL,yL=yL,MULTI=0.1,col=kol2)),60)
  save(Temp,file=paste(outDir,"/",cname,".tmp.Rdata",sep=""))
# }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(i in 1){
  #image.plot(lonx,laty,Temp[,,i],main=paste("Temperature",month.abb[i]),asp=1,col=tcol,breaks=tbrk)
  temp<-apply(Temp,c(1,2),mean)
  image.plot(lonx,laty,temp,main=paste("Mean monthly temperature 61-90"),asp=1,col=tcol)
  plot(basin,add=TRUE)
  plot(wat,add=TRUE,col="transparent",border="blue")
}
#title("TEMP",outer=TRUE);box();#par(mfrow=c(1,1))

#Shape to raster to Matrix##########################################################################################  
#Shape to raster at resolution diff(x) within xmn-xmx & ymn-ymx
  # Create a empty raster based on coordinates
  ras <- raster(ncols=length(lonx), nrows=length(laty))
  limRas<-c(min(lonx)-diff(lonx)[1]/2,max(lonx)+diff(lonx)[1]/2,min(laty)-diff(laty)[1]/2,max(laty)+diff(laty)[1]/2);#raster limits
  (limRas<-matrix(limRas,2,2,byrow=TRUE))
  (colnames(limRas)<-c("min","max"))
  (rownames(limRas)<-c("x","y"))
  (basin@bbox<-limRas)
  (extent(ras)<-extent(limRas))#raster extent
  R<-Rain[,,1]#raster data based on Rain
  ras[]<-1#t(R[,ncol(R):1]);
  bas<-rasterize(basin,ras)
  bas<-as.matrix(bas)#basin
  bas<-t(bas[nrow(bas):1,])
  #Matrix to remove outer cells
  bas[!is.na(bas)]<-1; 
  print(bas)#
  (ncells<-length(which(!is.na(bas)))) #cells in catchment
  #CELL AREA
  (cell.area<-sum(area)/ncells) #m2
  #LAKE CELLS
  Lake<-rasterize(wat,ras)#lakes
  Lake<-as.matrix(Lake)
  Lake[!is.na(Lake)]<-1
  Lake<-t(Lake[nrow(Lake):1,])
  Lake<-Lake*bas
  Lake[is.na(Lake)]<-0
 
  #Z-raster
  El<-crop(dem,ras)
  Z<-aggregate(El,fact=(dim(El)/dim(ras))[1])#Z<-resample(Z,ras)
  Z<-as.matrix(Z)
  Z<-t(Z[nrow(Z):1,])
  Z<-Z*bas
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fx.mpz(El,bas=basin)
image.plot(lonx,laty,Lake,main="Lake Cells",asp=1,col=(lcol),breaks=lbrk);  plot(basin,add=TRUE);plot(wat,add=TRUE,border=4);grid()
##dev.off()
#Plots######################################################################################################  
  #plots
    Y<-as.numeric(basin@bbox[2,])
    X<-as.numeric(basin@bbox[1,])
    lx<-seq(min(lonx)-diff(lonx)[1]/2,max(lonx)+diff(lonx)[1]/2,diff(lonx)[1])
    ly<-seq(min(laty)-diff(laty)[1]/2,max(laty)+diff(laty)[1]/2,diff(laty)[1])
    abl<-function(xxx,lx,ly){
      if(xxx==2){abline(v=lx,h=ly,lty=2,col="grey");axis(1,lx,lx);axis(2,ly,ly);box()}
      if(xxx==1)plot(expand.grid(lx,ly),xaxt="n",yaxt="n",ylab="lat",xlab="lon",pch="",asp=1)
    }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#pdf(paste(cname,"test.extraction.pdf",sep="."))
#par(mfrow=c(2,2),mar=c(4,4,4,4))
    abl(1,lx,ly);image.plot(lonx,laty,R,col=acol,ylim=Y,xlim=X,asp=1,add=TRUE);plot(basin,add=TRUE);abl(2,lx,ly);title("MATRIX")
    abl(1,lx,ly);plot(ras,col=acol,ylim=Y,xlim=X,asp=1,add=TRUE);plot(basin,add=TRUE);abl(2,lx,ly);title("RAS")
    abl(1,lx,ly);image.plot(lonx,laty,bas*R,col=acol,ylim=Y,xlim=X,asp=1,add=TRUE);plot(basin,add=TRUE);abl(2,lx,ly);title("bas:Ras")
    abl(1,lx,ly);image.plot(lonx,laty,bas*R,add=TRUE,col=acol);plot(basin,add=TRUE);abl(2,lx,ly);title("bas:Matrix")
##dev.off()

#RUNOFF#######################################################################################################################
#FLOW DATA
# 	$sname
# [1] "mtera.flow"
# # $flo
#         Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec
# 1940  167.2  212.2  305.0  312.7  171.5   71.1   34.0   18.6   10.9    6.7    4.8   60.6
#flofile<-"C:\\DATA\\HYDROMET\\FLOW\\EAPMP\\mtera.flo.rData" # TS 
Qobs<-get(load(flofile))$flo
Qobs<-zoo(Qobs,as.yearmon(time(Qobs)))
par(mar=c(5,2,2,2))
  plot(Qobs,xaxt="n",ylab="m3/s",xlab="",col=4);
  axis(1,(time(Qobs))[seq(1,length(Qobs),24)], (time(Qobs)[seq(1,length(Qobs),24)]),las=2,cex=0.2)
  abline(v=(time(Qobs))[seq(1,length(Qobs),24)],lty=2,col="grey",h=pretty(Qobs))
##dev.off()
#TIME#########################################################################################
#Combine Time periods & remove area outside catchment
TimeMet<-as.yearmon(time(ts(0,start=c(1901,1),end=c(2006,12),f=12))) #time of Rain/Temp records
TimeQ<-time(Qobs) #time of flow records
mIDX<-match(time(Qobs),TimeMet)
prc<-Rain[,,mIDX] #prec
tmp<-Temp[,,mIDX] # tmp

#Removing.array for outside cells###########################################################################################
library(abind)
rem.arr<-bas #array to exlude values in array of rsults that are ou of catchmnt
for(ii in 1:(dim(prc)[3]-1)){rem.arr<-abind(rem.arr,bas,along=3)}

#latdeg#####################################################################
(londeg<-matrix(lonx,dim(Rain)[1],dim(Rain)[2]))
(latdeg<-matrix(laty,dim(Rain)[1],dim(Rain)[2],byrow=TRUE))
Lat.arr<-latdeg
for(ii in 1:(dim(prc)[3]-1)){Lat.arr<-abind(latdeg,Lat.arr,along=3)}
#FC<-array(300,dim(prc)[1:2]) # estimate FC
################################################################################################################
(Years<-year(TimeQ))
(Months<-month(TimeQ))
Time.arr<-Lat.arr
for(i in 1:length(Months)){
  for(j in 1:nrow(latdeg)){
    for(k in 1:ncol(latdeg)){
       Time.arr[j,k,i]<-as.yearmon(TimeQ[i])}     
    }
  }
################################################################################################################
#Evaporation calc
library(SPEI) #thornthwaite
pet<-rem.arr
for(i in 1:nrow(latdeg)){
  for(j in 1:ncol(latdeg)){
    lat<-latdeg[i,j]
    tp<-ts(tmp[i,j,],start=c(Years[1],Months[1]),f=12)
    Pre<-ts(prc[i,j,],start=c(Years[1],Months[1]),f=12)
    pet[i,j,]<-SPEI::thornthwaite(Tave=tp,lat=lat)
    if(is.na(bas[i,j])){pet[i,j,]<-NA;}  
  }
}
################################################################################################################
################################################################################################################
################################################################################################################
#   	///////////DATA STRUCTURE/////////////////@  = QsimC 
# 		//		"nm" "nc" "pe" "wbm"
# 		//		 nm   nc   pe   wbm
# 		//       prec.tmpr.epot.d.h
# 		// C0 t0 ID1 =i*nm+t=to
# 		// C0 t1 ID2
# 		// ..........
# 		// c0 nm (tIMESTEPS).................
# 		///////////////////////////////////////////
#  		///   ir,ic,z,Lake,cellarea,latdeg 
# 		// c01
# 		// c02
# 		////////////////////////////////////////////
#  		//Qobs
# 		//to
# 		//t1

setwd(outDir)
#setwd("C:\\DATA\\STATIONS\\wbm_Rusumo\\Calib\\matlab")
(nm<-length(TimeQ))
(nc<-length(bas[!is.na(bas)]))
petype=3
wbmtype=3
init=cbind(nm,nc,petype,wbmtype)
write.table(init,file="Obs.ts")
write.table(init,file="init.txt",sep="\t",quote=FALSE)
#.............................................................................................................#
fx.array2mts<-function(Array,rem.arr){#turmns array in series by column, eliinates NA cells in rem.arr
  mts<-t(apply(Array,3,c)); mts<-mts[1:length(mts)]
  rem<-t(apply(rem.arr,3,c)); rem<-rem[1:length(rem)]
  valid<-which(!is.na(rem))
  mts<-mts[valid]
}
(prec<-fx.array2mts(prc,rem.arr))
(tmpr<-fx.array2mts(tmp,rem.arr))
(epot<-fx.array2mts(pet,rem.arr))
Dd<-fx.Monthdays.Daylight(Lat.arr,Time.arr)
(d<-array(Dd$d,dim(rem.arr)))
(d<-fx.array2mts(d,rem.arr))
(h<-Dd$D)
(h<-fx.array2mts(h,rem.arr))
inp<-cbind(prec,tmpr,epot,d,h)
inp<-apply(inp,2,function(x)round(x,2))
write.table(inp,file="Obs.ts",row.names=FALSE,append=TRUE,quote=FALSE,sep="\t")
write.table(inp,file="met.txt",row.names=FALSE,quote=FALSE,sep="\t")
#.............................................................................................................#
(iric=which(!is.na(bas),arr.ind=TRUE)) #valid cells indices
(ir<-iric[,1]);
(ic<-iric[,2])
(z<-Z[iric])
(lake=Lake[iric])
(cellarea<-rep(cell.area,nc))
(Lat<-latdeg[iric])
(Lon<-londeg[iric])
fxd<-cbind(ir,ic,z,lake,cellarea,Lat)
fxd<-apply(fxd,2,function(x)round(x,2))
write.table(fxd,file="Obs.ts",row.names=FALSE,append=TRUE,quote=FALSE,sep="\t")
write.table(fxd,file="fxd.txt",row.names=FALSE,quote=FALSE,sep="\t")
#.............................................................................................................#
write.table(as.data.frame(Qobs),file="Obs.ts",row.names=FALSE,append=TRUE,quote=FALSE,sep="\t")
Time<-as.character(as.Date(TimeQ,format="%d/%mm%yyy"))
Qobs<-data.frame(Time,Qobs)
write.table(as.data.frame(Qobs),file="Qobs.txt",row.names=FALSE,quote=FALSE,sep="\t")
#.............................................................................................................#
write.table(latdeg,file="lat.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(londeg,file="lon.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
  