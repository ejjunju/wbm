#Find and change into the script directory (when sourced)
(script.dir <- dirname(sys.frame(1)$ofile))
print(paste("script.dir <-",script.dir))
setwd(script.dir)

library(maptools)
gpclibPermit()
library(lattice)
library(rgdal)
library(topmodel)
library(akima)
library(clim.pact)
library(fields)
library(animation) #saveMovie
library(grDevices)#coloRampPalette
library(fields)
library(raster)
library(hydroTSM)
library(RFOC)
library(ggplot2)
library(reshape)#melt
library(gdata);#gpclibPermit()
library(zoo);library(tis);library(RSEIS) ;library(abind)#to Julian
library(fields);library(maptools);library(rgdal)
#Objects
kol<-colorRampPalette(c("aliceblue","cadetblue1","cornflowerblue","blue"))
load("myw.rda") #map data
load("world_countries.rda")

#
fx.addshp<-function(riv=riv,wat=wat,hep=hep,ctry=ctry,bas1=bas1,add=FALSE){
	#Add plots of read shapefiles for river, waterbodies, pouints(hep),country and basins.
	#can be the same shape 5 times
	#if(add==FALSE){sp::plot(riv,col="lightblue")} else{sp::plot(riv,add=TRUE,col="lightblue")}
	sp::plot(wat,add=TRUE,border="lightblue",col="lightblue")
	sp::plot(riv,add=TRUE,col="lightblue")
	sp::plot(ctry,add=TRUE,border="black")
	sp::plot(bas1,add=TRUE,border="green",lwd=3)
	plot(hep,add=TRUE,pch=16,cex=1.5,col=2)} 

#Compute trends on an array (monthly data)
fx.trend.array<-function(ax,y1=c(1961,1),freq=12,method=c("MK","SMK","zyp"),monx=12){
 trd<-sig<-array(dim=c(dim(ax)[1],dim(ax)[2]))
 if(method[1]=="MK"||method[1]=="SMK"){
    require(Kendall)
    for(i in 1:dim(sig)[1]){
      for(j in 1:dim(sig)[2]){
        x<-ts(ax[i,j,],y1,f=12)
        if(all(!is.na(x))){
          if(method[1]=="Kendall"){MK<-(MannKendall(x))} else {MK<-(SeasonalMannKendall(x))}
          trd[i,j]<-as.numeric(MK$tau)
          sig[i,j]<-as.numeric(MK$sl)
          if(sig[i,j]<=0.05){sig[i,j]<-0}else {sig[i,j]<-1}
          }
        }
      }
    }
  else{
    require(zyp)
    for(i in 1:dim(sig)[1]){
      for(j in 1:dim(sig)[2]){
        x<-ts(ax[i,j,],y1,f=freq)
        if(all(!is.na(x))){
          MK<-(zyp.zhang(x))
          trd[i,j]<-as.numeric(MK["trend"])*monx #per month*12=per year
          sig[i,j]<-as.numeric(MK["sig"])
          if(sig[i,j]<=0.05){sig[i,j]<-0}else {sig[i,j]<-1}
          }
        }
      }
    }
  out<-list(trd=trd,sig=sig)
  return(out)
  }

#uses above to give a spatial object. can also do for seasons
  fx.seas.spatial.trend<-function(Pc.sel,szn=c(12,13,14),nmons=360,div=1,monx=1){#monx=1 year/seas & 12 for months
    #avg=3 for mean =1 for sum
    if(szn[1]==0){pr.MAM<-Pc.sel} else {
      s1<-seq(szn[1],nmons,12);s1l<-length(s1)
      s2<-seq(szn[2],nmons,12);s2l<-length(s2)
      s3<-seq(szn[3],nmons,12);s3l<-length(s3)
      pr.Mar<-Pc.sel[,,seq(szn[1],nmons,12)[1:min(s1l,s2l,s3l)]]
      pr.Apr<-Pc.sel[,,seq(szn[2],nmons,12)[1:min(s1l,s2l,s3l)]]
      pr.Jun<-Pc.sel[,,seq(szn[3],nmons,12)[1:min(s1l,s2l,s3l)]]
      pr.MAM<-(pr.Mar+pr.Apr+pr.Jun)/div
    }
    pr.MAM.trend<-fx.trend.array(pr.MAM,method="zyp",monx=monx)
    sig<-pr.MAM.trend$sig; 
    sigxy<-which(sig==0,arr.ind=TRUE)
    trd<-pr.MAM.trend$trd*monx
    zL<-c((0-max(abs(range(trd,na.rm=TRUE)))),
          (0+max(abs(range(trd,na.rm=TRUE)))))
    out<-list(trd=trd,sig=sig,zL=zL,sigxy=sigxy,pr=pr.MAM)
    return(out)}

#Julian date

fx.julian<-function(date="1970-1-1")
{
    date<-as.character(as.Date(date))
    (Y<-as.numeric(strsplit(date,"-")[[1]][1]))
    (O<-as.Date(paste(Y,1,1,sep="-")))
    J<-julian(as.Date(date),origin=O)+1
    return(J)}

#MAP COUNTRIES    #ADD COUNTRIES TO A MAP     #RETURN AN EXPANDED GRID LON-LAT FOR GIVEN LATS LONS
fx.map.ctry<-function(ctry=c(213,211,105,178,23),lats=seq(-40,40,1),lons=seq(-18,50,1),clr="transparent",lwd=3,global=TRUE,kon=FALSE,newp=FALSE,world_countries=NA){
        if(is.na(world_countries)){world_countries<-get(load("world_countries_24.rda"))}
        library(rgdal); library(clim.pact)
        #print("This function Plots over the current Device if it exists")
        #print("Close any devices and re-plot if a new device is wanted")
        d1<-d2<-NULL
        if(newp==TRUE) #if a new plot is desired
        {
            d1 <- expand.grid(x=lons, y=lats);d2<-cbind(d1$x,d1$y) #grid lat-lon
            names(d2)<-NULL; d2<-as.data.frame(d2);names(d2)<-c("lon","lat");rm(d1)
            plot(d2,pch="",ylab="lat",xlab="lon");
            #add a grid
            lat.grd<-seq(floor(min(lats)),ceiling(max(lats)),1)
            lon.grd<-seq(floor(min(lons)),ceiling(max(lons)),1)
            for( j in 1: length(lat.grd)){abline(h=lat.grd[j],col="grey")}
            for( j in 1: length(lon.grd)){abline(v=lon.grd[j],col="grey")}
    
        }
       
       #Plot all countries
        if(global==T){for(i in 1:dim(world_countries)[1])
        {plot(world_countries[i, ],add=T,lwd=0.1)}}
    
        
        #plot specified countries
        if(kon==TRUE)
        {
            for(i in 1:length(ctry)) #plot for each wanted country
            {
                if(is.numeric(ctry)){ctry2plt<- world_countries[ctry[i], ]}
                if(!is.numeric(ctry)){ctry2plt<- world_countries[world_countries$names == ctry[i], ]}
                plot(ctry2plt,add=T,col=clr,lwd=lwd,ylab="lat",xlab="lon")
            }
        }
        #addland(col="grey")
        
        return(d2)}
    
fx.regrid<-function(x,y,z,xd,yd,res=0.5,plots=c(1,1),kol){
    x1=seq(xd[1],xd[length(xd)],res)
    y1=seq(yd[1],yd[length(yd)],res)
    #Gridded Bivariate Interpolation for Irregular Data
    library(akima)
    
    k<-expand.grid(x,y)
    k<-cbind(k,NA,NA,NA)
    kidx<-1:dim(k)[1]
    eval(parse(text=paste("k[",kidx,",3]<-which(x==k[",kidx,",1])",sep="")))
    eval(parse(text=paste("k[",kidx,",4]<-which(y==k[",kidx,",2])",sep="")))
    eval(parse(text=paste("k[",kidx,",5]<-z[k[",kidx,",3],k[",kidx,",4]]",sep="")))
    k<-k[,c(1,2,5)];names(k)<-c("x","y","z")
    z1<-interp(k$x,k$y,k$z,x1,y1)
    par(mfrow=plots)
    
    resx<-round(abs(mean(x[2:length(x)]-x[1:(length(x)-1)])),2)
    resy<-round(abs(mean(y[2:length(y)]-y[1:(length(y)-1)])),2)
    # image.plot(x,y,z,xlim=xd,ylim=yd,col=kol(20),
        # main=paste("Original", resx,"X",resy),xlab="lon",ylab="lat")
        # contour(x,y,z,add=TRUE);grid();map.ctry()
    # image.plot(z1$x,z1$y,z1$z,xlim=xd,ylim=yd,col=kol(20),
        # main=paste("Regridded", res, "X",res),xlab="lon",ylab="lat")
        # contour(z1$x,z1$y,z1$z,add=TRUE);grid();map.ctry()
    return(z1)}

#Read an esri ascii raster
fx.read.esri<-function(afile,con=TRUE,grd=FALSE)
{
    
    ashp<-as.matrix(read.table(afile,skip=6))
    #ashp[ashp<0]<-NA
    rnum<-dim(ashp)[1]
    cnum<-dim(ashp)[2]
    res<-as.numeric(scan(afile,what="character",skip=4,nlines=1)[2])
    #lon
    x11<-as.numeric(scan(afile, what="c", skip=2,nlines=1)[2])+res/2
    cx<-seq(x11,x11+(cnum-1)*res,res)
    #lat
    y11<-as.numeric(scan(afile, what="c", skip=3,nlines=1)[2])+res/2
    cy<-seq(y11,y11+(rnum-1)*res,res)
    ashp<-t(ashp)[,dim(t(ashp))[2]:1]
    out<-list(x=cx,y=cy,z=ashp)
    if(grd){out<-regrid(x=cx,y=cy,z=ashp)}
    if(con){out$z[out$z<0]<-NA;out$z[out$z>=0]<-1}
    return(out)
}

#plot study area
#coarse rasters on shapes
fx.plot.one<-function(bas.ras,ug,lks,bas,afr,nilb,nilr,ylim=c(-5,6),xlim=c(27,38)){
        par(mfrow=c(1,2))
        #location Map in Africa
        plot(afr,xlim=c(-18,55),ylim=c(-35,35),xaxs="i")
        box();axis(1);axis(2);axis(3);axis(4);grid()
        map.ctry()
        plot(bas,col="green",add=T,,border=2,lwd=1.5)
        plot(nilb,border="cyan1",add=T)
        plot(bas,col="green",add=T,,border=2,lwd=1.5)
        plot(nilr,col="blue",add=T)
        plot(ug,add=T)
        abline(h=0)
        mtext("lat", side=2, line=3, cex.lab=1,las=3, col="blue")
        mtext("lon", side=1, line=3, cex.lab=1,las=1, col="blue")
        
        #Local Basin
        image(bas.ras,col="aliceblue",ylim=ylim,xlim=xlim) 
        image(lks.ras,col="cornflowerblue",add=T)
        plot(bas,add=T,border=2)
        plot(lks,add=T,border=4) 
        plot(Vnile,col=4,add=T)
         map.ctry()
        yg<-ylim[1]:ylim[2]
        xg<-xlim[1]:xlim[2]
        grid(nx=2*(length(xg)-1),ny=2*(length(yg)-1))
        
        mtext("lat", side=2, line=3, cex.lab=1,las=3, col="blue")
        mtext("lon", side=1, line=3, cex.lab=1,las=1, col="blue")
        text(34,5.5,"SUDAN",cex=1.2);text(32,1,"UGANDA",cex=1.2);text(30,-2,"RWANDA",cex=1.2)
        text(30,-3,"BURUNDI",cex=1.2); text(29,1,"DRC",cex=1.2);text(36,3,"KENYA",cex=1.2);text(33,-4.5,"TANZANIA",cex=1.2)
        
        stxyz=fao.et()$xyz
        rn<-floor(runif(10,1,dim(stxyz)[1]))
        for(st in 1:dim(stxyz)[1])
        {
            points(stxyz$lon[st],stxyz$lat[st],pch=16,col=1)
            #if(length(which(rn==st))>0)
            #{text(stxyz$lon[st],stxyz$lat[st],tolower(stxyz$sname[st]),cex=0.8)}
        }
        leg.txt<-c("Lakes & Rivers","Lakes 0.5deg","Basin","Basin 0.5 deg","Met Stations")
        legend("bottomright",leg.txt,col=c(4,"cornflowerblue",2,"grey",1),pch=c(0,15,0,15,16),cex=1,bg="white")
        

        
    }; #plot.one(bas.ras=bas.ras,ug=ug,lks=lks,bas=bas,afr=afr,nilb=nilb,nilr=nilr)
    
    # improved list of objects
fx.ls.objects <- function (pos = 1, pattern, order.by, decreasing=FALSE, head=FALSE, n=5){
        napply <- function(names, fn) sapply(names, function(x)
                                             fn(get(x, pos = pos)))
        names <- ls(pos = pos, pattern = pattern)
        obj.class <- napply(names, function(x) as.character(class(x))[1])
        obj.mode <- napply(names, mode)
        obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
        obj.size <- napply(names, object.size)
        obj.prettysize <- sapply(obj.size, function(r) prettyNum(r, big.mark = ",") )
        obj.dim <- t(napply(names, function(x)
                            as.numeric(dim(x))[1:2]))
        vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
        obj.dim[vec, 1] <- napply(names, length)[vec]
        out <- data.frame(obj.type, obj.size,obj.prettysize, obj.dim)
        names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
        if (!missing(order.by))
            out <- out[order(out[[order.by]], decreasing=decreasing), ]
            out <- out[c("Type", "PrettySize", "Rows", "Columns")]
            names(out) <- c("Type", "Size", "Rows", "Columns")
        if (head)
            out <- head(out, n)
        out
    }
    # shorthand
fx.lsos <- function(..., n=10) {fx.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)}


#INPUT LON,LAT,FIELD & RESAMPLE TO REGULAR GRID AT RESOLUTION RX (GCM to 2.5)
fx.irreg2reg<-function(x,y,field,rx=2.5){
	# x=x-180 in one of my fuctions
    library(clim.pact)
    library(fields)
    yc<-xc<-field
    for( i in 1:dim(field)[1]){for(j in 1:dim(field)[2]){xc[i,j]<-x[i];yc[i,j]<-y[j]}}# series of coords as matrices
    zc<-c(field); xc<-c(xc);yc<-c(yc)
    toplot<-interp(xc,yc,zc, xo=seq(-180+rx/2,180-rx/2,rx),yo=seq(-90+rx/2,90-rx/2,rx)) #make make regular list with x,y,z[,] at rx=2.5
    
}

#CREATE A RASTER OF THE WORLD
#rgl package contains some bitmaps of the world with extents lon(-180,180) lat(-90.90)
#library(rgl) #path in rgl library
#use imagemagick to "> convert world.png world.txt" on cmd line
# file ouput
# ImageMagick pixel enumeration: 1989,994,255,rgb C * R
# 0,0: (255,255,255)  #FFFFFF  white
# 1988,993: (138,214,221)  #8AD6DD  rgb(138,214,221)
fx.bmp2myw<-function(infile,outfile)# infile="C:/B/DATA/ADATA/GISDATA/r-gis/world05.txt"
#outfile="C:/B/DATA/ADATA/GISDATA/r-gis/world05.Rdatadata"
{
    library(clim.pact)
    library(fields)
    #setwd("C:/B/DATA/ADATA/GISDATA/r-gis")
    wld<-scan(infile,skip=1,what="c",sep="\n") #Errors due to (  0,  0,  0) removed by(000,000,000)
    wcol<-vector(length=length(wld))
    for( i in 1:length(wld)){wcol[i]<-strsplit(wld[i]," ")[[1]][length(strsplit(wld[i]," ")[[1]])]}
    wcol[which(wcol=="white")]<-0#water
    wcol[!wcol==0]<-1 #Land
    wcol<-as.numeric(wcol)
    nCR<-strsplit(wld[length(wld)],":")[[1]][1]
    nR<-as.numeric(strsplit(nCR,",")[[1]][1])+1 #1989
    nC<-as.numeric(strsplit(nCR,",")[[1]][2])+1 #994
    myworld<-matrix(wcol,nR,nC)[,nC:1] #Rotated
    #Add coordinates
    lonD<-360/(dim(myworld)[1])#lon Resolution
    latD<-180/(dim(myworld)[2])# lat Resolution
    lonx<-seq(-180+lonD/2,180-lonD/2,lonD);laty<-seq(-90+latD/2,90-latD/2,latD)
    myw<-list(lonx=lonx,laty=laty,myworld=myworld)
    save(myw, file=outfile)
    # To load load("C:/B/DATA/ADATA/GISDATA/r-gis/world.ras.Rdata"); attach(myw)
    # Usage: library(fields); image(myw$lonx,myw$laty,myw$myworld,add=T,col=c("white","transparent")); box()
}


fx.plot.ocean<-function(col=c("white","transparent"),myw=NA)
{
    if(is.na(myw)){myw<-get(load("C:\\DATA\\GIS\\r-gis/world.ras.Rdata"))}
    image(myw$lonx,myw$laty,myw$myworld,add=T,col=col);
}

fx.grd<-function(x,y,res=0.5,col="grey",abso=FALSE,gc=TRUE) #Draw grids
    {
        if(abso==TRUE){x<-x-x%%res;y<-y-y%%res}
        #if(abso==FALSE & gc==TRUE)
#        {
#            xm<-mean(x[2:length(x)]- x[1:(length(x)-1)])
#            ym<-mean(y[2:length(y)]- y[1:(length(y)-1)])
#            x<-seq((x[1]-(0.5*xm)),(x[length(x)]+(0.5*xm)),xm)
#            y<-seq((y[1]-(0.5*ym)),(y[length(y)]+(0.5*ym)),ym)
#        }
        for(g in 1: length(x)){abline(v=x[g],col=col,lty="dashed")} 
        for(g in 1: length(y)){abline(h=y[g],col=col,lty="dashed")} 
    }


fx.plot.res<-function(x,y,z,xd=x,yd=y,ctr,xlim=c(xd[1],xd[length(xd)]),ylim=c(yd[1],yd[length(yd)]),title,af=FALSE,kol){
    
    
    zlim = range(z, finite = TRUE);nlevels=20
    #try(plot(bas,border=2,lwd=2),silent=T)
    image.plot(x,y,z,col=kol(20),zlim=zlim,xlab="lon",ylab="lat",xlim=xlim,ylim=ylim)
    
    try(plot(lks,add=T,border=4,lwd=2),silent=T) 
    try(plot(Vnile,col=4,add=T),silent=T)
    try(plot(bas,add=T,border=2,lwd=2),silent=T)
    try(plot(afr,add=T,border=1,lwd=2),silent=T)
    #map.ctry()
    yg<-ylim[1]:ylim[2]
    xg<-xlim[1]:xlim[2]
        
    #"mintcream"
    if(ctr==TRUE)
    {
        contour(x,y,z,col="grey",nlevels=40,levels = pretty(zlim,nlevels),add=T,labcex=1,, labels=" ",method="flattest",lwd=2)
        contour(x,y,z,col=1,nlevels=40,levels = pretty(zlim,nlevels),add=T,labcex=1,lty=0,method="flattest")
    }
    #map.ctry()
    title(title)
    axis(1,pretty(x))
    axis(2,pretty(y))

}
#plot.res(xd,yd,Mmy,title=paste(gcm[i], "\nEvaporation [mm]\n", yi))


#Draw grids at specified x & y coordinates
   fx.gridxy<-function(x,y,col="grey",lty="dashed",centreline=FALSE,grd.ax=TRUE,newp=FALSE){
    xc<-x;yc<-y
    if(newp==TRUE){plot(expand.grid(x,y),col="transparent",xaxt="n",yaxt="n",xlab="",ylab="",asp=1)}
    if(newp==TRUE & grd.ax==FALSE){plot(expand.grid(x,y),col="transparent",xlab="",ylab="",asp=1)}
    if(newp==FALSE & grd.ax==FALSE){lines(expand.grid(x,y),col="transparent",xlab="",ylab="",asp=1)}
    if(centreline==FALSE)
    {
        y<-y[2:length(y)]-diff(y)/2
        x<-x[2:length(x)]-diff(x)/2
    }
    for(i in 1:length(x)){abline(v=x,col=col,lty=lty)}
    for(i in 1:length(y)){abline(h=y,col=col,lty=lty)}
    if(grd.ax==TRUE & centreline==FALSE){axis(1,x);axis(2,y)}
    if(grd.ax==TRUE & centreline==TRUE){axis(1,xc);axis(2,yc)}
}


fx.grid.xy<-function(x,y){
    yg<-y[1]:y[length(y)]
    xg<-x[1]:x[length(x)]
    grid(nx=2*(length(xg)-1),ny=2*(length(yg)-1))
}
        
fx.plot.regrid.error<-function(kol=topo.colors(20)){
    a<-max(Mmd[,,1])
    b<-min(Mmd[,,1])
#    kol<-topo.colors(20)
    par(mfrow=c(2,2))
    plot(afr)
    image.plot(x,y,Md[,,1],add=T,xlim=c(xd[1],xd[length(xd)]),ylim=c(yd[1],yd[length(yd)]),zlim=c(b,a),xaxt="n",yaxt="n",col=kol)
    plot(afr,add=T)
    plot(afr)
    image.plot(xd,yd,Mmd[,,1],add=T,zlim=c(b,a),xaxt="n",yaxt="n",col=kol)
    plot(afr,add=T)
    image.plot(x,y,Md[,,1],col=kol,zlim=c(b,a),xlim=c(xd[1],xd[length(xd)]),ylim=c(yd[1],yd[length(yd)]),xaxt="n",yaxt="n")#,asp=1)
    plot(afr,add=T)
    plot(lks,border="grey",add=T)
    for(jj in 1: dim(Md[,,1])[1]){for(jk in 1: dim(Md[,,1])[2]){text(x[jj],y[jk],round(Md[jj,jk,1],2),cex=1.5,col=4)}}
    gridxy(x,y)
    abline(h=0,col=2,lwd=2)#;abline(h=4.5);abline(h=-2);
    box()
    image.plot(xd,yd,Mmd[,,1],col=kol,zlim=c(b,a),xlim=c(xd[1],xd[length(xd)]),ylim=c(yd[1],yd[length(yd)]),xaxt="n",yaxt="n")#,asp=1)
    plot(afr,add=T)
    plot(lks,border="grey",add=T)
    for(jj in 1: dim(Mmd[,,1])[1]){for(jk in 1: dim(Mmd[,,1])[2]){text(xd[jj],yd[jk],round(Mmd[jj,jk,1],2),cex=1.5,col=4)}}
    gridxy(xd,yd)
    abline(h=0,col=2,lwd=2)#abline(h=4.5);abline(h=-2);
    box()
}

fx.evap.plot<-function(x,y,z,xd,yd,zlim=range(z,na.rm=TRUE),cont=FALSE,title="",kol=PuBu,add=FALSE,grd=FALSE){
    if(add==FALSE)
    {
        
        image.plot(x,y,z,col=kol(20),xaxt="n",yaxt="n",main=title,
                    xlim=c(xd[1],xd[length(xd)]),ylim=c(yd[1],yd[length(yd)]),zlim=zlim)
        if(grd==TRUE){gridxy(x,y,grd.ax=T,centreline=T)}
#        plot(afr,add=T,lwd=2)
#        plot(lks,border=4,add=T)
#        plot(bas,border="green",lwd=2,add=T)
		fx.map.ctry()
        
        if(cont==T){contour(x,y,z,add=T,col="grey",nlevels=40,xaxt="n",yaxt="n")}
    }
    if(add==TRUE){if(cont==T){contour(x,y,z,add=T,col="grey",nlevels=40,xaxt="n",yaxt="n")}}
    
#    z.im<-im(t(Mmy),xcol=xd,yrow=yd)
#    z.im<-as.im(z.im,xcol=seq(29,39,0.5),yrow=seq(-5,5,0.5))
#    plot(z.im,col=heat.colors(20))
#    gridxy(xd,yd,grd.ax=T,centreline=T)
#    gridxy(xd,yd,grd.ax=F,centreline=F,col=1)
#    plot(afr,add=T,lwd=2)
#    plot(lks,border=4,add=T)
#    plot(bas,border="green",lwd=2,add=T)
#    contour(z.im,add=T,col=1,nlevels=20,xaxt="n",yaxt="n")

}


fx.pick.resample<-function(LON,LAT,Pc,xd,yd,plot=TRUE){
	#resample to xd,yd,PC[xd,yd,t] by picking corresponding value in bigger array LOn,LAT,Pc[LOn,LAT,t]
	PC<-array(dim=c(length(xd),length(yd),dim(Pc)[3]))
	for(tx in 1:dim(Pc)[3]){
		for(ix in 1: dim(PC)[1]){
			for(jx in 1: dim(PC)[2]){
				(ilon<-xd[ix])
				(ilat<-yd[jx])
				(lon.idx<-which(abs(ilon-LON)==min(abs(ilon-LON))))
				(lat.idx<-which(abs(ilat-LAT)==min(abs(ilat-LAT))))
				PC[ix,jx,tx]<-mean(Pc[lon.idx,lat.idx,tx],na.rm=TRUE)
			}
		}
		if(plot==TRUE){image.plot(xd,yd,PC[,,tx],col=PuBu,asp=1);fx.map.ctry()}
	}
	return(PC)
}



#Gridding a matrix with cords x,y,to a new resolution defined by cords xd,yd 
# without interpolating
fx.resample<-function(x,y,z,xd,yd,graph=TRUE)
{
    library(spatstat)
    # Convert z matrix to pixel-image
    z.im<-im(t(z),xcol=x,yrow=y)# 
        #Convert created pixel image to mtrix
        z.im.mat<-t(as.matrix(z.im))
        rownames(z.im.mat)<-x
        colnames(z.im.mat)<-y
        print(z.im.mat) 
    
    #Resample pixellated z.im  to zd.im using as.im and xd,yd
        zd.im<-as.im(z.im,xy=list(x=xd,y=yd))# the pixel image is rotated 
        ##Convert created pixel image to mtrix to check if it was ok
        zd<-t(as.matrix(zd.im))
        rownames(zd)<-xd
        colnames(zd)<-yd
        print(zd)
    if(graph==TRUE){par(mfrow=c(1,2));fx.evap.plot(x,y,z,xd,yd);fx.evap.plot(xd,yd,zd,xd,ydcont=F)}
       
    return(zd)
}


# resampe by picking a value at a point
fx.resample.strip<-function(x,y,z,xd,yd,graph=TRUE,title="",cont=FALSE,kol=PuBu,azlim=range(z,na.rm=TRUE)){
     xyzd<-cbind(expand.grid(xd,yd),NA); names(xyzd)<-c("xd","yd","zd")
     xbounds<-c(x[1],x[2:length(x)]-0.5*diff(x),x[length(x)])
     ybounds<-c(y[1],y[2:length(y)]-0.5*diff(y),y[length(y)])
     for( ii in 1:dim(xyzd)[1])
     {
        xc<-xyzd[ii,1];yc<-xyzd[ii,2]
        xyzd$zd[ii]<-z[min(length(x),max(which(xbounds<=xc))),min(length(y),max(which(ybounds<=yc)))]
     }
     zd<-matrix(xyzd$zd,length(xd),length(yd))
     if(graph==TRUE)
     {
        image.plot(xd,yd,zd,col=kol,asp=1)#par(mfrow=nplots)
        fx.map.ctry()# alim=azlim; #c(min(z,na.rm=T),max(z,na.rm=T))
        # evap.plot(x=xd,y=yd,z=zd,xd,yd,zlim=alim,cont=F,title=title,kol)
        # evap.plot(x,y,z,xd,yd,zlim=alim,title=title,cont=cont,kol,add=TRUE)
        
     }
     return(zd)

}

#remove Area outside basin if res=0.5 of basin raster bas.ras
    fx.rem<-function (xd,yd,z,bas.ras)    {
        if(diff(xd)[1]==0.5 & (diff(xd)[1]==diff(yd)[1]))
        {
            repp<-which(is.na(bas.ras$z),arr.ind=T)
            repp[,1][repp[,1]>dim(z)[1]]<-NA
            repp[,2][repp[,2]>dim(z)[2]]<-NA
            repp<-repp[!apply(repp,1,function(y) any(is.na(y))),]
            z[repp]<-NA
            return(z)
            #Mmy<-z
        }
    };

##Taylor Diagram
fx.taylor<-function(y.df,norm=T,mytitle="\nTaylor Plot",tofile=0)# ydf=data.frame(cbind(date,ref,models))
    {   
        library(plotrix)
        
        x<-y.df[,1] #reference. the first is a date
        if(tofile==1){png(width = 1000, height = 1000,file=paste(mytitle,".png",sep=""), bg="white")}
        
        #if(!is.data.frame(y.df)){y.df<-as.data.frame(y.df)}
        pchm<-c("?",LETTERS)
            
        nf <- layout(matrix(c(1,2),1,2,byrow=TRUE), c(6,2), c(6,6), TRUE)
        par(oma=c(1,1,2,1));#c(bottom, left, top, right)
        layout.show(nf)
        for( n in 1:ncol(y.df))
        {
            y<-y.df[,n]
            if(n==2){TF<-FALSE} else {TF<-TRUE} #Add or Not
            taylor.diagram(x,y,add=TF,col=1,pch=20,pos.cor=TRUE,xlab="",ylab="",main="",
                            show.gamma=F,ngamma=20,sd.arcs=10,ref.sd=T,
                            grad.corr.lines=seq(0.1,0.9,0.1),
                            pcex=1,normalize=norm,mar=c(1,1,2,1))
            R <- cor(x, y, use = "pairwise")
            sd.f <- sd(y); sd.r <-sd(x)
            if(norm == TRUE){sd.f<-sd.f/sd.r}
            text(sd.f*R, sd.f * sin(acos(R)), labels =pchm[n], cex = 0.8, col = 1,adj = -1,offset=4)
             lpos<-sd(x)
             
        }
       points(1,0,pch=10)
        par(mar=c(0,0,0,0))
        frame()
        legend("center",names(y.df),pch=pchm,col=1)
        mtext(mytitle, side=3, line=1, cex=1, col="blue", outer=TRUE)
        box("outer")
        if(tofile==1){dev.off()}
   }


#Some FAO climwat met stations and their reference ET
    fx.fao.et<-function(stfd="C:/B/PAPERS/Evaporation/Data/My_CLIMWAT_Files")
    {
        penf<-list.files(stfd,pattern="pen") #ETO files
        sname<-vector(length=length(penf)); alt<-lat<-lon<-sname
        et<-array(dim=c(length(penf),12))
        for(st in 1:length(penf))
        {
            filz<-paste(stfd,"/",penf[st],sep="")
            penin<-scan(filz,what="c",nlines=1,sep=",",quiet=T)
            lat[st]<-as.numeric(penin[4]) 
            lon[st]<-as.numeric(penin[6])
            sname[st]<-as.character(penin[2])
            alt[st]<-as.numeric(penin[3])
            et[st,]<-read.table(filz,skip=1)[,7]
        }
        colnames(et)<-month.abb
        obs<-cbind("obs",as.data.frame(matrix(colMeans(et),1,12)))
        colnames(obs)<-c("sname",month.abb)
        xyz<-cbind(sname,as.data.frame(cbind(lon,lat,alt,et)))
        et<-rbind(xyz[,c(1,5:16)],obs)
        xyz<-xyz[,1:4]
        xyz$sname<-as.character(xyz$sname)
        
        
        return(list(xyz=xyz,et=et))
    }
    
#plot fao stations
fx.plot.fao<-function(kol=2)
{
    stxyz=fao.et()$xyz
    rn<-floor(runif(10,1,dim(stxyz)[1]))
    for(st in 1:dim(stxyz)[1])
    {
         points(stxyz$lon[st],stxyz$lat[st],pch=16,col=1)
         text(stxyz$lon[st],stxyz$lat[st],tolower(stxyz$sname[st]),col=kol,cex=0.8)
    }
}


#anomaly
fx.anomaly<-function(x){xm<-mean(x,na.rm=T);xa<-x-xm;return(xa)}
fx.anomaly.ts<-function(x)
{
    if(class(x)=="zoo"){x<-aggregate(x,as.yearmon,FUN)}
}

#Change date format in a df[,1] from 1-dec-2007 to 2007-12-01

fx.s2d<-function(fl)
{    
    dates<-as.character(fl[,1])
    for(id in 1: dim(fl)[1])
    {
        txt<-as.character(dates[id])
        txt<-unlist(strsplit(txt,sep))
        dates[id]<-paste(txt[3],which(month.abb==txt[2]),txt[1],sep="-")
    }
    
    return(as.Date(dates))
}
   
#Convert daily or monthly to monthly or annual ts
fx.dmy.ts<-function(in.ts,daily.out="mo",fn="sum")# daily.out="yr"
{
    library(zoo)
    a<-0
    if(!(class(in.ts)=="zoo")){a<-1;mo<-in.ts;yr<-aggregate(mo,f=12,FUN="sum")}
    else {da<-in.ts;mo<-as.ts(aggregate(da, as.yearmon, FUN = "sum" ))
    yr<-aggregate(mo,f=12,FUN="sum")}
    if(a==0){if(daily.out=="mo"){return(mo)} else {return(yr)}}
    if(a==1){return(yr)}
    
}

fx.dmy.mts<-function(in.mts,daily.out="mo",fn=c("sum","mean","mean"))#mts prec-flow-tmpr
{
   prec<-dmy.ts(in.mts[,1])
   flow<-dmy.ts(in.mts[,2],fn="mean")
   tmpr<-dmy.ts(in.mts[,3],fn="mean")

}

kol2=colorRampPalette(c("azure2","cornflowerblue","blue","blue4"),space="Lab")
#plot in mm the map of flow or whatever in m
fx.img<-function(xos,yos,base.ann,kol,about="",X=1000,cont=3, cont.draw=T,
                xlim=c(min(xos),max(xos)),ylim=c(min(yos),max(yos)))
{
    m2p<-base.ann*X
    m2p[m2p==0]<-NA
    image.plot(xos,yos,m2p,xaxs="i",col=kol(20),xlab="lon",ylab="lat",xlim=xlim,ylim=ylim)
    nlevels=length(pretty(unique(sort(c(m2p)))))
    if(cont.draw==T){
    contour(xos,yos,m2p,add=T,nlevels=20,col=cont)
    contour(xos,yos,m2p,add=T,nlevels=20,lty=0,col=1)}
    map.ctry()
    grid()
    box()
    title(about)
}

#Use Point in Polygon to create a matrix representing a SPpolygonDF
#given xc and yc sereis of xoccrdinates
#given a spatialpolygondataframe
#given names of subs plys
#pnames<-c("benga","kwanza")
fx.poly2catch<-function(xo,yo,ang.shp,plots=FALSE,field="ID")
{
      library(grDevices)
      library(maptools)
      library(fields)
      library(raster)
      #Read basin shape(major basins)
      if(is.character(ang.shp)){bas1<-readShapePoly(ang.shp)} else {bas1<-ang.shp} #either read file or sp object
      #Shape limits
      xylim<-rbind(floor(apply(coordinates(bas1),2,min)),ceiling(apply(coordinates(bas1),2,max)))
      xL<-xylim[,1];yL<-xylim[,2]
      #rasterize
      r <- raster(ncols=length(xo), nrows=length(yo),xmn=min(xo),xmx=max(xo),ymn=min(yo),ymx=max(yo)) #empty raster
      R<-rasterize(bas1,r,field=field,silent=TRUE)
      M<-matrix(getValues(R),length(yo),length(xo),byrow=TRUE)
      M<-t(M[dim(M)[1]:1,]) #catchment areas in matrix format
      nkol<-length(unique(c(M)))
     if(plots==TRUE){image(xo,yo,M,col=sample(1:20,nkol))}#col=1:nkol)

#     rows<-length(xc)
#     cols<-length(yc)
#     xyk<-expand.grid(xc,yc)
#     pix<-xyk[,1]
#     piy<-xyk[,2]
#     
#     nshp<-length(ang.shp@polygons)
#     pip<-array(dim=c(rows,cols,nshp))
#     catch<-array(0,dim=c(rows,cols))
#     labxy<-array(dim=c(nshp,2))
#     image(xc,yc,catch,col="transparent")
#     for( i in 1: nshp)
#     {
#         polxy<-ang.shp@polygons[[i]]@Polygons[[1]]@coords
#         plx<-polxy[,1]
#         ply<-polxy[,2]
#         pINp<-point.in.polygon(pix,piy,plx,ply)
#         pINp[pINp==1]<-i
#         pip[,,i]<-matrix(pINp,rows,cols)
#         catch=catch+pip[,,i]
#         pip[,,i][pip[,,i]==0]<-NA
#         labxy[i,]<-c(mean(plx),mean(ply))
#         image(xc,yc,catch,col=c("transparent","azure2"),add=T)
#         try(text(labxy[i,1],labxy[i,2],labs[i]),silent=F)
#     }
#     plot(ang.shp,add=T,col="transparent")
#     grid()
#     catch[catch==0]<-NA
#     labxy<-cbind(as.data.frame(labxy),labs)
#     alist<-list(labxy=labxy,catch=catch)
      alist<-list(x=xo,y=yo,z=M)
      return(alist)
}

#impute with amelia

#given  a data.frame /matrix xwith t at x[,1]
fx.impute<-function(x,m=6,g=F,p2s=0)
{
	if(g==F){graphics.off()} else {print("graphics Off?")
	choice<-scan(n=1,what="c")
	if(choice=="y"|choice=="Y"){graphics.off()}}
	require(Amelia)
	a.out <- amelia(x, ts = names(x)[1],p2s=p2s)
	summary(a.out)
	#try(plot(a.out),silent=T)
	aL<-length(a.out$imputations)
	out<-array(dim=c(dim(x)[1],dim(x)[2]-1,aL))
	for(i in 1:aL){out[,,i]<- as.matrix(a.out$imputations[[i]][,2:dim(x)[2]])}
	out<-apply(out,c(1,2),mean)
	out<-data.frame(x[,1],out)
	names(out)<-names(x)
	dev.new()
	par(mfrow=c((ceiling((dim(x)[2]-1)/2)),2))
	for(i in 2:dim(out)[2]) {
		plot(x[,1],x[,i],lwd=5,type="l",col="grey")
		lines(out[,1],out[,i],col=1,lty="dashed")
		points(out[,1],out[,i])
		}
	return(out)
}


fx.month.days<-function(leap=0) 
{#eneter 366,365 to return month days or nothing for long term
if(leap==0){return(c(31,28.25,31,30,31,30,31,31,30,31,30,31))}
if(leap==366){return(c(31,29,31,30,31,30,31,31,30,31,30,31))}
if(leap==365){return(c(31,28,31,30,31,30,31,31,30,31,30,31))}
}
#source("C:/B/PAPERS/Evaporation/R/Evaporation/Evaporation/fx.gen.evap.r")


#east afruica fao et
fx.fao.eat<-function(stfd="C:/B/PAPERS/Evaporation/climwat")
    {
        penf<-list.files(stfd,pattern="pen") #ETO files
        sname<-vector(length=length(penf)); alt<-lat<-lon<-sname
        et<-array(dim=c(length(penf),12))
        for(st in 1:length(penf))
        {
            filz<-paste(stfd,"/",penf[st],sep="")
            penin<-scan(filz,what="c",nlines=1,sep=",",quiet=T)
            lat[st]<-as.numeric(penin[4]) 
            lon[st]<-as.numeric(penin[6])
            sname[st]<-as.character(penin[2])
            alt[st]<-as.numeric(penin[3])
            et[st,]<-read.table(filz,skip=1)[,7]
        }
        colnames(et)<-month.abb
        obs<-cbind("obs",as.data.frame(matrix(colMeans(et),1,12)))
        colnames(obs)<-c("sname",month.abb)
        xyz<-cbind(sname,as.data.frame(cbind(lon,lat,alt,et)))
        et<-rbind(xyz[,c(1,5:16)],obs)
        xyz<-xyz[,1:4]
        xyz$sname<-as.character(xyz$sname)
        
        
        return(list(xyz=xyz,et=et))
    }
	
fx.implot<-function(lon,lat,m2p,shp,xlim=c(28.5,42.5),ylim=c(-11.5,6.5),
				upper=2000,nl=20,main,lower=floor(min(m2p,na.rm=T)/100)*100)
		{
		ocean<- readShapeSpatial("C:/B/PAPERS/Evaporation/GIS/water_ocean",proj4string=pj4s)
		require(fields)
		#
		#plot a map of a variable within a given shapefile extent 
		#specifying maximum legend value for image plot
			#plot(shp,xlim=xlim,ylim=ylim)
			#lower<-floor(min(m2p,na.rm=T)/100)*100
			image.plot(lon,lat,m2p,main=main,zlim=c(lower,upper),xlim=xlim,ylim=ylim) 
			abline(h=0);plot(shp,add=T);box();grid();
			plot(ocean,col="white",border="grey",add=T,,xlim=xlim,ylim=ylim)
			contour(lon,lat,m2p,add=T,nlevels = nl,col="grey",labels="")
			contour(lon,lat,m2p,add=T,nlevels = nl,lty=0)
			image.plot(lon,lat,m2p,main=main,
			zlim=c(lower,upper),add=T,legend.only=T)
			#title(paste(gcm[i],scen,yearz[1],"-",max(yearz), unt,sep=" "))
		}

library(RColorBrewer)
#Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd 
brew<-c("Blues","BuGn","BuPu","GnBu","Greens","Greys","Oranges","OrRd","PuBu","PuBuGn","PuRd","Purples","RdPu","Reds","YlGn","YlGnBu","YlOrBr","YlOrRd")
for(i in 1: length(brew)){eval(parse(text=paste(brew[i],"<-brewer.pal(9,brew[i])")))}
# par(mfrow=c(3,6))
# for(i in 1: length(brew))
# {
# mypalette<-brewer.pal(9,brew[i])
# image.plot(matrix(1:100,10,10,byrow=T),col=mypalette,main=brew[i])
# }


fx.faotas<-function(stfd="C:/B/PAPERS/Evaporation/Data/My_CLIMWAT_Files")
    {
        penf<-list.files(stfd,pattern="pen") #ETO files
        sname<-vector(length=length(penf)); alt<-lat<-lon<-sname
        tas<-array(dim=c(length(penf),12))
        for(st in 1:length(penf))
        {
            filz<-paste(stfd,"/",penf[st],sep="")
            penin<-scan(filz,what="c",nlines=1,sep=",",quiet=T)
            lat[st]<-as.numeric(penin[4]) 
            lon[st]<-as.numeric(penin[6])
            sname[st]<-as.character(penin[2])
            alt[st]<-as.numeric(penin[3])
            tas[st,]<-rowMeans(read.table(filz,skip=1)[,1:2])
        }
        colnames(tas)<-month.abb
        obs<-cbind("obs",as.data.frame(matrix(colMeans(tas),1,12)))
        colnames(obs)<-c("sname",month.abb)
        xyz<-cbind(sname,as.data.frame(cbind(lon,lat,alt,tas)))
        tas<-rbind(xyz[,c(1,5:16)],obs)
        xyz<-xyz[,1:4]
        xyz$sname<-as.character(xyz$sname)
        
        
        return(list(xyz=xyz,tas=tas))
    }


fx.iplot<-function(x,y,z,labu=TRUE,kol=topo.colors(20),cex=0.8)
	{
		require(fields)
		#special image plot
		image.plot(x,y,z,zlim=c(min(z,na.rm=T),max(z,na.rm=T)),col=kol)
		if(labu==TRUE){} else{contour(x,y,z,zlim=c(min(z,na.rm=T),max(z,na.rm=T)),add=T)}
		grid(nx=dim(z)[1],ny=dim(z)[2])
		if(labu==TRUE){text(expand.grid(x,y)[,1],expand.grid(x,y)[,2],round(c(z),1),cex=cex)}
	}

fx.Rmean<-function(y){matrix(rep(apply(y,1,mean),dim(y)[2]),dim(y)[1],dim(y)[2])}
fx.Cmean<-function(y){matrix(apply(y,2,mean),dim(y)[1],dim(y)[2],byrow=T)}
fx.Mmean<-function(y){matrix(mean(y,na.rm=T),dim(y)[1],dim(y)[2],byrow=T)}

#xy matrix to yx ascii- esri grid
fx.xy2ascGRDyx<-function(xllcorner,yllcorner,cellsize,xy){
	file=paste(deparse(substitute(xy)),"asc",sep=".")
	hdr1<-c("ncols","nrows","xllcorner","yllcorner","cellsize","NODATA_value")
	yx<-t(xy[,dim(xy)[2]:1])#Rotate anti clockwise
	ncols<-dim(yx)[1]
	nrows=dim(yx)[2]
	hdr2<-c(ncols,nrows,xllcorner,yllcorner,cellsize,-9999)
	hdr<-cbind(hdr1,hdr2)
	write.table(hdr,file=file,quote=F,row.names=F,col.names=F)
	write.table(yx,file=file,quote=F,row.names=F,col.names=F,append=T)
}

#xy matrix to yx ascii- esri grid
fx.yx2ascGRDyx<-function(xllcorner,yllcorner,cellsize,yx){
	file=paste(deparse(substitute(xy)),"asc",sep=".")
	hdr1<-c("ncols","nrows","xllcorner","yllcorner","cellsize","NODATA_value")
	ncols<-dim(yx)[1]
	nrows=dim(yx)[2]
	hdr2<-c(ncols,nrows,xllcorner,yllcorner,cellsize,-9999)
	hdr<-cbind(hdr1,hdr2)
	write.table(hdr,file=file,quote=F,row.names=F,col.names=F)
	write.table(yx,file=file,quote=F,row.names=F,col.names=F,append=T)
}


#yx ascii- grid to xy matrix
fx.asc2xy<-function(filepath){
	require(maptools)
	dem<-readAsciiGrid(filepath)
	a<- unlist(strsplit(filepath,"/"))
	a<-a[length(a)]
	txt<-paste("dem<-matrix(dem@data$",a,", dem@grid@cells.dim)",sep="")#by col
	eval(parse(text=txt))
	dem<-dem[,dim(dem)[2]:1]
	return(dem)
}


#plot a matrix and its values in the way it looks
fx.plot.mat<-function(z,x=NULL,y=NULL,ret=T,main=deparse(substitute(z)),col=topo.colors(20),labs=F,ctr=F,ctr.col="white")
{
	require(fields)
	Z<-t(z[dim(z)[1]:1,])
	if(is.null(x)){x<-1:dim(Z)[1]}
	if(is.null(y)){y<-1:dim(Z)[2]}
	image.plot(x,y,Z,main=main,col=col)
	if(labs==T){for(i in 1:dim(Z)[1]){for(j in 1:dim(Z)[2]){text(x[i],y[j],Z[i,j])}}}
	mat<-list(x=x,y=y,z=Z)
	if(ctr==T){contour(x,y,Z,col=ctr.col,add=T)}
	if(ret==T){return(mat)}	
}

#plot with axes vertical
fx.plotWvxdate<-function(M,n=31,add=F){
	if(add==F){plot(M,xaxt="n",xlab="")} else {lines(M,xaxt="n",xlab="")}
	axis(1,time(M)[seq(1,length(M),n)],paste(year(time(M))[seq(1,length(M),n)],month.abb[month(time(M))[seq(1,length(M),n)]],sep="-"),las=2)}

#markov simulation
#P <- matrix( c(0.2, 0.5, 0.5, 0.1, 0.8, 0.1, 0.2, 0.1, 0.7), 3,3, byrow=TRUE)
 fx.simMarkov <- function( P,states, len=1000) {
        if(length(states)!=NROW(P)){states<-1:NROW(P)}
        result <- numeric(len)
        result[1] <- 1
        for (i in 2:len) {
            result[i] <- sample(states, 1, prob=P[ result[i-1], ])
        }
        result
 }
 
# Writes a data frame to an ASCII raster file, suitable for display in ArcView or ArcGIS.
fx.matrix2ascRas<-function(x,y,xyRC,file,cellsize){
	require(epiR)
	xllcorner<-x[1]-(0.5*cellsize)
	yllcorner<-y[1]-(0.5*cellsize)
	epi.asc(xyRC, file, xllcorner, yllcorner, cellsize, na = -9999)
	print(file)
}

#Compute the anomaly with 1961-1990 as reference period
fx.anom.daily<-function(x,FUN=mean){#x is a zoo object
	library(hydroTSM)
	library(tis)
	m<-month(time(x))
	base<-window(x,start=as.Date("1961-01-01"),end=as.Date("1990-12-31"))
	base.ltm<-monthlyfunction(base,FUN)
	base.x<-base.ltm[m]
	anom<-x-base.x
	return(anom)
}

fx.anom<-function(x,fun=sum){#x is a zoo object
	library(hydroTSM)
	library(tis)
	m<-month(time(x))
	base<-base.ltm<-0
	if(is.zoo(x)){
		base<-window(x,start=as.Date("1961-01-01"),end=as.Date("1990-12-31"))
		y<-daily2monthly(x,FUN=fun)
		base.ltm<-monthlyfunction(y,mean)
	}
	
	if(is.ts(x)) {
		base<-window(x,start=c(1961,1),end=c(1961,1))
		base.ltm<-monthlyfunction(base,mean)
	}
	
	base.x<-base.ltm[m]
	anom<-x-base.x
	return(anom)
}

#Confidence INtervak plot
fx.ciplot<-function(x,y,model="gam",ci=1.96,xlab="",ylab="",Title="",lposi=0,add=FALSE,pcol="grey",pch="o",lcol=1,xL= c(min(x),max(x)),border=pcol){
  if(model=="gam"){library(mgcv); model <- gam( y ~ s(x))}
  # get the estimated values of the fitted smoothing spline
  fit <- predict( model , se = TRUE )$fit
  # and pointwise standard errors
  se <- predict( model , se = TRUE)$se.fit
  # calculate the values of the upper and lower 95%-confidence limits
  lcl <- fit - ci * se
  ucl <- fit + ci * se
  # # plot an empty coordinate system with right scaling, labels , titles and so on. I prefer to include axes = FALSE and add the axes manually. See here for an example.
  if(add==TRUE){
    par(new=TRUE);
    plot( 0 , type = "n" , bty = "n" , xlab = "" , ylab = "" , main ="" , xlim = xL, ylim = c(min(y),max(y)),yaxt="n" )
    # Add axis for added plot
    axis(4)
    mtext(ylab,4,line=lposi,outer=TRUE)
  
  } else {
    par(oma=c(1,1,1,1))
    plot( 0 , type = "n" , bty = "n" , xlab = xlab , ylab = ylab , main = "" , xlim = xL, ylim = c(min(y),max(y)) )
  }
  # o it gets a bit tricky with sorting the coordinates in the right order. R provides a very effective way to achieve this, though it is not easily understandable at once. We create two indices on the dataset i.for and i.back which number the dataset according to the independend variable on the x-axis - forward and backward respectively:
  i.for <- order(x)
  i.back <- order(x , decreasing = TRUE )
  # The points of the shade polygon follow here first the upper confidence line (ucl) forward and then the lower confidence line backward. The coordinates of the shade-polygon are
  x.polygon <- c( x[i.for] , x[i.back] )
  y.polygon <- c( ucl[i.for] , lcl[i.back] )
  # # First plot the polygon with the confidence shade - otherwise the GAM plot vanishes behind it.
  if(add==FALSE){
    polygon( x.polygon , y.polygon , col = pcol , border = border )
    # now the mainplot afterwards
    lines( x[i.for] , fit[i.for], col = lcol , lwd = 3 )
    
    } else {
    polygon( x.polygon , y.polygon , col = pcol , border = border ,yaxt="n",ylab="n",xlim=xL)
    # now the mainplot afterwards
    lines( x[i.for] , fit[i.for], col = lcol , lwd = 3 ,yaxt="n",ylab="n",xlim=xL)
    title(Title,outer=TRUE)
    }  
    
  #plot original
  points(x,y,pch=pch,col=lcol,cex=0.5)
  return(model)
}

#Filling by regression by month to get trend line and sapling residuals(inout Time series)
fx.filltsReg<-function(RAIN){
  require(tis)
  X<-time(RAIN)
  for(m in 1:12){
    im<-which(month(X)==m)
    y<-RAIN[im]; im2<-which(!is.na(y)); y2<-y[im2]
    x<-X[im];x2<-x[im2]
    (reg<-lm(y2~x2))
    (res<-residuals(reg))
  }
  
}

#make a zoo object start from 1961-2010 and save a textfile
fx.zoo6190<-function(in.zoo,tdir,sname,vnam,LON="lon",LAT="lat",ALT="alt",ends=as.Date("2010-12-31")){
  require(tis)
  base<-zoo(NA,seq(as.Date("1961-01-01"),ends,1))
  out.zoo<-cbind(in.zoo,base)[,1]  
  #in.zoo lat
	attributes(out.zoo)$lon<-lon<-LON#lon
	attributes(out.zoo)$lat<-lat<-LAT#lat
	attributes(out.zoo)$alt<-alt<-ALT
	#SAVE
  (outfile=paste(tdir,sname,"_",lon,"_",lat,"_",alt,"_",1961,"_",year(ends),"_",vnam,".dat",sep=""))
  Arua.out<-out.zoo; Arua.out[is.na(Arua.out)]<-(-999)
	write.table(Arua.out, file=outfile,col.names=FALSE,row.names=FALSE)
  print("**************")
  print(unlist(strsplit(outfile,"/"))[length(unlist(strsplit(outfile,"/")))] )
  print("**************")
  #"print(list.files(tdir,pattern=".dat"))
  plot(out.zoo)
  return(out.zoo)
  }
# 
fx.pointfromarray<-function(x,y,X,Y,Z){
	  #Extract the series of value at the given coordinates of an array/matrix
	  #given x,y and the matrix coordinates X,Y and values Z
	  nearest<-function(x,X){return(ix<-which(abs(X-x)==min(abs(X-x)))[1])}
	  xi<-nearest(x,X)
	  yi<-nearest(y,Y)
	  z<-0
	  if(length(dim(Z))==3){z<-Z[xi,yi,]} else {z<-Z[xi,yi]}
	  return(z)	
  }
#
fx.extract.season.array<-function(Pc.sel,DJF,fun=sum,zlim=c(0.1,1000),lmar=1,bordcol="black",indeje){
	#Given an array with jan--, and an index in 2 col matrix form with each year index in a row e.g for DJF 12|14, 24:26 etc
	#exract a sum for each season each year
	library(abind)
	#indeje<-readShapePoly("C:/DATA/Workspace/ICH_SA/GIS/Indeje_shp.shp")
	for(i in 1:dim(DJF)[1]){
		if(i==1){Pc.sel.DJF<-apply(Pc.sel[,,DJF[i,1]:DJF[i,2]],c(1,2),fun)}
		if(i>1){
			temp<-apply(Pc.sel[,,DJF[i,1]:DJF[i,2]],c(1,2),fun)
			Pc.sel.DJF<-abind(Pc.sel.DJF,temp,along=3)}
	}
	MAIN<-deparse(substitute(DJF))
	PcDJFann<-apply(Pc.sel.DJF,c(1,2),mean)
	if(lmar==1){image.plot(lon,lat,PcDJFann,col=c("white",PuBu),nlevel=10,axis.args = list(cex.axis = 0.8),cex=0.5,axes=FALSE,xlab="",ylab="",main=MAIN,zlim=zlim,asp=1);}
	if(lmar==0){image(lon,lat,PcDJFann,col=c("white",PuBu),cex=0.5,axes=FALSE,xlab="",ylab="",main=MAIN,zlim=zlim,asp=1);}
	fx.map.ctry();sp::plot(indeje,add=TRUE,border=bordcol,lty=2,lwd=1);box()
	return(Pc.sel.DJF)
}


fx.plot.mts<-function(zoots,lpos="top",kol=1:ncol(zoots),Lwd=1:ncol(zoots),Lty=1:ncol(zoots),unt="DegC",DL=31){
      #Plot a multiple time series
      #zoots is the multiple time series
      #kol=colours
      #Lwd=line widths
      #Lty=line types
      #unt=y-axis label
      #DL=spacing of x-axis lable in days
      #Usage
      # fx.plot.mts(z.day,Lwd=c(2,2,2,3,3,1),unt="mm",kol=c(1,2,3,4,"Orange"))
      # fx.plot.mts(z.mo,DL=1)
      #lpos="topright"#legend position

      plot(zoots[,1],col=kol[1],lty=Lty[1],lwd=Lwd[1],ylim=c(min(zoots,na.rm=TRUE),max(zoots,na.rm=TRUE)),xaxt="n",ylab=unt,xlab="")

      for(i in 2:dim(zoots)[2]){lines(zoots[,i],col=kol[i],lty=Lty[i],lwd=Lwd[i])}

      (xL<-time(zoots)[seq(1,length(time(zoots)),DL)])

      library(tis)

      (txt<-paste(month.abb[month(xL)],year(xL),sep="-"))

      axis(1,xL, txt,las=2)

      abline(v=xL,lty=2,col="grey")

      abline(h=pretty(zoots),col="grey",lty=2)

      legend(lpos,names(zoots),col=kol,lty=Lty,lwd=Lwd,bg="white")

}

 cru.retrieve<-function(fil, #matrix of Data rep(nLat X nLon)*nMonths
                      lonx,laty,xL,yL,MULTI,
                      y1=1901,y2=2006,
                      xCRU=seq(-179.75,180,0.5),#longitudes of data
                      yCRU=seq(-89.75,90,0.5),#latitudes of data
                      nais=-999,
                      nMonths=1272,
                      col=rgb( colorRamp(c("azure", "blue"))(seq(0, 1, length = 100)), max = 255),
                        plt=FALSE)
   {
  library(clim.pact)
  f<-(get(load(fil))) #34 secs
  f<-matrix(f,457920,720,byrow=TRUE)
  f[f==nais]<-NA;f<-f*MULTI
  xy<-expand.grid(lonx,laty)
  nLat<-length(yCRU)
  nLon<-length(xCRU)
  nROWS<-nLat*nMonths
  r1<-seq(1, 457920,nLat)#sequence to define addresses of first rows
  r2<-seq(nLat, 457920,nLat)#sequence to define addresses of last rows
  labs<-paste(sort(rep(y1:y2,12)),month.abb[rep(1:12,length(y1:y2))])
  if(plt==TRUE){
    image.plot(xCRU,yCRU,t(f[r1[1]:r2[1],]),col=col,main=labs[1],
        zlim=c(0.01,max(f[r1[1]:r2[1],],na.rm=TRUE)),xlim=xL,ylim=yL)
  addland()
  }
  VAL<-array(dim=c(length(lonx),length(laty),nMonths))
  
  for( k in 1:dim(xy)[1]){ #for each point
    val<-vector(length=nMonths)#months in dataset
    lon<-xy[k,1]
    lat<-xy[k,2]
    if(plt==TRUE){points(lon,lat,pch=".")}
    (Ri<-myk(lat,yCRU)) #Row number of point in first month of f[457920,720]
    (Ci<-myk(lon,xCRU)) # Col number of point
    (Ris<-seq(Ri,457920,360))
    idx<-1:1272
    eval(parse(text=paste("val[",idx,"]<-f[Ris[",idx,"],Ci]",sep="")))
    val<-ts(val,start=c(1901,1),frequency=12)
    VAL[which(lonx==lon),which(laty==lat),]<-val
  }
  return(VAL)
}


#setwd(d1)
#
##make Package
# package.skeleton("fxgenr")
#setwd(d1)

#hypsographic curve given z
fx.hypso<-function(z){
  z<-sort(z);L<-length(z);diffz<-length(z)/10; 
  s2<-floor(seq(diffz,L,diffz));H<-data.frame(seq(10,100,10),z[s2])
  names(H)<-c("%","masl")
  plot(H,type="l",main=paste(cname,"\nHypsographic curve"),xlab="%",ylab="masl");grid();
  return(H)
}

#WBM Functions
Favg<-function(Qobs,Qsim,rtype=1){
  #Objective function, which is the average (Favg) of four measures modified from Zhang et al. (2008)(13)
  Qobs[Qobs<=0]<-NA;Qsim[Qsim<=0]<-NA;
  library(topmodel)
  F1<-NSeff(Qobs,Qsim) #Nash Sutcliffe (NS) gives highest weight to large errors (often but not always happen during periods of high flow). 
  F2<-0
  
  #F2 is the NS of logarithmically transformed flows with greater weight given to errors during low flows. F3 is the Pearson correlation coefficient, which measures the covariability of the simulated and observed without penalizing for bias. F4 is a symmetric measure of the match between the average simulation and average observation (i.e. a Bias Score). All of the measures are "skill scores" (Wilks, 2006) in that they measure the relative accuracy of the simulations, with respect to a baseline (in this case, the long term average). Throughout this article, the terms "performance" and "skill" are used interchangeabl
  
  
}
sArray<-function(x,fname){
  Tim<-dim(x)[3]
  write.table(x[,,1],file=fname,quote=FALSE)
  for(i in 2:Tim){
    write.table(matrix("--",1,dim(x)[2]),file=fname,quote=FALSE,append=TRUE,col.names=FALSE)
    write.table(x[,,i],file=fname,quote=FALSE,append=TRUE,col.names=FALSE)
  }
}
#CRU RetRIEVE
myk<-function(myc,x)  {#given lat or lon outputs index in matrix
  
  bx<-c(x-diff(x)[1]/2, max(x)+diff(x)[1]/2)
  res<-which(x==mean(bx[c(max(which(bx<myc)),min(which(bx>=myc)))]))
  return(res)
}
#Min & max functions for apply
maxfun<-function(x){max(x,na.rm=TRUE)}
minfun<-function(x){min(x,na.rm=TRUE)}
options(width = 500)  
fx.Monthdays.Daylight<-function(latd,TQ){ #givena latotude and a month-date
  j=pi/180*latd
  #DATE
  (mDAY<-as.Date(as.yearmon(TQ))+14)
  (Yr<-tis::year(mDAY))
  (Mon<-tis::month(mDAY))
  Ym<-cbind(Yr,Mon)
  #Days in Month
  d1<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  d2<-c(31,29,31,30,31,30,31,31,30,31,30,31)
  mdays<-function(x) if(x[1]%%4==0){return(d2[x[2]])}else{return(d1[x[2]])}
  (d<-apply(Ym,1,mdays))
  if(is.matrix(TQ)){(d<-matrix(d,dim(TQ)))}
  #(d<-as.POSIXlt(seq(as.Date(TQ), by="month",length.out=12)-1)$mday[tis::month(as.Date(TQ))])
  #Dayligh hrsin units of 12hrs
  (J1<-tojul(tis::year(mDAY),1,1)) #january 1st Julian Day)
  (J<-tojul(tis::year(mDAY),tis::month(mDAY),tis::day(mDAY))-J1+1)# Julian Day
  #inverse relative distance EadSh-Sun, dr,
  (dr<-1+0.033*cos(2*pi/365*J))
  #solar declination, d
  (decl<-0.409*sin(2*pi/365*J-1.39))
  #sunset hour angle, w s
  (X<-max(0.00001,1-(tan(j))^2 *(tan(decl))^2))
  (ws<-pi/2-atan(-tan(j)*tan(decl)/(X^0.5)))
  #Daylight Hours
  (D<-24/pi*ws)
  #maximum possible duration of sunshine
  out<-list(d=d,D=D)
}
fx.mpz<-function(dem,eaf,basin=NULL,xlim=NA,ylim=NA,main=cname,colf=rainbow){#plot a terrain dataset
  if(is.na(xlim[1])){xlim<-as.numeric(bbox(dem)[1,])}
  if(is.na(ylim[1])){ylim<-as.numeric(bbox(dem)[2,])}
  slope <- terrain(dem, opt='slope')
  aspect <- terrain(dem, opt='aspect')
  hill <- hillShade(slope, aspect, 40, 270)
  plot(hill, col=grey(0:100/100), legend=FALSE, main=main,xlim=xlim,ylim=ylim)
  plot(dem, col=colf(32, alpha=0.35), add=TRUE)#  plot(dem,col=terrain.colors(32));
  plot(eaf,border="transparent",add=TRUE);
  plot(wat,add=TRUE,col="blue",border="transparent",,xlim=xlim,ylim=ylim);#col=rgb(0,0,1,alpha=0.2)
  plot(ctry,add=TRUE);#plot(cat,col="grey",add=TRUE)
  if(!is.null(basin)){plot(basin,border="cyan",lwd=2,col="transparent",add=TRUE)}
  grid()
}


fx.areasqm=function(basin){#SPPDF object ##Area of basin
  basin.UTM=spTransform(basin,CRS("+proj=utm +zone=36N +datum=WGS84"))
  area<-sapply(basin.UTM@polygons, function(x) x@area)
  return(area)
}

#read flow files and create ts -objects
fx.read.flo=function(flofile="Akagera-Rusumo-SL70003.xls",sname,kagera=TRUE){
	if(kagera==TRUE){
		#RUNOFF (RUSUMO)
		Q<-read.xls(paste("C:/DATA/HYDROMET/FLOW/Kagera/",flofile,sep=""),sheet=1); 
		Q[,1]<-as.character(Q[,1])
		names(Q)[1:4]<-c("Date","Min","Mean","Max")
		head(Q)
		
		#Qdate format= mm/yyyy
		if(substr(Q$Date[1],3,3)=="/"){
			date1<-as.numeric(unlist(strsplit(Q$Date[1],"/")))#Kagera
			month1<-date1[1]; year1<-date1[2]
		}
		
		#date format= yyyy-mm-dd
		if(substr(Q$Date[1],5,5)=="-"){
			(date1<-as.Date(Q$Date)[1])
			library(tis)
			(month1<-tis::month(date1))
			(year1<-tis::year(date1))
		}
		
		Qobs<-ts(Q[,3],start=c(year1,month1),f=12)
		} else {
			#RUNOFF (EAPMP)
			Qin<-read.table(paste("C:\\DATA\\HYDROMET\\FLOW\\EAPMP\\",flofile,sep=""),sep="|",header=TRUE)
			T1<-c(Qin[1,1],1)
			Qobs<-ts(data=c(t(as.matrix(Qin[,-1]))),start=T1,frequency=12)
		}
	
	out=list(sname=sname,flo=Qobs)
	save(out,file=paste(sname,".flo.rdata",sep=""))
}

fx.array2mts<-function(Array,rem.arr){#turmns array in series by column, eliinates NA cells in rem.arr
  mts<-t(apply(Array,3,c)); mts<-mts[1:length(mts)]
  rem<-t(apply(rem.arr,3,c)); rem<-rem[1:length(rem)]
  valid<-which(!is.na(rem))
  mts<-mts[valid]
}

#find the closets value to a number within a sereis
fx.closest<-function(ext,lon){a<-abs(lon-ext);i<-which(a==min(a));return(i)} #find closest value in series

#gibe a shapefile, lom,lat sereis bouding and time/3rd dim turn shapefile to raster
fx.polygon2raster2array<-function(shp=ind1,xo,yo,cTim){#shape,lon,lat,time-length or array-3rd dimension
	ras <- raster(ncols=length(xo), nrows=length(yo))
	limRas<-c(min(xo)-diff(xo)[1]/2,max(xo)+diff(xo)[1]/2,min(yo)-diff(yo)[1]/2,max(yo)+diff(yo)[1]/2);#raster limits
	(limRas<-matrix(limRas,2,2,byrow=TRUE))
	(colnames(limRas)<-c("min","max"))
	(rownames(limRas)<-c("x","y"))
	(shp@bbox<-limRas)
	(extent(ras)<-extent(limRas))#raster extent
	ras[]<-1#t(R[,ncol(R):1]);
	bas<-rasterize(shp,ras)
	bas<-as.matrix(bas)#basin
	bas<-t(bas[nrow(bas):1,])
	#Matrix to remove outer cells
	bas[!is.na(bas)]<-1;
	bas<-array(bas,c(dim(bas),cTim))
	return(bas)
}


#retrieve.nc to a list (lon,lt,Tim,val[lon,lat,Tim])
#Retrievs to model resolution
fx.retrieve.nc=function(file, #nc-file
	LONX, #longitude required
	LATX,#latitude required
	shp=NULL,#shapefile to plot
	iplot=1, #Timesteps to plot
	vname=NA, scen="", gcm="",#variable name (optional), scenario,model
	multip=c(30*24*3600,12,1), #multiply to get units (kg/m2/s to mm/month), temp=1, hum=?,wind=?,12 prcp 3 seas prc else all=1
	kol1=rgb( colorRamp(c("azure", "blue"))(seq(0, 1, length = 100)), max = 255),
	kol2=rgb( colorRamp(c("yellow","red")) (seq(0, 1, length = 100)), max = 255),
	resample=TRUE){ #Resampling to LON,LAT by interp
	#Extents
	(lonEx<-c(min(LONX),max(LONX)))
	(latEx<-c(min(LATY),max(LATY)))
	#libraries
	library(ncdf)
	library(zoo)
	library(fields)
	library(raster)
	#open file
	d<-open.ncdf(file);
	if(is.na(vname)){vname=names(d$var)[length(names(d$var))]}
	#Dimensions of GCM
	LON<-d$dim$lon;
	lon<-LON$vals;nx=length(lon)
	LAT<-d$dim$lat;
	lat<-LAT$vals;ny=length(lat)
	TIM<-d$dim$time;
	tim<-TIM$vals;
	torg<-as.Date(unlist(strsplit(TIM$units,"since "))[2])
	tim<-as.yearmon(torg+tim)
	#variable
	#time limits
	start.end=as.yearmon(c(as.Date("1961-1-1"),as.Date("1990-12-31")))
	IDt<-match(start.end,tim); 
	Tim<-tim[IDt[1]:IDt[2]] #subset Time
	cTim=length(Tim) #count Time
	#Lon /lat limits
	fx.closest<-function(ext,lon){a<-abs(lon-ext);i<-which(a==min(a));return(i)} #find closest value in series
	(IDlon<-c(fx.closest(lonEx[1],lon),fx.closest(lonEx[2],lon)))
	(IDlat<-c(fx.closest(latEx[1],lat),fx.closest(latEx[2],lat)))
	lonx<-lon[(IDlon[1]:IDlon[2])]; clonx=length(lonx) #subset latitude
	laty<-lat[(IDlat[1]:IDlat[2])];claty=length(laty) #subset latitude
	#get value ##Val[lon,lat,time]
	start=c(IDlon[1],IDlat[1],IDt[1])
	count=c(clonx,claty,cTim)
	val<-get.var.ncdf(d,start=start,count=count)*multip[1]; #Val[lon,lat,time]
	val[val<0]<-0
	out=list(lon=lonx,lat=laty,tim=Tim,val=val) #At GCM resolution
	
	#Resample to a given LONX,LATY,& cut to the map
	if(resample==TRUE){
		s=array(0,dim=c(length(LONX),length(LATY)))
		s<-raster(s)# new resolution
		val<-brick(val)
		val<-raster::resample(val,s,method="bilinear")
		val<-as.array(val)
		
	
		
		if(!is.null(shp)){
			#Shape to raster to Matrix##########################################
			#Shape to raster at resolution diff(x) within xmn-xmx & ymn-ymx
			# Create a empty raster based on coordinates
			ras <- raster(ncols=length(LONX), nrows=length(LATY))
			limRas<-c(min(LONX)-diff(LONX)[1]/2,max(LONX)+diff(LONX)[1]/2,
				min(LATY)-diff(LATY)[1]/2,max(LATY)+diff(LATY)[1]/2);#raster limits
			(limRas<-matrix(limRas,2,2,byrow=TRUE))
			(colnames(limRas)<-c("min","max"))
			(rownames(limRas)<-c("x","y"))
			(shp@bbox<-limRas)
			(extent(ras)<-extent(limRas))#raster extent
			ras[]<-1#t(R[,ncol(R):1]);
			bas<-rasterize(shp,ras)
			bas<-as.matrix(bas)#basin
			bas<-t(bas[nrow(bas):1,])
			#Matrix to remove outer cells
			bas[!is.na(bas)]<-1;
			bas<-array(bas,c(dim(bas),cTim))
			val<-bas*val#remove areas outside
		}
		
		out=list(lon=LONX,lat=LATY,tim=Tim,val=val,cut=bas) #Resampled
	}
	
	#Seasonal values
	fx<-function(x) mean(x,na.rm=TRUE)
	out$ann<-multip[2]*apply(out$val,c(1,2),fx)
	out$DJF<-multip[3]*apply(val[,,c(seq(12,cTim,12),seq(1,cTim,12),seq(2,cTim,12))],c(1,2),fx)
	out$MAM<-multip[3]*apply(val[,,c(seq(3,cTim,12),seq(4,cTim,12),seq(5,cTim,12))],c(1,2),fx)
	out$JJA<-multip[3]*apply(val[,,c(seq(6,cTim,12),seq(7,cTim,12),seq(8,cTim,12))],c(1,2),fx)
	out$SON<-multip[3]*apply(val[,,c(seq(9,cTim,12),seq(10,cTim,12),seq(11,cTim,12))],c(1,2),fx)
	szlim=range(c(out$DJF,out$MAM,out$JJA,out$SON),na.rm=TRUE)
	
	if(vname=="pr"){kol=kol1}else {kol=kol2}
	fx.seas.plot<-function(X,tt,kol,szlim){
		image.plot(out$lon,out$lat,X,asp=1,xlab="lon",ylab="lat",col=kol,zlim=szlim,horizontal=TRUE);
			if(!is.null(shp)){plot(shp,add=TRUE)}
			title(toupper(paste(tt)),cex=0.5)
	}
	
	##Plot if iplot is given
	if(!is.na(iplot)){
		for(k in iplot){
			zl<-range(out$val,na.rm=TRUE)
			image.plot(out$lon,out$lat,out$val[,,k],asp=1,xlab="lon",ylab="lat",col=kol,zlim=zl,horizontal=TRUE);
			if(!is.null(shp)){plot(shp,add=TRUE)}
			title(toupper(paste(vname,Tim[k],scen,gcm,sep="|")))
		}
		#se3asonal plots
		fx.seas.plot(out$ann,paste(vname,gcm,"\nANNUAL 1961-1990"),kol,szlim=range(out$ann,na.rm=TRUE))
		par(mfrow=c(2,2))
		fx.seas.plot(out$DJF,"DJF 1961-1990",kol,szlim)
		fx.seas.plot(out$MAM,"MAM 1961-1990",kol,szlim)
		fx.seas.plot(out$JJA,"JJA 1961-1990",kol,szlim)
		fx.seas.plot(out$SON,"SON 1961-1990",kol,szlim)
		title(toupper(paste("\n",vname,scen,gcm,sep="|")),outer=TRUE)
		par(mfrow=c(1,1))
	}
	
	return(out)
}


