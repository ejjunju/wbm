x11()
library(zoo)
#WBM FUNCTIONS
source('C:/DATA/Rscript/04_hydModel/wbm/wbmfxEPOT.R')
source('C:/DATA/Rscript/04_hydModel/wbm/SNOWROUTINE.r')
source('C:/DATA/Rscript/04_hydModel/wbm/wbmThornthwaite.r')


#load('testData.mat') Load data if it exists
setwd('C:\\DATA\\STATIONS\\wbm_Rusumo\\Calib\\matlab')

#DATA######################################################
###########################################################
head(init<-read.table("init.txt",header=TRUE))
eval(parse(text=paste(names(init),"<-init$",names(init))))
rm(init)

#FLOW####################################################
head(obs<-read.table("Qobs.txt",header=TRUE))
eval(parse(text=paste(names(obs),"<-obs$",names(obs))))
Qobs<-as.ts(zoo(Qobs,as.yearmon(Time)))
rm(obs)
#subset Qobs
#(Can subset nm to nms)
qobs=window(Qobs,start=c(1961,1),end=c(1970,1))
(startsel=as.Date(as.yearmon(time(qobs1)[1])));
(Timestr=as.Date(as.yearmon(time(Qobs))));
(id=match(time(qobs),time(Qobs)))
(newtime=time(Qobs)[id]);
(nms=length(qobs))

#MET##############################################
head(MET<-read.table("met.txt",header=TRUE))
met=vector("list",length=nrow(MET)/length(Qobs))
mt<-array(dim=c(nm,5,nc))
(ia<-seq(1,nm*nc-1,nm))
(ib<-seq(ia[2]-1,nm*nc,nm))
for(i in 1:nc){
  met[[i]]<-MET[ia[i]:ib[i],]
  mt[,,i]=as.matrix(MET[ia[i]:ib[i],])
}
mt=as.data.frame(apply(mt,c(1,2),mean))
names(mt)=names(MET)
head(mt)

#Gridded latitude/longitude##################################
LAT=read.table('lat.txt');
LON=read.table('lon.txt');
nrow=nrow(LAT);
ncol=nrow(LON)

#PARAMS#####################################################
par<-read.table("Free_Par.txt",header=TRUE)
lower<-par$lower
upper<-par$upper
params<-par$params


#allocate size of outputs
#gridded out put##########################################
PETg=array(0, c(nrow,ncol,nms));

#output for valid cells
#PET=zeros(nms,nc);SN=PET;RunOff=PET;ST=PET;S=PET;AET=PET;

#test run with given parameters########################################
RES=wbmmodel(params,nc,nm,nms,id,prec,tmpr,epot,D,H,lake,cellarea,petype,wbmtype,qobs,newtime,startsel);

#Optimise
x0=t(params);
bl=t(lower);
bu=t(upper);
maxn=10000;
kstop=5;
pcento=0.01;
peps=0.0001;
ngs=2*13;
iseed=1969;
iniflg=0;
NS<-wbmmodel(x0,nc,nm,nms,id,prec,tmpr,epot,D,H,lake,cellarea,petype,wbmtype,qobs,newtime,startsel);
#[bestx,bestf] = sceua(x0,bl,bu,maxn,kstop,pcento,peps,ngs,iseed,iniflg, nc,nm,nms,id,prec,tmpr,epot,D,H,lake,cellarea,petype,wbmtype,qobs,newtime,startsel);

optim(x0,fn="wbsmodel",gr=NULL,
	,nc,nm,nms,id,prec,tmpr,epot,D,H,lake,cellarea,petype,
	wbmtype,qobs,newtime,startsel,
	lower=bl,upper=bu,control(fnscale=-1))
