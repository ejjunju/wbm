summ=function(inp,newtime,typ){#sum/or mean data array along (t)
res =apply(inp,3,mean);# mean along time
if (typ==1){ res=sum(inp,2);}
res=zoo(res,newtime)
}

zeros=function(x,y){z=matrix(0,x,y);return(z)} #matrix of zeros

wbmmodel=function(params,nc,nm,nms,id,prec,tmpr,epot,D,H,lake,cellarea,petype,wbmtype,qobs,newtime,startDate){
#params
Train=params(1);Tsnow=params(2);drofrac=params(3);meltmax=params(4);
STC=params(5);rfactor=params(6);A1=params(7);A2=params(8);B=params(9);
K=params(10);SNo=params(11);STo=params(12);So=params(13);

PET=zeros(nms,nc);SN=PET;RunOff=PET;ST=PET;S=PET;AET=PET;

for (i in 1:nc){ # For each grid cell
    #get data for cell
    #disp(i)
    idi1=(nm*i-nm+1);idi2=(nm*i); idi=idi1:idi2; idi=t(idi);
    prci=prec[idi];
    tmpi=tmpr[idi];
    peti=epot[idi];
    di=D[idi];
    hi=H[idi];
    
    
    #subset
    prci=prci[id];#prci=timeseries[prci,newtime);
    tmpi=tmpi[id];#tmpi=timeseries[tmpi, id);
    peti=peti[id];#peti=timeseries[peti, id);
    di=di[id];#di=timeseries[di,id);
    hi=hi[id];#hi=timeseries[hi,id);
    
    #fixed
    Lake=lake[i];
    CArea=cellarea[i];
    
    #Loop through time 
    for (t in 1:nms){
        #initial States
        if(t==1){
            snostor=SNo;
            ST=STo;
            S=So;
        }
        if(t>1){
            snostor=SN((t-1),i);
            ST=ST((t-1),i);
            S=S((t-1),i);
        }
        #data
        prc=prci[t];tmp=tmpi[t];pet=peti[t];d=di[t];h=hi[t];
        #Potential Evaporation
        pot=fxEPOT(petype,tmp,pet,di,h);
        #Snow routine
        S3=SNOWROUTINE(prc,tmp,snostor,Tsnow,Train,meltmax);
        Prain=S3$Prain
        snostor=S3$snostor
        #wbmmodel
        if (wbmtype==3){
            #wapaba
        	W3=WAPABAwbm(Prain,pot,S,ST,A1,A2,B,K,STC,CArea,d);
        	ro<-W3$ro;st<-W3$st;s=W3$s;aet=W3$aet;
        }
        if (wbmtype!=3){
            #thornthwaite
            W3=THORNTHWAITEwbm(Prain,pot,S,ST,drofrac,rfactor,STC,CArea,d,Lake);
        	ro<-W3$ro;st<-W3$st;s=W3$s;aet=W3$aet;
        } 
        
        #To be used in summary
        SN[t,i]=snostor;
        PET[t,i]=pot;
        RunOff[t,i]=ro;
        ST[t,i]=st;
        S[t,i]=s;
        AET[t,i]=aet;        
    }
}

##summary for all grids
PETts =summ(PET,newtime,0);
AETts =summ(AET,newtime,0);
STts =summ(ST,newtime,0);
Sts =summ(S,newtime,0);#
QSim =summ(RunOff,newtime,1);
QObs =timeseries(qobs,newtime);


qsim=sum(RunOff,2);
T=1:nms;T=t(T);

library(topmodel)
NS=NSeff(qobs,qsim);
out=list(NS = NS,PETts = PETts,AETts = AETts,STts = STts,Sts = Sts,QObs = QObs,QSim = QSim,ST = ST,S = S,PET = PET,AET = AET,RunOff = RunOff)
if(optim==TRUE){return(NS)} else {return(OUT)}
}