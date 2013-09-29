WAPABAwbm=function(Prain,pet,S,ST,A1,A2,B,K,STC,cellarea,d){
	##stage1:  Catchment rainfall retention
	if (Prain>0){
		Xot=max(0, (pet+STC-ST));
		Xt = Prain + Xot - power((power(Prain,A1) + power(Xot,A1)),(1/A1));
    }
    if (Prain<=0 ){
        Xt=0;
    }

	Yt=Prain-Xt;
	
	##Stage2 :
	#Water availability for ET 
	Wt=ST+Xt;
	##actual evapotranspiration
	if (Wt>0){
		AET = Wt + pet - power((power(Wt,A2) + power(pet,A2)),(1/A2));
    } 
    
    if (Wt<=0){
        AET=0;
    }
	#####//std::cout<<"AET["<<i<<","<<j<<"]="<<AET<<"\t";//For debugging							
	 ST=(Wt-AET); ##Water in the soil store at the } of time step t

	#Deep drainage (recharge) and direct runoff
	Yt=max(0,Yt); #Ensures that Yield (Yt=Prain-Xt) is not negative

	if (ST>STC){ #(to resolve condition when st > Smax)
		Yt=Yt-STC+ST;# Yield takes any extra water over STC
		ST=STC;
    }

	Rt=B*Yt; #recharge to the groundwater store
	DRO=max(0,(Yt-Rt));# #surface runoff
	
	##Stage4
	# #Groundwater store
	##The groundwater store is assumed to follow a simple linear storage-outflow relationship
	Texp=(1 - exp(-K));#=(1-e-T/K) (T/K) has been represented as One parameter K.
	RO=S*(Texp) +Rt*(1-Texp/K);#Base flow
	
	#Update Ground water storage
	 S=S+Rt-RO;
	###//std::cout<<"S["<<i<<","<<j<<"]="<< S<<"\t";#
	st=ST;
    s=S;
	##Stage5
	#RUNOFF
	RunOffmm=DRO+RO;# #total flow
	RunOff=RunOffmm/1000*cellarea/(d*24*3600);#monthly in m3/s
	out<-list(RunOff=RunOff,st=st,s=s,AET=AET)
return(out)
}
