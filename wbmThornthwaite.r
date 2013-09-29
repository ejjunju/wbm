THORNTHWAITEwbm=function(Prain,pet,S,ST,drofrac,rfactor,STC,cellarea,d,Lake){

	#DIRECT RUNOFF
	DRO=drofrac*Prain;
	Prain=Prain-DRO;
	##SOIL MOISTURE AND AET
	##STW (withdrawable water FROM SOIL)
	STW =  ST * (1 - abs((Prain- pet)/STC));
	##ConditionS FOR EVAPORATION AND storage
	##SOIL CELL
	if( Lake==0 ){
	
		if( Prain<  pet){  
            AET=Prain+STW; ST= ST-STW; 
        }
		if( Prain>= pet){   
			AET= pet;
			 ST= ST+Prain- pet; 
		}
		if( ST>STC){
			 S= S+ ST-STC;
			 ST=STC;
        }
	 }

	####std::cout<<"S["<<i<<","<<j<<"]="<< S<<"\t";#
	##LAKE CELL
	##LAKE EVAPORATION IN A DISTRIupperTED MODEL (AET= pet)
	if(Lake==1){
	
		AET= pet;
		 ST= ST+Prain- pet;
		 S= S+Prain- pet;
		if(  S<0) { 
            S=0;
        }##Lake surplus dep}s on Lake depth (use 40m as for Lake vic)
    }
    st=ST;
    s=S;
	#-----------------------------------------------------------
	##RUNOFF
	RO=rfactor* S;
	S= S-RO; ##next T-step
	##std::cout<<"RO["<<t<<"]="<< RO<<"\t";#
	##std::cout<<"S["<<t<<"]="<< S<<"\n";#
	##RUNOFF
	RunOffmm= DRO+RO;
	RunOff=RunOffmm/1000*cellarea/(d*24*3600);#monthly in m3/s
out<-list(RunOff=RunOff,st=st,s=s,AET=AET)
return(out)
}