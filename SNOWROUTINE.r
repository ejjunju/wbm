functionSNOWROUTINE=function(prc,tmp,snostor,Tsnow,Train,meltmax){#SNOW ACCUMULATION


	if (tmp<Tsnow){
		Psnow= prc; 
	}
	if (tmp>=Tsnow){
		if((tmp>=Tsnow) && ( tmp<=Train)){
			
			Psnow= prc*(Train- tmp)/(Train-Tsnow);
		}
		if (tmp>Train){  
			Psnow=0;
		}
	}
	
	Prain= prc-Psnow;
	SMF=( tmp-Tsnow)/(Train-Tsnow)*meltmax;
	
	if(SMF>meltmax){
		SMF=meltmax;
	}
	SM = snostor*SMF; ##Snow melt
	snostor= snostor-SM;##Update state
	
	if( snostor<0){ 
		snostor=0;##Update state
	}
	
	Prain=Prain+SM; #melt + remain (Sent to next routine)
out=list(Prain=Prain,snostor=snostor)
return(out)
}
