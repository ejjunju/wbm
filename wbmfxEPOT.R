fxEPOT=function(petype,tmpr,pet,d,h,rhum=0,wind=0){
	if (petype==1){
        # Hamon
		#Wt is a saturated water vapor density term
		Wt =((4.95/100 )* exp(0.062* tmpr));
		#Potential Evap
		pet=(13.97 * d * power((h/12),2) * Wt); #mm/month
    }
    
    if (petype==2){
        #Method2  Thornthwaite
		 pet=(d*2.1*power(h,2)*0.611*exp(17.3* tmp/( tmp+237.3))/( tmp+273.2));
		if (tmp<0 ){
            pet=0;
        }
    }
	
	if (petype==3){
        pet=pet;# From Datamm/d*days
    }
	
	if(petype==3){#penman monteith/priestly/Odin
		
	}
return(pet)
}