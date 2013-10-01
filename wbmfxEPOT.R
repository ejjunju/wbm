fxEPOT=function(petype,tmpr,pet,d,h,rhum=0,wind=0,tmon=0,tmax=0,lat){
	if (petype==1){        # Hamon
		#Wt is a saturated water vapor density term
		Wt =((4.95/100 )* exp(0.062* tmpr));
		pet=(13.97 * d * power((h/12),2) * Wt); #mm/month
    }
    
    if (petype==2){#Method2  Thornthwaite
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


#library(SPEI)
# 
# thornthwaite(Tave, lat, na.rm = FALSE)
# # hargreaves(Tmin, Tmax, Ra = NA, lat = NA, Pre = NA, na.rm = FALSE)
# # penman(Tmin, Tmax, U2, Ra = NA, lat = NA, Rs = NA, tsun = NA,
#     CC = NA, ed = NA, Tdew = NA, RH = NA, P = NA, P0 = NA,
#     z = NA, crop='short', na.rm = FALSE)