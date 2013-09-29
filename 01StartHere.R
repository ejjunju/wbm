#SOURCE THIS FILE
#Find the script directory
(script.dir <- dirname(sys.frame(1)$ofile))
print(paste("script.dir <-",script.dir))
setwd(script.dir)

#Load general functions
source("fx.gen.r")
graphics.off();
#Choose Modelling exercise
print("ENTER 1=wbm; 2=Evap; 3=hydroMapp; 4=nMAG\n")
#opts=scan(n=1)
opts<-1

#WBM
if(opts==1){print("WBM")
	#BASIN
##########################################################################################################
# ID GRIDCODE BAS catch            Name          Region US done        Area
# 28       21  21     0        Turkwell            <NA>  0    0  5951199954
# 32       49  49     0          Muzizi            <NA>  0    0  3617074775
# 43       20  20     0     Sondo Miriu            <NA>  0    0  3361406991
# 53       22  22     0          Rusumo            <NA>  0    0 30622842529
# 55       40  40     0 Nyumba ya Mungu            <NA>  0    0 13386599626
# 56       31  31     1 Pangani Falls 1 Pangani Falls 2  0    0    23857322
# 57       38  38     1 Pangani Falls 2            Hale  0    0    82193049
# 58       39  39     1            Hale Nyumba ya Mungu  0    0 25522520293
# 62       37  37     0           Mtera            <NA>  0    0 69633727012
	#
	cat<-readShapeSpatial("GIS/EastAfrica/catchments.shp",proj4string=pj4s)
	(basin<-cat[cat@data$ID =="53",]);
	(cname=as.character(basin@data$Name))
	(Dxy<-0.5)#resolution
	flofile="Rusumo.flo.rdata"
	(flofile=paste("C:\\DATA\\HYDROMET\\FLOW\\EAPMP\\",flofile,sep=""))
	source("wbmDataGridded.R") #create data for wbm model
	source("wbmCalib.R") #calibrate Model
	#Discuss fit
	source("wbmSim20c3m.R") #Simulate for 1961-1990 with GCM
	#Discuss error comapred to Actual data
	source("wbmPredictGCM.R") #PREDICT FOR CLIMATE CHNAGE
	#SOURCE("WBMfdr.R") #fLOW DURATION ANALYSIS
	#NMAG PREPARE
	#HYDROCALC
	#CHANGE ANALYSIS
	
	

	
}