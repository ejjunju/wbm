%load('testData.mat') Load data if it exists
cd 'C:\DATA\STATIONS\wbm_Rusumo\Calib\matlab'

%DATA IMPORT
%Control/initialising 
[nm,nc,petype,wbmtype] = importInit('init.txt');%Time & Qobs
[Time,Qobs] = importTimeQobs('Qobs.txt');
[prec,tmpr,epot,d,h] = importMet('met.txt');%met data (nrows=nm*nC)[prec,tmpr,epot,d,h]
[ir,ic,z,lake,cellarea,Lat] = importFxd('fxd.txt');%Fixed
[cp,params,lower1,upper1,FLG] = importPar('Free_par.txt');%Parameters

Qobs=timeseries(Qobs,Time);
Qobs.TimeInfo.Format = 'mmm dd, yy';  % Set format for display on x-axis.
plot(Qobs);