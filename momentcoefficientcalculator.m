%Time point in minutes since the powerup of the recording of data


indexmeasurements = find(diff(find(flightdata.measurement_running.data))~=1);

indexes = find(flightdata.measurement_running.data);

start = [indexes(1),indexes(indexmeasurements(1));
    indexes(indexmeasurements(1:end-1)+1),indexes(indexmeasurements(2:end)); 
    indexes(indexmeasurements(end)+1), indexes(end)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geospecs = [30,2.0569]; %in metric, surface & mac

bem=9165;
xcgbem=292.18;
fuelloaded=4050;

%Index of the start of the test

indexcgshift1 = start(15,1):start(15,2);
indexcgshift2 = start(16,1):start(16,2);

indexspeedrun1 = start(8,1):start(8,2);
indexspeedrun2 = start(9,1):start(9,2);
indexspeedrun3 = start(10,1):start(10,2);
indexspeedrun4 = start(11,1):start(11,2);
indexspeedrun5 = start(12,1):start(12,2);
indexspeedrun6 = start(13,1):start(13,2);
indexspeedrun7 = start(14,1):start(14,2);

disp(['Elevator Trim Test 1: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun1(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun1(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 2: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun2(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun2(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 3: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun3(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun3(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 4: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun4(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun4(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 5: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun5(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun5(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 6: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun6(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun6(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 7: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun7(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun7(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['CG Shift Test 1: Start:', num2str((flightdata.Gps_utcSec.data(indexcgshift1(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexcgshift1(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['CG Shift Test 2: Start:', num2str((flightdata.Gps_utcSec.data(indexcgshift2(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexcgshift2(end))-flightdata.Gps_utcSec.data(1))/60)])


%atmospheric

cgshiftatmospheric = [];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexcgshift1(1), indexcgshift1(end));

cgshiftatmospheric = [cgshiftatmospheric; mean(rho), mean(tas)];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexcgshift2(1), indexcgshift2(end));

cgshiftatmospheric = [cgshiftatmospheric; mean(rho), mean(tas)];

%%%

speedrunatmospheric = [];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun1(1), indexspeedrun1(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun2(1), indexspeedrun2(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun3(1), indexspeedrun3(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun4(1), indexspeedrun4(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun5(1), indexspeedrun5(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun6(1), indexspeedrun6(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun7(1), indexspeedrun7(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas)]];


%Calculation of CG Shift Elevator Angle


deltae = [mean(flightdata.delta_e.data(indexcgshift1));
mean(flightdata.delta_e.data(indexcgshift2))];


%Plot Data


speedrunplot = [];

speedrunplot = [speedrunplot; mean(flightdata.vane_AOA.data(indexspeedrun1)), mean(flightdata.delta_e.data(indexspeedrun1)); 
    mean(flightdata.vane_AOA.data(indexspeedrun2)), mean(flightdata.delta_e.data(indexspeedrun2));
    mean(flightdata.vane_AOA.data(indexspeedrun3)), mean(flightdata.delta_e.data(indexspeedrun3));
    mean(flightdata.vane_AOA.data(indexspeedrun4)), mean(flightdata.delta_e.data(indexspeedrun4));
    mean(flightdata.vane_AOA.data(indexspeedrun5)), mean(flightdata.delta_e.data(indexspeedrun5));
    mean(flightdata.vane_AOA.data(indexspeedrun6)), mean(flightdata.delta_e.data(indexspeedrun6));
    mean(flightdata.vane_AOA.data(indexspeedrun7)), mean(flightdata.delta_e.data(indexspeedrun7));];


stickforcesplot = [diff(flightdata.column_fe.data(10000:20000))./diff(flightdata.Ahrs1_VertAcc.data(10000:20000)),diff(flightdata.delta_e.data(10000:20000))./diff(flightdata.Ahrs1_VertAcc.data((10000:20000)))];

%Matrix with values




% %Weight obtaining the weight and cg at particular time point
% 
% fuelused = (flightdata.lh_engine_FU.data + flightdata.rh_engine_FU.data)*0.45359237;
% weightmetric = rampweightmetric(1) - fuelused(index:indexend);
% 
% %cg = [(rampweightmetric(1)*rampweightmetric(2)-fuelused(index)*geoposmetric(3,1))/weightmetric(1),0]
% [ow,xcg,t] = cgcomp(bem,xcgbem,index,flightdata.lh_engine_FU.data(index),flightdata.rh_engine_FU.data(index),payload,fuelloaded);



%Plot Elevator Trim Curve vs alpha

figure(1)
plat = polyfit(speedrunplot(:,1),speedrunplot(:,2),1);
scatter(speedrunplot(:,1),speedrunplot(:,2))
axis('ij')
hold on
xplt=[0:1:13];
yplt=polyval(plat,xplt);
%P = polyfit(flightdata.vane_AOA.data(index:indexend),flightdata.delta_e.data(index:indexend),1)
plot(xplt,yplt)


%Plot Elevator Trim Curve vs tas

figure(2)
blyat = polyfit((speedrunatmospheric(:,2).^(-2)),speedrunplot(:,2),1);
xplt=[35:1:150];
ypltav=polyval(blyat,xplt.^(-2));
scatter(speedrunatmospheric(:,2),speedrunplot(:,2))
hold on
plot(xplt,ypltav)
axis ij



%[ow,xcg,t] = cgcomp(bem,xcgbem,index,flightdata.lh_engine_FU.data(index),flightdata.rh_engine_FU.data(index),payload,fuelloaded);

[owref1,xcgref1,t1] = cgcomp(bem,xcgbem,round(mean(indexcgshift1)),flightdata.lh_engine_FU.data(round(mean(indexcgshift1))),flightdata.rh_engine_FU.data(round(mean(indexcgshift1))),payloadref,fuelloaded);
[owref2,xcgref2,t2] = cgcomp(bem,xcgbem,round(mean(indexcgshift2)),flightdata.lh_engine_FU.data(round(mean(indexcgshift2))),flightdata.rh_engine_FU.data(round(mean(indexcgshift2))),payloadrefshifted,fuelloaded);
deltacg = xcgref2 - xcgref1;
owrefmean=(owref1+owref2)/2;

%cmde(W,V,rho,deltae,deltacg,S,cbar)

cmde = cmdee(owrefmean*4.44822,mean(cgshiftatmospheric(:,2)*0.5144),mean(cgshiftatmospheric(:,1)),diff(deltae)*pi()/180,deltacg*0.0254,geospecs(1),geospecs(2));

%Cmalpha
dealpha = plat(1);
cmalpha = cmde * (-1) * dealpha;

disp(['Cmalpha is: ',num2str(cmalpha), ' and Cmdeltae is:', num2str(cmde)])

