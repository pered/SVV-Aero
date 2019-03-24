%%% Variables %%%

firstspeedrunindex = 8;
firstcgshiftindex = 15;

Ws = 60500;
mdotfs = 0.048;
cmtc = -0.0064;

geospecs = [30,2.0569]; %in metric, surface & mac
bem=9165;
xcgbem=292.18;
fuelloaded=4050;


%%% Detection of Measurement Periods %%%


indexmeasurements = find(diff(find(flightdata.measurement_running.data))~=1);

indexes = find(flightdata.measurement_running.data);

start = [indexes(1),indexes(indexmeasurements(1));
    indexes(indexmeasurements(1:end-1)+1),indexes(indexmeasurements(2:end)); 
    indexes(indexmeasurements(end)+1), indexes(end)];


%%% Indexing of the Tests %%%

indexcgshift1 = start(firstcgshiftindex,1):start(firstcgshiftindex,2);
indexcgshift2 = start(firstcgshiftindex+1,1):start(firstcgshiftindex+1,2);

indexspeedrun1 = start(firstspeedrunindex,1):start(firstspeedrunindex,2);
indexspeedrun2 = start(firstspeedrunindex+1,1):start(firstspeedrunindex+1,2);
indexspeedrun3 = start(firstspeedrunindex+2,1):start(firstspeedrunindex+2,2);
indexspeedrun4 = start(firstspeedrunindex+3,1):start(firstspeedrunindex+3,2);
indexspeedrun5 = start(firstspeedrunindex+4,1):start(firstspeedrunindex+4,2);
indexspeedrun6 = start(firstspeedrunindex+5,1):start(firstspeedrunindex+5,2);
indexspeedrun7 = start(firstspeedrunindex+6,1):start(firstspeedrunindex+6,2);

disp(['Elevator Trim Test 1: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun1(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun1(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 2: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun2(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun2(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 3: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun3(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun3(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 4: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun4(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun4(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 5: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun5(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun5(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 6: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun6(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun6(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['Elevator Trim Test 7: Start:', num2str((flightdata.Gps_utcSec.data(indexspeedrun7(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexspeedrun7(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['CG Shift Test 1: Start:', num2str((flightdata.Gps_utcSec.data(indexcgshift1(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexcgshift1(end))-flightdata.Gps_utcSec.data(1))/60)])
disp(['CG Shift Test 2: Start:', num2str((flightdata.Gps_utcSec.data(indexcgshift2(1))-flightdata.Gps_utcSec.data(1))/60) ,' End:',num2str((flightdata.Gps_utcSec.data(indexcgshift2(end))-flightdata.Gps_utcSec.data(1))/60)])


%%% Atmospheric Parameters %%%

cgshiftatmospheric = [];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexcgshift1(1), indexcgshift1(end));

cgshiftatmospheric = [cgshiftatmospheric; mean(rho), mean(tas)];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexcgshift2(1), indexcgshift2(end));

cgshiftatmospheric = [cgshiftatmospheric; mean(rho), mean(tas)];

speedrunatmospheric = [];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun1(1), indexspeedrun1(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas), mean(eas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun2(1), indexspeedrun2(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas), mean(eas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun3(1), indexspeedrun3(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas), mean(eas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun4(1), indexspeedrun4(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas), mean(eas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun5(1), indexspeedrun5(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas), mean(eas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun6(1), indexspeedrun6(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas), mean(eas)]];

[alt,pressure,sat,rho,tas,eas] = atmospheric(flightdata, indexspeedrun7(1), indexspeedrun7(end));
speedrunatmospheric = [speedrunatmospheric;[mean(rho), mean(tas), mean(eas)]];


%%% Plot Data %%%

speedrunplot = [];

speedrunplot = [speedrunplot; mean(flightdata.vane_AOA.data(indexspeedrun1)), mean(flightdata.delta_e.data(indexspeedrun1)); 
    mean(flightdata.vane_AOA.data(indexspeedrun2)), mean(flightdata.delta_e.data(indexspeedrun2));
    mean(flightdata.vane_AOA.data(indexspeedrun3)), mean(flightdata.delta_e.data(indexspeedrun3));
    mean(flightdata.vane_AOA.data(indexspeedrun4)), mean(flightdata.delta_e.data(indexspeedrun4));
    mean(flightdata.vane_AOA.data(indexspeedrun5)), mean(flightdata.delta_e.data(indexspeedrun5));
    mean(flightdata.vane_AOA.data(indexspeedrun6)), mean(flightdata.delta_e.data(indexspeedrun6));
    mean(flightdata.vane_AOA.data(indexspeedrun7)), mean(flightdata.delta_e.data(indexspeedrun7));];

%%% dDeltae/dalpha best fit %%%

plat = polyfit(speedrunplot(:,1),speedrunplot(:,2),1);


%%% Calculation of Cmde and Cmalpha %%%

deltae = [mean(flightdata.delta_e.data(indexcgshift1));
mean(flightdata.delta_e.data(indexcgshift2))];

[owref1,xcgref1] = cgcomp(bem,xcgbem,indexcgshift1,flightdata.lh_engine_FU.data,flightdata.rh_engine_FU.data,payloadactual,fuelloaded);
[owref2,xcgref2] = cgcomp(bem,xcgbem,indexcgshift2,flightdata.lh_engine_FU.data,flightdata.rh_engine_FU.data,payloadactualshifted,fuelloaded);

cmde = cmdee((owref1+owref2)/2*4.44822,mean(cgshiftatmospheric(:,2)),mean(cgshiftatmospheric(:,1)),diff(deltae)*pi()/180,(xcgref2 - xcgref1)*0.0254,geospecs(1),geospecs(2));

%
dealpha = plat(1);
cmalpha = cmde * (-1) * dealpha;

disp(['Cmalpha is: ',num2str(cmalpha), ' and Cmdeltae is:', num2str(cmde)])


%%% Reductions %%%

%Plot reduced Elevator Trim Curve

owlist = [cgcomp(bem,xcgbem,indexspeedrun1,flightdata.lh_engine_FU.data,flightdata.rh_engine_FU.data,payloadactual,fuelloaded);
    cgcomp(bem,xcgbem,indexspeedrun2,flightdata.lh_engine_FU.data,flightdata.rh_engine_FU.data,payloadactual,fuelloaded);
    cgcomp(bem,xcgbem,indexspeedrun3,flightdata.lh_engine_FU.data,flightdata.rh_engine_FU.data,payloadactual,fuelloaded);
    cgcomp(bem,xcgbem,indexspeedrun4,flightdata.lh_engine_FU.data,flightdata.rh_engine_FU.data,payloadactual,fuelloaded);
	cgcomp(bem,xcgbem,indexspeedrun5,flightdata.lh_engine_FU.data,flightdata.rh_engine_FU.data,payloadactual,fuelloaded);
	cgcomp(bem,xcgbem,indexspeedrun6,flightdata.lh_engine_FU.data,flightdata.rh_engine_FU.data,payloadactual,fuelloaded);
	cgcomp(bem,xcgbem,indexspeedrun7,flightdata.lh_engine_FU.data,flightdata.rh_engine_FU.data,payloadactual,fuelloaded)];

speedrunatmosphericreduced = [speedrunatmospheric(:,1), speedrunatmospheric(:,2)./sqrt(owlist(:,1)*4.44822).*sqrt(Ws) , speedrunatmospheric(:,3)];


%Thrust Reductions

atmospherethrust = [mean(flightdata.Dadc1_alt.data(indexspeedrun1))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun1)), mean(flightdata.Dadc1_tat.data(indexspeedrun1)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun1))),mean(flightdata.lh_engine_FMF.data(indexspeedrun1))*0.4535923666/3600,mean(flightdata.rh_engine_FMF.data(indexspeedrun1))*0.4535923666/3600;
    mean(flightdata.Dadc1_alt.data(indexspeedrun2))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun2)), mean(flightdata.Dadc1_tat.data(indexspeedrun2)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun2))),mean(flightdata.lh_engine_FMF.data(indexspeedrun2))*0.4535923666/3600,mean(flightdata.rh_engine_FMF.data(indexspeedrun2))*0.4535923666/3600;
    mean(flightdata.Dadc1_alt.data(indexspeedrun3))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun3)), mean(flightdata.Dadc1_tat.data(indexspeedrun3)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun3))),mean(flightdata.lh_engine_FMF.data(indexspeedrun3))*0.4535923666/3600,mean(flightdata.rh_engine_FMF.data(indexspeedrun3))*0.4535923666/3600;
    mean(flightdata.Dadc1_alt.data(indexspeedrun4))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun4)), mean(flightdata.Dadc1_tat.data(indexspeedrun4)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun4))),mean(flightdata.lh_engine_FMF.data(indexspeedrun4))*0.4535923666/3600,mean(flightdata.rh_engine_FMF.data(indexspeedrun4))*0.4535923666/3600;
    mean(flightdata.Dadc1_alt.data(indexspeedrun5))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun5)), mean(flightdata.Dadc1_tat.data(indexspeedrun5)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun5))),mean(flightdata.lh_engine_FMF.data(indexspeedrun5))*0.4535923666/3600,mean(flightdata.rh_engine_FMF.data(indexspeedrun5))*0.4535923666/3600;
    mean(flightdata.Dadc1_alt.data(indexspeedrun6))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun6)), mean(flightdata.Dadc1_tat.data(indexspeedrun6)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun6))),mean(flightdata.lh_engine_FMF.data(indexspeedrun6))*0.4535923666/3600,mean(flightdata.rh_engine_FMF.data(indexspeedrun6))*0.4535923666/3600;
    mean(flightdata.Dadc1_alt.data(indexspeedrun7))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun7)), mean(flightdata.Dadc1_tat.data(indexspeedrun7)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun7))),mean(flightdata.lh_engine_FMF.data(indexspeedrun7))*0.4535923666/3600,mean(flightdata.rh_engine_FMF.data(indexspeedrun7))*0.4535923666/3600];

atmospherethruststandardised = [mean(flightdata.Dadc1_alt.data(indexspeedrun1))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun1)), mean(flightdata.Dadc1_tat.data(indexspeedrun1)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun1))),mdotfs,mdotfs;
    mean(flightdata.Dadc1_alt.data(indexspeedrun2))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun2)), mean(flightdata.Dadc1_tat.data(indexspeedrun2)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun2))),mdotfs,mdotfs;
    mean(flightdata.Dadc1_alt.data(indexspeedrun3))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun3)), mean(flightdata.Dadc1_tat.data(indexspeedrun3)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun3))),mdotfs,mdotfs;
    mean(flightdata.Dadc1_alt.data(indexspeedrun4))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun4)), mean(flightdata.Dadc1_tat.data(indexspeedrun4)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun4))),mdotfs,mdotfs;
    mean(flightdata.Dadc1_alt.data(indexspeedrun5))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun5)), mean(flightdata.Dadc1_tat.data(indexspeedrun5)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun5))),mdotfs,mdotfs;
    mean(flightdata.Dadc1_alt.data(indexspeedrun6))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun6)), mean(flightdata.Dadc1_tat.data(indexspeedrun6)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun6))),mdotfs,mdotfs;
    mean(flightdata.Dadc1_alt.data(indexspeedrun7))*0.3048,mean(flightdata.Dadc1_mach.data(indexspeedrun7)), mean(flightdata.Dadc1_tat.data(indexspeedrun7)-(15-1.98./1000.*flightdata.Dadc1_alt.data(indexspeedrun7))),mdotfs,mdotfs];

dlmwrite('matlab.dat', atmospherethrust, ' ');

fileexe_path = which ('thrust.exe');
system ([fileexe_path]);

thrust = sum(load('thrust.dat'),2);

dlmwrite('matlab.dat', atmospherethruststandardised, ' ');

fileexe_path = which ('thrust.exe');
system ([fileexe_path]);

thruststandardised = sum(load('thrust.dat'),2);

thrustcoefactual = thrust./(0.5.*speedrunatmospheric(:,1).*speedrunatmospheric(:,2).^2.*geospecs(1));

thrustcoefstandardised = thruststandardised./(0.5.*speedrunatmospheric(:,1).*speedrunatmospheric(:,2).^2.*geospecs(1));

deltaereduced = speedrunplot(:,2)-1./cmde.*cmtc.*(thrustcoefstandardised-thrustcoefactual);


%Reduced Elevator Control Force Curve 


column_fe = [mean(flightdata.column_fe.data(indexspeedrun1));
    mean(flightdata.column_fe.data(indexspeedrun2));
    mean(flightdata.column_fe.data(indexspeedrun3));
    mean(flightdata.column_fe.data(indexspeedrun4));
    mean(flightdata.column_fe.data(indexspeedrun5));
    mean(flightdata.column_fe.data(indexspeedrun6));
    mean(flightdata.column_fe.data(indexspeedrun7))];

column_fereduced = column_fe.*Ws./(owlist*4.44822);

%%% Extra Plots %%%

%Plot Elevator Trim Curve vs tas

figure(1)
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
scatter(speedrunatmosphericreduced(:,2),speedrunplot(:,2))
plot(xplt,ypltav)
axis ij

clc





