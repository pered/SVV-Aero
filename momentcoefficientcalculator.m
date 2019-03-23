%Time point in minutes since the powerup of the recording of data


indexmeasurements = find(diff(find(flightdata.measurement_running.data))~=1)

indexes = find(flightdata.measurement_running.data)

start = [indexes(1),indexes(indexmeasurements(1));
    indexes(indexmeasurements(1:end-1)+1),indexes(indexmeasurements(2:end)); 
    indexes(indexmeasurements(end)+1), indexes(end)]

%%%CG Shift Test%%%


timeoffcgshift = [51 + 02/60, 5/60;
    52+46/60, 5/60]; %End time and duration


%%%Elevator Trim Tests%%%

%Speedrun 1-7

timeoffspeedrun = [37+19/60, 10/60; %End time and duration
    39+11/60,10/60;
    41+24/60,10/60;
    42+56/60,10/60;
    45+41/60,10/60;
    47+20/60,10/60;
    48+40/60,10/60];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geopos = [[4.2, 2.7];[6.7,3.8];[285.56,0]]; %1st wing position x and y, 2nd horizontal stabiliser position x and y and fuel position in inches

geospecs = [30,2.0569]; %in metric, surface & mac

bem=9165;
xcgbem=292.18;
fuelloaded=4050;
load payloadvals;

rampweight = [bem+fuelloaded+sum(payload(:,2)),xcgbem,3]; %OW, x cg pos, z cg pos in ibs and inches

%Unit conversion

rampweightmetric = [rampweight(1)*0.45359237, rampweight(2)*0.0254, rampweight(3)*0.0254];
geoposmetric = [[geopos(1,1)*0.0254, geopos(1,2)*0.0254];[geopos(2,1)*0.0254,geopos(2,2)*0.0254];[geopos(3,1)*0.0254, geopos(3,2)*0.0254]];



%Index of the start of the test

indexcgshift1 = start(17,1):start(17,2);
indexcgshift2 = start(18,1):start(18,2);

indexspeedrun1 = start(9,1):start(9,2);
indexspeedrun2 = start(10,1):start(10,2);
indexspeedrun3 = start(11,1):start(11,2);
indexspeedrun4 = start(12,1):start(12,2);
indexspeedrun5 = start(13,1):start(13,2);
indexspeedrun6 = start(14,1):start(14,2);
indexspeedrun7 = start(15,1):start(15,2);

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
plat = polyfit(speedrunplot(:,1),speedrunplot(:,2),1)
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

[owref1,xcgref1,t1] = cgcomp(bem,xcgbem,mean(round(mean(indexcgshift1))),flightdata.lh_engine_FU.data(round(mean(indexcgshift1))),flightdata.rh_engine_FU.data(round(mean(indexcgshift1))),payloadref,fuelloaded);
[owref2,xcgref2,t2] = cgcomp(bem,xcgbem,round(mean(indexcgshift2)),flightdata.lh_engine_FU.data(round(mean(indexcgshift2))),flightdata.rh_engine_FU.data(round(mean(indexcgshift2))),payloadrefshifted,fuelloaded);
deltacg = xcgref2 - xcgref1;
owrefmean=(owref1+owref2)/2;

%cmde(W,V,rho,deltae,deltacg,S,cbar)

cmde = cmdee(owrefmean*4.44822,mean(cgshiftatmospheric(:,2)*0.5144),mean(cgshiftatmospheric(:,1)),diff(deltae)*pi()/180,deltacg*0.0254,geospecs(1),geospecs(2))


%Cmalpha

dealpha = plat(1);
cmalpha = cmde * (-1) * dealpha

