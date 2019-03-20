%Time point in minutes since the powerup of the recording of data

%%%CG Shift Test%%%


timeoffcgshift = [51 + 02/60, 5/60;
    52+46/60, 5/60]; %End time and duration


%%%Elevator Trim Tests%%%

%Speedrun 1-7

timeoffspeedrun = [37+19/60, 10/60; %End time and duration
    0,0;
    0,0;
    0,0;
    0,0;
    0,0;
    0,0];


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

indexcgshift1 = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + (timeoffcgshift(1,1)-timeoffcgshift(1,2))*60 & flightdata.Gps_utcSec.data <= flightdata.Gps_utcSec.data(1) + timeoffcgshift(1,1)*60);
indexcgshift2 = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + (timeoffcgshift(2,1)-timeoffcgshift(2,2))*60 & flightdata.Gps_utcSec.data <= flightdata.Gps_utcSec.data(1) + timeoffcgshift(2,1)*60);

indexspeedrun1 = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(1,1)-timeoffspeedrun(1,2)).*60 & flightdata.Gps_utcSec.data <= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(1,1)).*60);
indexspeedrun2 = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(2,1)-timeoffspeedrun(2,2)).*60 & flightdata.Gps_utcSec.data <= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(2,1)).*60);
indexspeedrun3 = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(3,1)-timeoffspeedrun(3,2)).*60 & flightdata.Gps_utcSec.data <= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(3,1)).*60);
indexspeedrun4 = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(4,1)-timeoffspeedrun(4,2)).*60 & flightdata.Gps_utcSec.data <= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(4,1)).*60);
indexspeedrun5 = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(5,1)-timeoffspeedrun(5,2)).*60 & flightdata.Gps_utcSec.data <= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(5,1)).*60);
indexspeedrun6 = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(6,1)-timeoffspeedrun(6,2)).*60 & flightdata.Gps_utcSec.data <= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(6,1)).*60);
indexspeedrun7 = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(7,1)-timeoffspeedrun(7,2)).*60 & flightdata.Gps_utcSec.data <= flightdata.Gps_utcSec.data(1) + (timeoffspeedrun(7,1)).*60);

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


%Matrix with values




% %Weight obtaining the weight and cg at particular time point
% 
% fuelused = (flightdata.lh_engine_FU.data + flightdata.rh_engine_FU.data)*0.45359237;
% weightmetric = rampweightmetric(1) - fuelused(index:indexend);
% 
% %cg = [(rampweightmetric(1)*rampweightmetric(2)-fuelused(index)*geoposmetric(3,1))/weightmetric(1),0]
% [ow,xcg,t] = cgcomp(bem,xcgbem,index,flightdata.lh_engine_FU.data(index),flightdata.rh_engine_FU.data(index),payload,fuelloaded);
% 
% 
% 
% 
% 
% %Plot Elevator Trim Curve vs alpha
% 
% plat = polyfit(flightdata.vane_AOA.data(index:indexend),flightdata.delta_e.data(index:indexend),1)
% scatter(flightdata.vane_AOA.data(index:indexend),flightdata.delta_e.data(index:indexend))
% axis([0 5 0 5],'ij')
% hold on
% xplt=[0:1:20];
% yplt=polyval(plat,xplt);
% %P = polyfit(flightdata.vane_AOA.data(index:indexend),flightdata.delta_e.data(index:indexend),1)
% plot(xplt,yplt)
% 
% %Plot Elevator Trim Curve vs tas
% %Still have to pick the right data to approximate this
% % blyat = polyfit(tas,flightdata.delta_e.data(index:indexend),2)
% % xplt=[0:1:150];
% % ypltav=polyval(blyat,xplt);
% % scatter(tas,flightdata.delta_e.data(index:indexend))
% % hold on
% % plot(xplt,ypltav)
% % axis ij
% 


%Cmdeltae calculation from cg shift

cmde = cmdee(11000*4.448,mean(cgshiftatmospheric(:,2)),mean(cgshiftatmospheric(:,1)),diff(deltae),3.,geospecs(1),geospecs(2))

indexcgshift1 = 30531;
indexcgshift2 = 31571;

[ow,xcg,t] = cgcomp(bem,xcgbem,index,flightdata.lh_engine_FU.data(index),flightdata.rh_engine_FU.data(index),payload,fuelloaded);
cmde = cmdee(11000*4.448,145,1.,2,3.,30,2.)

%cmde(W,V,rho,deltae,deltacg,S,cbar)
% 
% %Cmalpha
% 
% dealpha = plat(1);
% cmalpha = cmde * (-1) * dealpha 

