thrust = [2200,5.7,3.2]; %thrust force, position x, position y

geopos = [[4.2, 2.7];[6.7,3.8];[285.56,0]]; %1st wing position x and y, 2nd horizontal stabiliser position x and y and fuel position in inches

geospecs = [30,2.0569]; %in metric, surface & mac

%Time point in minutes since the powerup of the recording of data

timeoffset = 19+50/60;
duration = 1;

timeoffsetcgshift = 31;
duration = 1;


bem=9165;
xcgbem=292.18;
fuelloaded=4050;
load payloadvals;

rampweight = [bem+fuelloaded+sum(payload(:,2)),xcgbem,3]; %OW, x cg pos, z cg pos in ibs and inches

%Unit conversion

rampweightmetric = [rampweight(1)*0.45359237, rampweight(2)*0.0254, rampweight(3)*0.0254];
geoposmetric = [[geopos(1,1)*0.0254, geopos(1,2)*0.0254];[geopos(2,1)*0.0254,geopos(2,2)*0.0254];[geopos(3,1)*0.0254, geopos(3,2)*0.0254]];

%Index of the start of the test

index = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + timeoffset*60);
index = index(1);
indexend = index + duration*112*60;

%Weight obtaining the weight and cg at particular time point

fuelused = (flightdata.lh_engine_FU.data + flightdata.rh_engine_FU.data)*0.45359237;
weightmetric = rampweightmetric(1) - fuelused(index:indexend);

%cg = [(rampweightmetric(1)*rampweightmetric(2)-fuelused(index)*geoposmetric(3,1))/weightmetric(1),0]
[ow,xcg,t] = cgcomp(bem,xcgbem,index,flightdata.lh_engine_FU.data(index),flightdata.rh_engine_FU.data(index),payload,fuelloaded);


%Velocity of the aircraft in m/s

tas = flightdata.Dadc1_tas.data(index:indexend)*0.5144447;

%Pressure and density in metric units

alt = flightdata.Dadc1_alt.data(index:indexend);
pressure = 101325.*(288.15./(288.15-0.0019812.*(alt))).^((9.80665*0.0289644)/(8.3144598*-0.0065));
sat = flightdata.Dadc1_sat.data(index:indexend);
rho = pressure./(8.3144598./0.0289644.*(273.15+sat));

%Equivalent velocity of aircraft

eas = tas.*sqrt(rho./1.225);



%Plot Elevator Trim Curve vs alpha

% plat = polyfit(flightdata.vane_AOA.data(index:indexend),flightdata.delta_e.data(index:indexend),1)
% scatter(flightdata.vane_AOA.data(index:indexend),flightdata.delta_e.data(index:indexend))
% axis([0 5 0 5],'ij')
% hold on
% xplt=[0:1:20];
% yplt=polyval(plat,xplt);
% %P = polyfit(flightdata.vane_AOA.data(index:indexend),flightdata.delta_e.data(index:indexend),1)
% plot(xplt,yplt)

%Plot Elevator Trim Curve vs tas
%Still have to pick the right data to approximate this
% blyat = polyfit(tas,flightdata.delta_e.data(index:indexend),2)
% xplt=[0:1:150];
% ypltav=polyval(blyat,xplt);
% scatter(tas,flightdata.delta_e.data(index:indexend))
% hold on
% plot(xplt,ypltav)
% axis ij

%Cmdeltae calculation from cg shift

%cmde=cmde(11000*4.448,145,1.,2,3.,30,2.)
%cmde(W,V,rho,deltae,deltacg,S,cbar)

%Cmalpha



