thrust = [2200,5.7,3.2]; %thrust force, position x, position y

rampweight = [12000,292,3]; %OW, x cg pos, z cg pos in ibs and inches

coefficients = [0.21;0.04]; %1st Cl, 2nd Cd

geopos = [[4.2, 2.7];[6.7,3.8];[285.56,0]]; %1st wing position x and y, 2nd horizontal stabiliser position x and y and fuel position in inches

geospecs = [30,2.0569]; %in metric


%Unit conversion

rampweightmetric = [rampweight(1)*0.45359237, rampweight(2)*0.0254, rampweight(3)*0.0254]
geoposmetric = [[geopos(1,1)*0.0254, geopos(1,2)*0.0254];[geopos(2,1)*0.0254,geopos(2,2)*0.0254];[geopos(3,1)*0.0254, geopos(3,2)*0.0254]]

%Time point in minutes since the powerup of the recording of data

timeoffset = 30;

%Index of the start of the test

index = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + timeoffset*60);
index = index(1)

%Weight obtaining the weight and cg at particular time point

fuelused = (flightdata.lh_engine_FU.data + flightdata.rh_engine_FU.data)*0.45359237;
weightmetric = rampweightmetric(1) - fuelused(index)

cg = [(rampweightmetric(1)*rampweightmetric(2)-fuelused(index)*geoposmetric(3,1))/weightmetric(1),0]

%Velocity of the aircraft in knots

tas = flightdata.Dadc1_tas.data(index)*0.5144447

%Pressure and density in metric units

alt = flightdata.Dadc1_alt.data(index)
pressure = 101325*(288.15/(288.15-0.0019812*(alt)))^((9.80665*0.0289644)/(8.3144598*-0.0065))
sat = flightdata.Dadc1_sat.data(index)
rho = pressure/(8.3144598/0.0289644*(273.15+sat))

%













