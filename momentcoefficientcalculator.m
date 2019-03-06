thrust = [2200,5.7,3.2]; %thrust force, position x, position y

rampweight = [12000,292,3]; %OW, x cg pos, z cg pos

coefficients = [0.21;0.04]; %1st Cl, 2nd Cd

geopos = [[4.2, 2.7];[6.7,3.8];[285.56,0]]; %1st wing position x and y, 2nd horizontal stabiliser position x and y and fuel position

geospecs = [30,2.0569];

rho = 0.885;

v = 100;


%Time point in minutes since the powerup of the recording of data

timeoffset = 30;

%Index of the start of the test

index = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(1) + timeoffset*60);
index = index(1)

%Weight obtaining the weight and cg at particular time point

fuelused = (flightdata.lh_engine_FU.data + flightdata.rh_engine_FU.data);
weight = rampweight(1) - fuelused(index)

cg = [(rampweight(1)*rampweight(2)-fuelused(index)*geopos(3,1))/weight(1),0]

%Velocity of the aircraft

tas = flightdata.Dadc1_tas












