function [alt,pressure,sat,rho,tas,eas] = atmoshperic(flightdata, index, indexend)

%Pressure and density in metric units

alt = flightdata.Dadc1_alt.data(index:indexend);
pressure = 101325.*(288.15./(288.15-0.0019812.*(alt))).^((9.80665*0.0289644)/(8.3144598*-0.0065));
sat = flightdata.Dadc1_sat.data(index:indexend);
rho = pressure./(8.3144598./0.0289644.*(273.15+sat));

%Equivalent velocity of aircraft

tas = flightdata.Dadc1_tas.data(index:indexend)*0.5144447;
eas = tas.*sqrt(rho./1.225);

