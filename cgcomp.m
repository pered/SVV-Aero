%function [ow,xcg,t] = cgcomp(bem,xcgbem,t,lfu,rfu,payload,fuelloaded)
%xcgbem=292.18
%ow=bem+payload
momentarms=[131 131 214 214 251 251 288 288 170];
payload=payloadfun();
fuel=fuelloaded-lfu-rfu;
fuelmoment = fuelslope*(
fuelslope=2.8570048899755501222493887530562