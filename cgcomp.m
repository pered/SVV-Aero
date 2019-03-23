function [ow,xcg,t] = cgcomp(bem,xcgbem,t,lfu,rfu,payload,fuelloaded)
    %xcgbem=292.18
    %luggage not included. Payload format as output of payloadfun.m
    %Typical: fuelloaded=4050, bem=9165
    momentarms=[131 131 214 214 251 251 288 288 170];
    
    payloadmoment=momentarms*payload(:,2)*2.20462;
    
    fuel=fuelloaded-lfu-rfu;
    fuelslope=2.8570048899755501222493887530562;
    fuelmoment = (fuelslope*(fuel-100)+298.16)*100; %inch-pounds
    fuelmomentinit = (fuelslope*(fuelloaded-100)+298.16)*100;
    
    ow=bem+sum(payload(:,2))*2.20462+fuel;
    
    %xcgs:
    zfm=bem+sum(payload(:,2))*2.20462
    xcgzfm=(xcgbem*bem+payloadmoment)/zfm
    rampm=zfm+fuelloaded
    xcgramp=(xcgzfm*zfm+fuelmomentinit)/rampm
    
    xcg=(bem*xcgbem+payloadmoment+fuelmoment)/ow