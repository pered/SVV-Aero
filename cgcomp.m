function [ow,xcg,t] = cgcomp(bem,xcgbem,t,lfu,rfu,payload,fuelloaded)
    %xcgbem=292.18
    %luggage not included
    momentarms=[131 131 214 214 251 251 288 288 170];
    %payload=payloadfun();
    payloadmoment=momentarms*payload(:,2);
    
    fuel=fuelloaded-lfu-rfu;
    fuelslope=2.8570048899755501222493887530562;
    fuelmoment = (fuelslope*(fuel-100)+298.16)*100; %inch-pounds
    
    
    ow=bem+sum(payload(:,2))+fuel;
    
    %xcgs:
    zfm=bem+sum(payload(:,2))
    xcgzfm=(xcgbem*bem+payloadmoment)/zfm
    rampm=zfm+fuelloaded
    xcgramp=(xcgzfm*zfm+fuelmoment)/rampm
    
    xcg=(bem*xcgbem+payloadmoment+fuelmoment)/ow