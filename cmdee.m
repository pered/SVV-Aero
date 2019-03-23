function cmde = cmdee(W,V,rho,deltae,deltacg,S,cbar)
    cmde = W/(0.5*rho*V*V*S)/(-deltae)*deltacg/cbar;
end