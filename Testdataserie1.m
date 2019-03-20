% Inputs
m     = 6689.174;        % kg (fuel+aircraft+passengers)
v_ind = 221     ;        % m/s
T0    = 288.15  ;        % K (temperature on day of test flight)
g0    = 9.81    ;        % m/s^2
S     = 30      ;        % m^2
rho0  = 1.225   ;        % kg/m^3
A     = 8.4387  ;        % -
tCLa  = 0.1293  ;
tCD0  = 0.0246  ;
te    = 0.5856  ;
tmeasure = [1157,1297,1426,1564,1787]; %Time at which measurements were taken in sec
Th12 = dlmread('thrust.dat','');
Tht  = sum(Th12,2);

j = 1;
for i = tmeasure
    %Index of the start of the test
    index = find(flightdata.Gps_utcSec.data >= flightdata.Gps_utcSec.data(20) + (i-100));
    index = index(1);

    %Pressure and density in metric units
    alt = flightdata.Dadc1_alt.data(index);
    pressure = 101325*(288.15/(288.15-0.0019812*(alt)))^((9.80665*0.0289644)/(8.3144598*-0.0065));
    sat = flightdata.Dadc1_sat.data(index);
    rho = pressure/(8.3144598/0.0289644*(273.15+sat));
    
    % CL-alpha curve
    v_tas = flightdata.Dadc1_tas.data(index)*0.514444444;
    W = (m-((flightdata.lh_engine_FU.data(index)+flightdata.rh_engine_FU.data(index))*0.45359))*9.81;
    CL = W/(0.5*rho*v_tas^2*S);
    alpha = flightdata.vane_AOA.data(index);
    CLalpha(1,j) = CL;
    CLalpha(2,j) = alpha;
    
    % CL-CD curve
    CD = Tht(j)/(0.5*rho*v_tas^2*S);
    CD2 = tCD0 + CL^2/(pi()*A*te);
    CLCD(1,j) = CL;
    CLCD(2,j) = CD;
    CLCD(3,j) = CL^2;
    CLCD(4,j) = CD2;
    
    % Input for thrust calculation
    h_p  = flightdata.Dadc1_alt.data(index) * 0.3;
    M = flightdata.Dadc1_mach.data(index);
    Mf1 = flightdata.lh_engine_FMF.data(index) * 0.000126;
    Mf2 = flightdata.rh_engine_FMF.data(index) * 0.000126;
    T_isa = T0 + (-0.0065*h_p);
    tat = flightdata.Dadc1_tat.data(index);  %not sure about this one
    dT_isa = (tat+273.15) - T_isa;
    Tdata(j,1) = h_p;
    Tdata(j,2) = M;
    Tdata(j,3) = dT_isa;
    Tdata(j,4) = Mf1;
    Tdata(j,5) = Mf2;
    %Tdata(j,6) = flightdata.Dadc1_cas.data(index);
    j = j + 1;
end

% Write data to matlab.dat file
dlmwrite('matlab.dat', Tdata, 'delimiter', ' ', 'newline' , 'pc');

% Plot CL-alpha curve
xa = linspace(-5,15,10);
p2 = polyfit(CLalpha(2,:),CLalpha(1,:),1);
plot(xa,polyval(p2,xa))
hold on
scatter(CLalpha(2,:),CLalpha(1,:))
CLa = p2(1,1)

% Plot CD-CL^2 curve
% xp = linspace(0,1,10);
% p1 = polyfit(CLCD(3,:),CLCD(2,:),1);
% plot(xp,polyval(p1,xp))
% hold on
% scatter(CLCD(3,:),CLCD(2,:))
% CD0 = polyval(p1,0)
% e = (p1(1,1)*pi()*A)^(-1)

% Plot CL-CD curve
% xp = linspace(0,1,10);
% p3 = polyfit(CLCD(3,:),CLCD(2,:),1);
% plot(xp,polyval(p1,xp))
% hold on
% scatter(CLCD(3,:),CLCD(2,:))
% CD0 = polyval(p1,0)
% e = (p1(1,1)*pi()*A)^(-1)



