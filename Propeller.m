function[outputs] = Propeller(P,D)

P = ;   %power absorbed by propeller, ft-lb/s or Hp
n = ;   %propeller speed, rps
D = ;   %propeller diameter,ft
new = ; %propulsive efficiency of propeller engine unit
J = V/(n*D);  
V = ;   %air speed fps

%%%Air Constants
    gamma = 1.4;
    Rair = 287; %(J/KgK)
    rho = 1.18883546;   %avg. air density (kg/m^3), range of 400-1700ft (121.92-518.16 m) above sea level
    Tair = 339.77598;   %(K), range of 400-1700ft (121.92-518.16 m) above sea level
    a = sqrt(gamma*Rair*Tair);  %speed of sound (m/s)

CS = V*(rho/(P*n^2))^0.2;   %something coefficient
CP = P/(rho*n^4*D^5);   %pressure coefficient
CT = T/(rho*n2*D^4);    %thrust coefficient
Mt = (V/a)*sqrt(1+(pi/J)^2); %mach number at tip of propeller blade


