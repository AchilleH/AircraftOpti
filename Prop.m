function[CS, CP, CT] = Prop(V,Pe)

gamma = 1.4;
Rair = 287;     %(J/KgK)
rho = 1.18883546;     %avg. air density (kg/m^3), range of 400-1700ft (121.92-518.16 m) above sea level
Tair = 339.77598;        %(K), range of 400-1700ft (121.92-518.16 m) above sea level
a = sqrt(gamma*Rair*Tair);      %speed of sound (m/s)

%Pe= ;   %power absorbed by prop, W 
n = 2975/60;   %prop speed rps
new = 0.9; %propulsive efficiency of prop unit

for D = 1:15   %prop diameter, m
CS = V.*(rho/.Pe/n^2)^0.2; %some coefficient
CP = Pe./rho/n^3/D^5; %pressure coefficient
%CT = T/rho/n^2/D^4; 
end
end





