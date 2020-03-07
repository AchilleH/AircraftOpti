function[CS, CP, CT] = Prop(V, Pmotor)

% I hope this worked


% gamma = 1.4;
% Rair = 287;
% rho = 1.18883546; % avg. air density (kg/m^3), range of 400-1700ft (121.92-518.16 m) above sea level
% Tair = 339.77598; % (K), range of 400-1700ft (121.92-518.16 m) above sea level
% a = sqrt(gamma*Rair*Tair); % speed of sound (m/s)
% M = V/a; % mach number

pitcharr = [10:10:70]; % range of propeller pitch angles


%Pmotor =    % power supplied by motor
n = 2975/60; % propeller speed [rev/s]
eta_prop = 0.9; % propeller efficiency

Pin = Pmotor*eta_prop; % power used by propeller


for proppitch = 10:10:70 % plots for individual pitch angles
for D = 1:15   % different diameters, basically x-axis
    
    J = V/(D*n); % advance ratio

    CT = T/(rho*n^2*D^4); % thrust coefficent
    CP = Pin/(rho*n^3*D^5); % power coefficent
    CS = V^5*sqrt(rho/(Pin*n^2)); % speed coefficient
    eta_p = CT*J/CP; % propulsive efficency
    
end
end





