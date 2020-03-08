function[CS, CP, CT, eta_p, J] = Prop(V, Pmotor, T)

% I hope this worked

% gamma = 1.4;
% Rair = 287;
rho = 1.18883546; % avg. air density (kg/m^3), range of 400-1700ft (121.92-518.16 m) above sea level
% Tair = 339.77598; % (K), range of 400-1700ft (121.92-518.16 m) above sea level
% a = sqrt(gamma*Rair*Tair); % speed of sound (m/s)
% M = V/a; % mach number

%T =         % amount of thrust applied/needed
%Pmotor =    % power supplied by motor
n = 2975/60; % propeller speed [rev/s]
eta_prop = 0.9; % propeller efficiency

Pin = Pmotor*eta_prop; % power used by propeller
% for proppitch = 10:10:70 % plots for individual pitch angles

i = 1;
for D = .9:.1:15   % different diameters, basically x-axis
    
    J(i) = V/(D*n); % advance ratio

    CT(i) = T/(rho*n^2*D^4); % thrust coefficent
    CP(i) = Pin/(rho*n^3*D^5); % power coefficent
    CS(i) = (V^5*rho/(Pin*n^2))^0.2; % speed coefficient
    eta_p(i) = CT(i)*J(i)/CP(i); % propulsive efficency
    i = i+1;
    
    % refer to NACA report 640 for plot comparison in selecting blade and
    % pitch type
    
end
    figure(1)
    plot(J,CT)
    xlabel('J')
    ylabel('CT')
    figure(2)
    plot(J,CP)
    xlabel('J')
    ylabel('CP')
    figure(3)
    plot(J,CS)
    xlabel('J')
    ylabel('Cs')
    figure(4)
    plot(J,eta_p)
    xlabel('J')
    ylabel('eta_p')
end





