%function[CS, CP, CT, eta_p, J] = Prop(V, Pmotor, T)
clear all; clc; close all

% I hope this worked
% gamma = 1.4;
% Rair = 287;
rho = 1.18883546; % avg. air density (kg/m^3), range of 400-1700ft (121.92-518.16 m) above sea level
% Tair = 339.77598; % (K), range of 400-1700ft (121.92-518.16 m) above sea level
% a = sqrt(gamma*Rair*Tair); % speed of sound (m/s)
% M = V/a; % mach number

T = 268.25;        % amount of thrust applied/needed, Nm
Pmotor = 37*10^3;    % power supplied by motor, W
n = 2975/60; % propeller speed [rev/s]
%Q = 0:12:120; %torque of engine Nm
% for proppitch = 10:10:70 % plots for individual pitch angles

i = 1;
figure
 D = 2;
    r = D/2; A = pi*r^2;
for V = 0:100
    eta_p(i) = (0.5+sqrt(0.25+T/2/rho/A/(V^2)))^(-1);
    J(i) = V/n/D;
    beta(i) = atan(V/n/r)*(180/pi);  %pitch angle in degrees
    CT(i) = T/(rho*n^2*D^4); % thrust coefficent
    Pin = Pmotor*eta_p(i); % power used by propeller
    %Pin = 2*pi*n*Q(i); % power used by propeller
    CP(i) = Pin/(rho*n^3*D^5); % power coefficent
    %eta_p(i) = CT(i)*J(i)/CP(i); % propulsive efficency
    i = i+1;
end
hold on 
yyaxis right
plot(J,eta_p)
ylabel('eta_p')
yyaxis left
plot(J, CT)
plot(J, CP)
xlabel('J')
ylabel('CT, CP')
legend



%% stuff ryan did
% i = 1;
% V = 30; %speed, m/s
% for D = .9:.1:15   % different diameters, basically x-axis
%     
%     J(i) = V/(D*n); % advance ratio
%     CT(i) = T/(rho*n^2*D^4); % thrust coefficent
%     CP(i) = Pin/(rho*n^3*D^5); % power coefficent
%     CS(i) = (V^5*rho/(Pin*n^2))^0.2; % speed coefficient
%     eta_p(i) = CT(i)*J(i)/CP(i); % propulsive efficency
%     i = i+1;
%     hold on 
%     plot(J,eta_p(i))
%     xlabel('J')
%     ylabel('eta_p')
%     
%     % refer to NACA report 640 for plot comparison in selecting blade and
%     % pitch type
%     
% end
%     figure(1)
%     plot(J,CT)
%     xlabel('J')
%     ylabel('CT')
%     figure(2)
%     plot(J,CP)
%     xlabel('J')
%     ylabel('CP')
%     figure(3)
%     plot(J,CS)
%     xlabel('J')
%     ylabel('Cs')
%     figure(4)
%     plot(J,eta_p)
%     xlabel('J')
%     ylabel('eta_p')
%end






