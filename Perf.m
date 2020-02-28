function [Sto,Sl,C,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,TRmin] = Perf(rho,T,V,D,S,L,k,nu,Clg,Cdg,Cdo,Pav,Preq,W,Vhead,Vstall)
%general Performance notes
%Tr = W/(L/D)
%Pr = Tr*V = D * V
%Tr=D for level flight
%Battery Related Notes
%Rt[hr] is the time that the battery of choice was evaluated over(i.e 1Ah = 1
%amp for 1 hour so Rt is 1 hour.
%n is a discharge parameter specific to each battery and changes throughout
%it's lifespan
%C is battery capacity in ampere hours

%Below are Variable values used/assumptions made
g = 9.81; %metric g
Vlof = 1.1*Vstall;
Vtd = 0.7*Vstall;
Va = 1.3*Vstall;
V2 = 1.2*Vstall;
in1 = find(V==round(1.1*Vstall)); %gets the index of the closest V for avg takeoff V
in2 = find(V==round(1.25*Vstall)); %same but landing
in3 = find(V==round(0.7*Vtd)); %for touchdown acceleration assumption
in4 = find(V==round(Vstall));
T1 = T(in1); D1 = T(in1);
T2 = T(in2); D2 = D(in2);
T3 = T(in3); D3 = D(in3); L3 = L(in3);

mu = 0.35; %hard braking friction
q = .5.*rho.*V.^2.*S;
Rt = 1; %For Ah batteries rated over 1 hour.
n = 1.3; %Typical Lithium Polymer battery value
C = 0; %this program will generate batt capacity that satisfies endurance spec


%Climbing Performance
%maybe take array inputs to get a spectrum and with it our min and max
%Cl(maxRC) = sqrt(3*Cdo*pi*AR*e); %For Propeller Aircraft
RC = (Pav(in4:end) - Preq(in4:end))./W;
RCmin = min(RC);
RCmax = max(RC);
gam = (Pav(in4:end) - Preq(in4:end))./(W.*V(in4:end));
gamMin = min(gam);
gamMax = max(gam);

%To/L
%Sg is Ground roll for lift off, denominator is evaluated at V = 0.7Vlof
%Sa is air distance: use average T-D, V2 = 1.2 Vstall
%Vlof is 1.1*Vmu, Vmu is the min speed aircraft could TO w/o stiking tail
%on ground
%TO and landing calc not working, TO goest to inf, landing hovers around 0
Sg = ((Vlof - Vhead)^2 / (2*g*(T1/W - mu - (Cdg - mu*Clg)*q(in1)/(W/S))));
Sa = W / (T1 - 0.01*D1) * floor((V2^2 - Vlof^2)/(2*g) + 15.24); %Relies on faa screen height converted to meters
Sto = abs(Sg) + abs(Sa);

%landing, a is avg acceleration at 0.7Vtd
a = T3-D3-mu*(W-L3); %use
Slg = (Vtd - Vhead)^2/(2*a);
Sa = W / (T2-D2) * floor((Va^2 - Vtd^2)/(2*g) + 15.24); %15.24 is faa screen height in meters
Sl = abs(Slg) + abs(Sa);

%Endurance and Range
%metric units, km for distance
%t = Rt/i^n *(c/Rt)^n; %battery endurance
Emax = 0;
while(Emax < 2.5)
Emax = Rt^(1-n)*((nu(in4)*Vstall*C)/((2/sqrt(rho*S))*Cdo^.25*(2*W*sqrt(k/3)^1.5)))^n;
C = C+.1;
end
Rmax = Rt^(1-n)*(nu(in4)*Vstall*C/((2/sqrt(rho*S))*Cdo^.25*(2*W*sqrt(k)^1.5)))^n * sqrt(2*W/(rho*S) * sqrt(k/Cdo)) * 3.6;

%Turning
%nlc = .5*rho*V^2 * (Clmax/(W/S)); %lift constrained load factor
n = L./W; %fyi
R = 2.*q./(g.*rho.*sqrt(n.^2 - 1));
TRmin = min(R);
end