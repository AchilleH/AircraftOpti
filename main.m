clc; clear; close all;

%% CONSTANTS
rho = 1.18; % Air density kg/m^3
V = 0:47;       %aircraft velocity (m/s)

%% Variables
%AIRFOIL DATA
%NACA 1412 w/flap 8deg aoa default
c = 1;       %chord of aircraft wing(m)
Clwa = 0.1; %WING 2d lift coefficient and vs alpha
Clta = Clwa; % TAIL 2d lift coef. and vs alpha
Aw = 11; At = Aw; % WING & TAIL aspect ratio
Sw = 0; St = 0;     % WING & TAIL planform area (m^2)
tapw = 0; tapt = tapw; %Taper for wing and tail
phiw = 0; phit = phiw; %Sweep from horizontal/LE
bw = 10.5; bt = bw/2; %Wingspans
Clwo = 0; Clto = 0; %Cl at 0 aoa(y axis offset)
T = 0.1*c; %Airfoil thickness
aoarange = -8:1:8; %range of aoa to evaluate over

%AIRCRAFT DATA
Wb = 400; % weight of the body(minus engine) [kg]
ModNum = 26;    %Num of modules needed to filter 10^6 cf2h based on modules capable of 330cfm
Wprefilter = 1.2224;
Wcarbonfilter = 2.5;
Whyperhepafilter = 2.1044;
Wfilters = ModNum*(Wprefilter+Wcarbonfilter+Whyperhepafilter); %weight of all filters in kg
%engine weight seems to inc by ~53 every 10 hp increase. 365kg for 72 hp
h = 1.8;       %height of aircraft (m)
Pe = 20*10^3:5*10^3:75*10^3; %engine power [W]
We = 192.5:27.5:495; %engine weight
xdivc = .7;   %chord wise position of max thickness (m)
Sref = Sw;   %Reference Surface Area
Df = 1;      %diameter of fuselage
t = Df; %approximate max horizontal thickness along the vertical (m)
Vstall = 15.5;%m/s
Vhead = 0; %headwind
%Stability Specific Variables
sspan = 10;
fuselageL = 33; % Length of the fuselage from tip USED IN NEUTRAL POINT CALC

% WEIGHT DISTRIBUTION
% Order: x position, z position, weight
nosearr = [0 0 10];
avioarr = [13 0 10];
filtarr = [16 0 Wfilters];
fusearr = [20 0 98.85];
htailarr = [33 0 5.416];
engarr = [32 0 We];
vtailarr = [33 0 3.2625];
wingarr = [12 0 125.45];
% Order:                     avionics                     fueselage                      engine                           wing
Xarmarray =   [ nosearr(1)   avioarr(1)     filtarr(1)    fusearr(1)     htailarr(1)     engarr(1)       vtailarr(1)      wingarr(1)];      % x position from nose of masses's listed in weight array USED FOR INERTIAL AND CG CALC
Zarmarray =   [ nosearr(2)   avioarr(2)     filtarr(2)    fusearr(2)     htailarr(2)     engarr(2)       vtailarr(2)      wingarr(2)];          % z position from aircraft centerline along bottom of fuselage of masses listed in weight array USED FOR INERTIAL AND CG CALC
weightarray = [ nosearr(3)   avioarr(3)     filtarr(3)    fusearr(3)     htailarr(3)     engarr(3)       vtailarr(3)      wingarr(3)];     % masses of subsystems in aircraft USED FOR INERTIAL AND CG CALC
% Order:        nose                        filters                      h. tail                         vertical tail
downwash = 0; %downwash effect on tail
tailac = wingarr(1)-htailarr(1)-.25*c; % Position of tail AC relative to wing LE USED IN INERTIAL AND NEUTRAL POINT CALC

%% Calculation Control
j = 1; %var to control engine choice
W = Wb + We(j) + Wfilters; %Total Weight
PE = Pe(j); %engine power
C = 0; %Battery Capacity
%Recalculating Swet for changing fuselage diameter
Swet = pi*(Df/2)^2 + T*bw + T*bt;
%planform surface area of plane
S = Sw+St; %assuming only wings provide lift
%k calc for performance
k = 1/(pi*Aw*S);

%% Lift
[Sw,St,CLwa,CLta,CLw,CLt,CL,CLmax,Lw,Lt,bw,bt,Aw,At,ew,et] = lift(rho,Clwa,Clta,Clwo,Clto,Sw,St,Sref,bw,bt,Aw,At,Df,tapw,tapt,phiw,phit,V,aoarange);

%% Drag, second Swet is input for Sref in DPT
[CDi,Q,CDo,CD,D,Di,Do,Tr,np,Pav,Tav,Pr] = DPT(PE,CL,W,Swet,Swet,Aw,Sw,At,St,ew,et,t,c,phiw,CLw,CLt,xdivc,V);

%% Performance
[Sto,Sl,C,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,Rmin,Vstall] = Perf(length(aoarange),C,rho,Tav,V,D,S,L,k,np,CL,CD,CDo,Pav,Pr,W,Vhead,Vstall);

%% CG
[XCG,ZCG,Wtotal] = CG_calc(Xarmarray,Zarmarray,weightarray);

%% Neutral Point
[hn] = neutral_point(0.25*c, tailac, St, Sw, CLta, CLa, downwasheffect);

%% Static Margin
staticmargin = XCG/c-hn;

%% Wing, Engine, Fuselage Inertias
[Ixw, Iyw, Izw, Ixzw] = wing_inertia(wingweight,type, leadingsweep,trailingsweep,roottaper,tiptaper,chord,wingstart,span,fuselageradius,planetoCG);
[Ixe, Iye, Ize, Ixze] = engine_inertia(engineweight,nacelleradius,XZtoengineCG,XYtoengineCG,enginelength);
[Ixf, Iyf, Izf, Ixzf] = fuselage_inertia(fuselageradius,noseconeweight,tailconeweight,mainfuselageweight,xyplanetocenterline,noseconelength,tailconelength,mainfuselagelength);
Ixarray = [Ixw,Ixe, Ixf];
Iyarray = [Iyw, Iye, Iyf];
Izarray = [Izw, Ize, Izf];
Ixzarray = [Ixzw, Ixze, Ixzf];
[Ixcg,Iycg,Izcg,Ixzcg] = inertial_calc(Ixarray,Iyarray,Izarray,Ixzarray,weightarray,Xarmarray,Zarmarray);

%% Spec Verification
%TODO: Have a Verification Test for all concatenated data, then output histograms

%Old Spec Verification
% if(Emax >= 2  && Sto < 121 && Sl < 121 && Sw <= bw*c && St <= bt*c)
% fprintf('Success! \n') %in case we do generate one then add portion to display variables
% fprintf('j =')
% disp(j)
% fprintf('tapw')
% disp(tapw)
% fprintf('Df')
% disp(Df)
% fprintf('Capacity of battery Ah')
% disp(C)
% fprintf('Sto =')
% disp(Sto)
% fprintf('Sl =')
% disp(Sl)
% fprintf('static margin =')
% disp(staticmargin)
% pause;
% close all;
% else
% fprintf('end of test \n') %for debug
% close all;
% end
