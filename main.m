clc; clear; close all;

%% CONSTANTS
rho = 1.18; % Air density kg/m^3
V = [0:47];       %aircraft velocity (m/s)

%% Variables
%Variables that are assigned strictly by functions
CLw = 0; CLt = 0;
arCL = [];
L = [];
Sto = 0; Sl = 0;
C = 0;
Emax = 0; Rmax = 0;
RCmin = 0; RCmax = 0;
gamMin = 0; gamMax = 0;
Rmin = 0;

%AIRFOIL DATA
%NACA 1412 w/flap 8deg aoa default
c = 1;       %chord of aircraft wing(m)
Clw = 1.2; Clwa = 0.1; %WING 2d lift coefficient and vs alpha
Clt = Clw; Clta = Clwa; % TAIL 2d lift coef. and vs alpha
ew = 0; et = ew; %Oswald's for wing and tail
Aw = 11; At = Aw; % WING & TAIL aspect ratio
Sw = 0; St = 0;     % WING & TAIL planform area (m^2)
tapw = 0; tapt = tapw; %Taper for wing and tail
phiw = 0; phit = phiw; %Sweep from horizontal/LE
bw = 10.5; bt = bw/2; %Wingspans
aoaw = 5; aoat = aoaw; %angle of attack for wing and tail
Clwo = 0; Clto = 0; %Cl at 0 aoa(y axis offset)
T = 0.1*c; %Airfoil thickness
aoarange = -8:1:8; %range of aoa to evaluate over

%AIRCRAFT DATA
Wb = 400; % weight of the body(minus engine) [kg]
%engine weight seems to inc by ~53 every 10 hp increase. 365kg for 72 hp
h = 1.8;       %height of aircraft (m)
t = bw;       %approximate max horizontal thickness along the vertical (m)
Pe = [20*10^3:5*10^3:75*10^3]; %engine power [W]
We = [192.5:27.5:495]; %engine weight
xdivc = .7;   %chord wise position of max thickness (m)
Swet = 0;    %wetted area (m^2) CALC AGAIN BELOW
Sref = Sw;   %Reference Surface Area
Df = 1;      %diameter of fuselage
Vstall = 15.5;%m/s
Vhead = 0; %headwind
%Stability Specific Variables
sspan = 10;
fuselageL = 10; % Length of the fuselage from tip USED IN NEUTRAL POINT CALC
noseX = 0; noseL = 1; % Position of start of nose cone from tip; Length of nose cone USED IN INERTIAL AND NEUTRAL POINT CALC
wingX = 1; wingL = c; % Position of start of wing from tip; Chord of wing USED IN INERTIAL AND NEUTRAL POINT CALC
tailX = 8; tailL = c; % Position of start of tail from tip; Chord of tail USED IN INERTIAL AND NEUTRAL POINT CALC
tailconeX = 9; tailconeL = 1; % Position of start of tail cone from tip; Chord of tail cone USED IN INERTIAL AND NEUTRAL POINT CALC
tailac = wingX-tailX-.25*c; % Position of tail AC relative to wing LE USED IN INERTIAL AND NEUTRAL POINT CALC
% WEIGHT DISTRIBUTION
Xarmarray = []; % x position from nose of masses's listed in weight array USED FOR INERTIAL AND CG CALC
Zarmarray = []; % z position from aircraft centerline along bottom of fuselage of masses listed in weight array USED FOR INERTIAL AND CG CALC
weightarray = []; % masses of subsystems in aircraft USED FOR INERTIAL AND CG CALC
downwash = 0; %downwash effect on tail

%% Calculation Control
j = 1; %var to control engine choice
W = Wb + We(j);
PE = Pe(j);
%Recalculating Swet for changing fuselage diameter
Swet = pi*(Df/2)^2 + T*bw + T*bt;

%% Lift
[W,Sw,St,CLwa,CLta,CLw,CLt,CL,CLmax,Lw,Lt,bw,bt,Aw,At,ew,et] = lift(rho,Clwa,Clta,Clwo,Clto,W,Sw,St,Sref,bw,bt,Aw,At,Df,tapw,tapt,phiw,phit,V,aoarange);

%% Drag, second Swet is input for Sref in DPT
[Wfilters,CDiw,CDit,CDi,tdivc,Q,K,Cf,CDmisc,CDleak,CDprot,CDo,CD,q,D,Di,Do,Tr,np,Pav,Tav,Pr] = DPT(PE,CL,W,Swet,Swet,Aw,Sw,At,St,ew,t,c,phiw,CLw,CLt,xdivc,h,V);

%% Performance
%planform surface area of plane
S = Sw+St; %assuming only wings provide lift
%k calc
kw = 1/(pi*Aw*ew);
[Sto,Sl,C,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,Rmin] = Perf(rho,Tav,V,D,S,L,kw,np,CL,CD,CDo,Pav,Pr,W,Vhead,Vstall);

%% CG
[XCG,ZCG,Wtotal] = CG_calc(Xarmarray,Zarmarray,weightarray);

%% Neutral Point
[hn] = neutral_point(0.25*c, tailac, St, Sw, CLta, CLa, downwasheffect);

%% Static Margin
staticmargin = XCG/c-hn;

%% Spec Verification
if(Emax >= 2  && Sto < 121 && Sl < 121 && Sw <= bw*c && St <= bt*c)
fprintf('Success! \n') %in case we do generate one then add portion to display variables
fprintf('j =')
disp(j)
fprintf('tapw')
disp(tapw)
fprintf('Df')
disp(Df)
fprintf('Capacity of battery Ah')
disp(C)
fprintf('Sto =')
disp(Sto)
fprintf('Sl =')
disp(Sl)
fprintf('static margin =')
disp(staticmargin)
pause;
close all;
else
fprintf('end of test \n') %for debug
close all;
end
