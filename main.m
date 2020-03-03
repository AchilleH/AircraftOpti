clc; clear; close all;

%% CONSTANTS
rho = 1.18; % Air density kg/m^3
V = 0:47;       %aircraft velocity (m/s)

%% Variables
%AIRFOIL DATA
%NACA 1412 w/flap 8deg aoa default
c = 1;       %chord of aircraft wing(m)
cac = .25*c; %AC of wing chord
htailc = c; %chord of horizontal tail(m)
htailac = htailc*.25; %AC of horizontal tail
vtailc = c; %chord of vertical tail
vtailac = vtailc*.25; %AC of vertical tail
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
Tc = 1/12; %max thickness ratio, given by Niccolai

%AIRCRAFT DATA
Wb = 400; % weight of the body(minus engine) [kg]
ModNum = 26;    %Num of modules needed to filter 10^6 cf2h based on modules capable of 330cfm
Wprefilter = 1.2224;
Wcarbonfilter = 2.5;
Whyperhepafilter = 2.1044;
Wfilters = ModNum*(Wprefilter+Wcarbonfilter+Whyperhepafilter); %weight of all filters in kg
h = 1.8;       %height of aircraft (m)
C = 10; %Battery Capacity (assuming LiPo Batteries)
Pe = 20*10^3:5*10^3:75*10^3; %engine power [W]
We = 192.5:27.5:495; %engine weight, engine weight seems to inc by ~53 every 10 hp increase. 365kg for 72 hp
xdivc = .7;   %chord wise position of max thickness (m)
Sref = Sw;   %Reference Surface Area
Df = 1;      %diameter of fuselage
depthf = 0.02; %thickness of fuselage
t = Df; %approximate max horizontal thickness along the vertical (m)
Vstall = 15.5;%m/s
Vhead = 0; %headwind
%Stability Specific Variables
sspan = 10;
fuselageL = 10; % Length of the fuselage from tip USED IN NEUTRAL POINT CALC
V_max = 50; % max velocity USED IN NICCOLAI, NEEDS TO BE RECONSIDERED!!!!!!!!!!!!!

% WEIGHT DISTRIBUTION

% Order: x position, z position, weight
W_guess= Wb;
W_avionics = 7;
W_landgear = 7;

%Trial Variables to Save data
n = 500; %number of trials to run
 %Can't seem to find a way to preallocate for structs

%The Loop to calculate all our data
for i = 1:n
    %% Calculation Control
    %if you add to here, also add to save
    j = ceil(rand()*length(We)); %Chooses a random engine
    W = Wb + We(j) + Wfilters; %Total Weight
    PE = Pe(j); %engine power
    Df = rand()*5 + 1; %Df range control
    %Recalculating Swet for changing fuselage diameter
    Swet = pi*(Df/2)^2 + T*bw + T*bt;
    %% Niccolai Estimate

    [W, Ww, Wf, Wht, Wvt, We] = weight_viability(W_guess*2.205,Wfilters*2.205,W_avionics*2.205,W_landgear*2.205,We(j)*2.205,Aw,Sw*3.281^2,St*3.281^2,St*3.281^2,bt*3.281,bt*3.281,tapw,Tc,V_max*3.281,wingchord*3.281, cac*3.281, htailac*3.281, vtailac*3.281, htailc*3.281, vtailc*3.281,fuselageL*3.281,Df*3.281,depthf*3.281);
    W = W*.4536; %translating from lb to kg
    Ww = Ww*.4536; %translating from lb to kg
    Wf = Wf*.4536; %translating from lb to kg
    Wht = Wht*.4536; %translating from lb to kg
    Wvt = Wvt*.4536; %translating from lb to kg
    We = We*.4536; %translating from lb to kg
    nosearr = [0 0 10];
    avioarr = [13 0 W_avionics];
    filtarr = [16 0 Wfilters];
    fusearr = [20 0 Wf];
    htailarr = [33 0 Wht];
    engarr = [32 0 We];
    vtailarr = [33 0 Wvt];
    wingarr = [12 0 Ww];
    geararr = [wingarr(1) 0 W_landgear];
    % Order:        land gear     nose         avionics       filters       fueselage      h. tail         engine          v. tail          wing
    Xarmarray =   [ geararr(1)    nosearr(1)   avioarr(1)     filtarr(1)    fusearr(1)     htailarr(1)     engarr(1)       vtailarr(1)      wingarr(1)];     % x position from nose of masses's listed in weight array USED FOR INERTIAL AND CG CALC
    Zarmarray =   [ geararr(2)    nosearr(2)   avioarr(2)     filtarr(2)    fusearr(2)     htailarr(2)     engarr(2)       vtailarr(2)      wingarr(2)];     % z position from aircraft centerline along bottom of fuselage of masses listed in weight array USED FOR INERTIAL AND CG CALC
    weightarray = [ geararr(3)    nosearr(3)   avioarr(3)     filtarr(3)    fusearr(3)     htailarr(3)     engarr(3)       vtailarr(3)      wingarr(3)];     % masses of subsystems in aircraft USED FOR INERTIAL AND CG CALC
    downwash = 0; %downwash effect on tail
    tailac = wingarr(1)-htailarr(1)-.25*c; % Position of tail AC relative to wing LE USED IN INERTIAL AND NEUTRAL POINT CALC

    
    %% Lift
    [Sw,St,CLwa,CLta,CLw,CLt,CL,CLmax,Lw,Lt,bw,bt,Aw,At,ew,et] = lift(rho,Clwa,Clta,Clwo,Clto,Sw,St,Sref,bw,bt,Aw,At,Df,tapw,tapt,phiw,phit,V,aoarange);
    L = Lw + Lt;
    %planform surface area of plane
    S = Sw+St; %assuming only wings provide lift
    %k calc for performance input
    k = 1/(pi*Aw*S);

    %% Drag, second Swet is input for Sref in DPT
    [CDi,CDo,CD,D,Di,Do,Tr,np,Pav,Tav,Pr] = DPT(length(aoarange),PE,CL,W,Swet,Swet,Aw,Sw,At,St,ew,et,t,c,phiw,CLw,CLt,xdivc,V);
    %% Performance
    [Sto,Sl,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,Rmin,Vstall] = Perf(length(aoarange),C,rho,Tav,V,D,S,L,k,np,CL,CD,CDo,Pav,Pr,W,Vhead);
    % %% CG
    % [XCG,ZCG,Wtotal] = CG_calc(Xarmarray,Zarmarray,weightarray);
    %
    % %% Neutral Point
    % [hn] = neutral_point(0.25*c, tailac, St, Sw, CLta, CLa, downwasheffect);
    %
    % %% Static Margin
    % staticmargin = XCG/c-hn;
    %
    % %% Wing, Engine, Fuselage Inertias
    % [Ixw, Iyw, Izw, Ixzw] = wing_inertia(wingweight,type, leadingsweep,trailingsweep,roottaper,tiptaper,chord,wingstart,span,fuselageradius,planetoCG);
    % [Ixe, Iye, Ize, Ixze] = engine_inertia(engineweight,nacelleradius,XZtoengineCG,XYtoengineCG,enginelength);
    % [Ixf, Iyf, Izf, Ixzf] = fuselage_inertia(fuselageradius,noseconeweight,tailconeweight,mainfuselageweight,xyplanetocenterline,noseconelength,tailconelength,mainfuselagelength);
    % Ixarray = [Ixw,Ixe, Ixf];
    % Iyarray = [Iyw, Iye, Iyf];
    % Izarray = [Izw, Ize, Izf];
    % Ixzarray = [Ixzw, Ixze, Ixzf];
    % [Ixcg,Iycg,Izcg,Ixzcg] = inertial_calc(Ixarray,Iyarray,Izarray,Ixzarray,weightarray,Xarmarray,Zarmarray);

    %% Saving the Data, considering
    %Change the static stab. var when ryan's functions function
    XCG = 0;
    ZCG = 0;
    Wtotal = W;
    hn = 0;
    staticmargin = -1; %the passing cond
    UAV = Save(Df,We(j),Pe(j),Sw,St,CLwa,CLta,CLw,CLt,CL,CLmax,Lw,Lt,bw,bt,Aw,At,ew,et,CDi,CDo,CD,D,Di,Do,Tr,np,Pav,Tav,Pr,Sto,Sl,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,Rmin,Vstall,XCG,ZCG,Wtotal,hn,staticmargin);
    Data(i) = UAV;
end

%% Spec Verification
for i = 1:n
    Data(i).result = test(Data(i));
end
%% Result Plotting
%Preallocate arrays to hold histogram data here
%Successful Data arrays
HDf = []; %empty arrays bc. i dont want it to saturate the 0 mark

%All Data Arrays
PDf = zeros(1,n);
for i = 1:n
    if Data(i).result == true
        HDf(i) = Data(i).Df; 
    end
    PDf(i) = Data(i).Df;
end
% Add more figures following the format to plot other data
figure()
hold
histogram(HDf);
histogram(PDf);
legend('Successful','All');


