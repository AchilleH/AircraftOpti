clc; clear; close all;

%% CONSTANTS
rho = 1.18; % Air density kg/m^3
V = 0:100;       %aircraft velocity (m/s)

%% Variables
%AIRFOIL DATA
%NACA 1412 w/flap 8deg aoa default
c = 1;       %chord of aircraft wing(m)
acw = .25*c; %AC of wing chord
ct = c; %chord of horizontal tail(m)
act = ct*.25; %AC of horizontal tail
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
T = 0.1*c; %Airfoil max thickness
aoarange = -8:1:8; %range of aoa to evaluate over

%AIRCRAFT DATA
Wb = 400; % weight of the body(minus engine) [kg]
ModNum = 26;    %Num of modules needed to filter 10^6 cf2h based on modules capable of 330cfm
Wprefilter = 1.2224;
Wcarbonfilter = 2.5;
Whyperhepafilter = 2.1044;
Wfilters = ModNum*(Wprefilter+Wcarbonfilter+Whyperhepafilter); %weight of all filters in kg
h = 1.8;       %height of aircraft (m)
C = 10; %Battery Capacity (assuming LiPo Batteries)
xdivc = .7;   %chord wise position of max thickness (m)
Sref = Sw;   %Reference Surface Area
Df = 1;      %diameter of fuselage
depthf = 0.02; %thickness of fuselage
t = Df; %approximate max horizontal thickness along the vertical (m)
Vstall = 15.5;%m/s
Vhead = 0; %headwind
%Motor info
Motors = ['180M-2','200L1-2','200L2-2','225M-2','250M-2','280S-2','280M-2','315S-2'];
Pe = [22000,30000,37000,45000,55000,75000,90000,110000]; %engine power [W]
We = [165,218,230,280,365,495,565,890]; %engine weight, engine weight seems to inc by ~53 every 10 hp increase. 365kg for 72 hp
Rmot = .5* [.455,.505,.505,.560,.615,.680,.680,.845]; %Radius of motor
Lmot = [.7,.77,.77,.815,.910,.985,1.035,1.185]; %(m) Length of motor
%Stability Specific Variables
%Weights
Wnc = 0; %weight of nosecone [kg]
Wtc = 0; %weight of tailcone [kg]
W_guess= Wb;
W_avionics = 7;
W_landgear = 7;
%Lengths
Lnc = 0;%length of nosecone [m]
Ltc = 0;%length of tailcone [m]
fuselageL = 10; % Length of the fuselage from tip USED IN NEUTRAL POINT CALC [m]
V_max = max(V); % max velocity USED IN NICCOLAI, NEEDS TO BE RECONSIDERED!
Rnac = Rmot + .1; %radius of nacelle (motor housing) ASSUMING: 10cm larger diameter than motor
XZmotor = 0; %distance between motor CG and XZplane
XYmotor = 0; %distance between motor CG and XY plane
XYfuselage = 0; %distance between XY plane and fuselage centerline
Xwing = 0.4*fuselageL; %Distance between fuselage tip and tip of wing [m]
Xtail = 0.9*fuselageL; %distance between fuselage tip and tip of tail [m]
%Trial Variables to Save data
%Can't seem to find an easy/worthwhile way to preallocate for structs
n = 500; %number of trials to run

%The Loop to Run Trials
for i = 1:n
    %% Calculation Control
    %if you add to here, also add to save
    j = ceil(rand()*length(We)); %Chooses a random engine
    PE = Pe(j); %engine power
    Df = rand()*2 + 1; %Df range control
    bw = rand()*5 + 10; %wingspan randomizer
    bt = bw/2; %tail wingspan
    YZmotor = fuselageL-Lmot(j)*0.5; %distance between motor CG and YZ plane (total length minus half motor length, assuming motor CG is 1/2way)

    %Recalculating Swet for changing fuselage diameter
    Swet = pi*(Df/2)^2 + T*bw + T*bt;
    %% Niccolai Estimate
    [W, Ww, Wf, Wht, Wvt, Weng] = weight_viability(W_guess*2.205,Wfilters*2.205,W_avionics*2.205,W_landgear*2.205,We(j)*2.205,Aw,Sw*3.281^2,St*3.281^2,St*3.281^2,bt*3.281,bt*3.281,tapw,T/c,V_max*3.281,c*3.281, acw*3.281, act*3.281, vtailac*3.281, ct*3.281, vtailc*3.281,fuselageL*3.281,Df*3.281,depthf*3.281);
    %conversion from lb to kg
    W = W*.4536; %total weight
    Ww = Ww*.4536; %wing weight
    Wf = Wf*.4536; %fuselage weight
    Wht = Wht*.4536; %horizontal tail weight
    Wvt = Wvt*.4536; %vertical tail weight
    Weng = Weng*.4536; %weight of propulsion system(eng and air intake etc)
    nosearr = [0 0 Wnc]; %x position, z position, weight for nose cone
    avioarr = [2 0 W_avionics]; %x position, z position, weight for avionics
    filtarr = [2 0 Wfilters]; %x position, z position, weight for filters
    fusearr = [fuselageL/2 0 Wf]; %x position, z position, weight for fuselage
    htailarr = [9 0 Wht]; %x position, z position, weight for horizontal tail
    engarr = [8 0 Weng]; %x position, z position, weight for engine
    vtailarr = [9 0 Wvt]; %x position, z position, weight for vertical tail
    wingarr = [3 0 Ww]; %x position, z position, weight for wing
    geararr = [wingarr(1) 0 W_landgear]; %x position, z position, weight for landing gear
    tailarr = [9.5 0 Wnc]; %x position, z position, weight for landing gear
    % Order:        land gear     nose         avionics       filters       fueselage      h. tail         engine          v. tail          wing
    Xarmarray =   [ geararr(1)    nosearr(1)   avioarr(1)     filtarr(1)    fusearr(1)     htailarr(1)     engarr(1)       vtailarr(1)      wingarr(1)      tailarr(1)];     % x position from nose of masses's listed in weight array USED FOR INERTIAL AND CG CALC
    Zarmarray =   [ geararr(2)    nosearr(2)   avioarr(2)     filtarr(2)    fusearr(2)     htailarr(2)     engarr(2)       vtailarr(2)      wingarr(2)      tailarr(1)];     % z position from aircraft centerline along bottom of fuselage of masses listed in weight array USED FOR INERTIAL AND CG CALC
    weightarray = [ geararr(3)    nosearr(3)   avioarr(3)     filtarr(3)    fusearr(3)     htailarr(3)     engarr(3)       vtailarr(3)      wingarr(3)      tailarr(1)];     % masses of subsystems in aircraft USED FOR INERTIAL AND CG CALC
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
    %% CG
    [XCG,ZCG,Wtotal] = CG_calc(Xarmarray,Zarmarray,weightarray);
    
    %% Neutral Point
    [hn] = neutral_point(acw, act, St, Sw, CLta, CLwa, downwash);
    
    %% Static Margin
    staticmargin = XCG/c-hn;
    
    %% Wing, Engine, Fuselage Inertias
    %Wing inertia (type is bool, 0 = wing, 1 = tail)
    %Wing inertia seems to need at least 2 calls, maybe 4?
    %somehow they need to find their way into the final array
    %how that will affect the final calc who tf knows
    
    %Assumptions
    %assuming we need 4 calls (that the function only evaluates over 1 wing), so for wing weight val, dividing by 2
    %assuming planetoCG only cares about y dist, since were aiming for 0
    %sweep, going to assume that half the single wing span plus fuselage
    %radius is sufficient for that variable - (bw/4 + Df/4) is equivalent
    %2nd call has plane2cg as negative bc I assume we need that to balance
    %symmetric forces
    %They don't seem to account for z differences but output z values?
    %in light of previous comment I will not account for vert. tail since
    %it should be inline with fuselage centerline
    [Ixw, Iyw, Izw, Ixzw] = wing_inertia(Ww/2,false,phiw,phiw,tapw,tapw,c,Xwing,(bw-Df)/2,Df/2,bw/4 + Df/4);
    [Ixw2, Iyw2, Izw2, Ixzw2] = wing_inertia(Ww/2,false,phiw,phiw,tapw,tapw,c,Xwing,(bw-Df)/2,Df/2,-bw/4 - Df/4);
    [Ixt, Iyt, Izt, Ixzt] = wing_inertia(Wht/2,true,phit,phit,tapt,tapt,ct,Xtail,(bt-Df)/2,Df/2,-bt/4 - Df/4);
    [Ixt2, Iyt2, Izt2, Ixzt2] = wing_inertia(Wht/2,true,phit,phit,tapt,tapt,ct,Xtail,(bt-Df)/2,Df/2,-bt/4 - Df/4);
    Ixw = Ixw + Ixw2 + Ixt + Ixt2; %summation of 4 wing intertial contributions
    Iyw = Iyw + Iyw2 + Iyt + Iyt2;
    Izw = Izw + Izw2 + Izt + Izt2;
    Ixzw = Ixzw + Ixzw2 + Ixzt + Ixzt2;
    
    [Ixe, Iye, Ize, Ixze] = engine_inertia(Weng,Rnac(j),XZmotor,XYmotor,YZmotor,Lmot(j));
    [Ixf, Iyf, Izf, Ixzf] = fuselage_inertia(Df/2,Wnc,Wtc,Wf,XYfuselage,Lnc,Ltc,fuselageL);
    Ixarray = [Ixw,Ixe, Ixf];
    Iyarray = [Iyw, Iye, Iyf];
    Izarray = [Izw, Ize, Izf];
    Ixzarray = [Ixzw, Ixze, Ixzf];
    [Ixcg,Iycg,Izcg,Ixzcg] = inertial_calc(Ixarray,Iyarray,Izarray,Ixzarray,weightarray,Xarmarray,Zarmarray);

    %% Saving the Data, considering
    %Change the static stab. var when ryan's functions function
    Data(i) = Save(Df,Motors(j),Rmot(j),Lmot(j),We(j),Pe(j),Rnac(j),Sw,St,CLwa,CLta,CLw,CLt,CL,CLmax,Lw,Lt,bw,bt,Aw,At,ew,et,CDi,CDo,CD,D,Di,Do,Tr,np,Pav,Tav,Pr,Sto,Sl,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,Rmin,Vstall,XCG,ZCG,Wtotal,hn,staticmargin);
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


