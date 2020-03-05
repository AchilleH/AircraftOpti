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
vtailh = 1.5*vtailc;
htailc = c;
Clwa = 0.1; %WING 2d lift coefficient and vs alpha
Clta = Clwa; % TAIL 2d lift coef. and vs alpha
Aw = 11; At = Aw; % WING & TAIL aspect ratio
Sw = 0; St = 0;     % WING & TAIL planform area (m^2)
tapw = 1; tapt = tapw; %Taper for wing and tail
phiw = 0; phit = phiw; %Sweep from horizontal/LE
bw = 10.5; bt = bw/2; %Wingspans
Clwo = 0; Clto = 0; %Cl at 0 aoa(y axis offset)
T = 0.1*c; %Airfoil max thickness
aoarange = -8:1:8; %range of aoa to evaluate over
wingtapfac = 1; % "wing taper factor", or what fraction of the root thickness is the tip thickness
tailtapfac = 1; % "tail taper factor", or what fraction of the root thickness is the tip thickness
wingroot_T = T; % thickness of wing airfoil root
wingtip_T = wingtapfac*wingroot_T; % thickness of wing airfoil tip
tailroot_T = T; % thickness of tail airfoil root
tailtip_T = tailtapfac*tailroot_T; % thickness of tail airfoil tip

%AIRCRAFT DATA
Wb = 400; % weight of the body(minus engine) [kg]
ModNum = 26;    %Num of modules needed to filter 10^6 cf2h based on modules capable of 330cfm
Wprefilter = 1.2224;
Wcarbonfilter = 2.5;
Whyperhepafilter = 2.1044;
Wfilters = ModNum*(Wprefilter+Wcarbonfilter+Whyperhepafilter); %weight of all filters in kg
doubledup = 0; % tracker to see how many times filter doubling up happened
h = 1.8;       %height of aircraft (m)
C = 10; %Battery Capacity (assuming LiPo Batteries)
xdivc = .7;   %chord wise position of max thickness (m)
Sref = Sw;   %Reference Surface Area
Df = 1;      %diameter of fuselage
depthf = 0.02; %thickness of fuselage
t = bw; %approximate max horizontal thickness along the vertical (m)
Vstall = 15.5;%m/s
Vhead = 0; %headwind
%Motor info
Motors = ["160M1-2","160M2-2","160L-2","180M-2","200L1-2","200L2-2","225M-2","250M-2","280S-2","280M-2","315S-2"];
Pe = [11000,15000,18500,22000,30000,37000,45000,55000,75000,90000,110000]; %engine power [W]
We = [110,120,135,165,218,230,280,365,495,565,890]; %engine weight, engine weight seems to inc by ~53 every 10 hp increase. 365kg for 72 hp
Rmot = .5* [.420,.420,.420,.455,.505,.505,.560,.615,.680,.680,.845]; %Radius of motor
Lmot = [.615,.615,.670,.7,.77,.77,.815,.910,.985,1.035,1.185]; %(m) Length of motor

%Stability Specific Variables
%Weights
Wnc = 0; %weight of nosecone [kg]
Wtc = 0; %weight of tailcone [kg]
W_guess= Wb;
W_landgear = 7;
Wbat = 128; %weight of the batteries [kg]
W_avionics = 20 + Wbat;
%Lengths
Lnc = 0;%length of nosecone [m]
Ltc = 0;%length of tailcone [m]
fuselageL = 15; % Length of the fuselage from tip USED IN NEUTRAL POINT CALC [m]
V_max = max(V); % max velocity USED IN NICCOLAI, NEEDS TO BE RECONSIDERED!
Rnac = Rmot + .1; %radius of nacelle (motor housing) ASSUMING: 10cm larger diameter than motor
XZmotor = 0; %distance between motor CG and XZplane
XYmotor = 0; %distance between motor CG and XY plane
XYfuselage = 0; %distance between XY plane and fuselage centerline
Xwing = 0.4*fuselageL; %Distance between fuselage tip and tip of wing [m]
Xtail = 0.9*fuselageL; %distance between fuselage tip and tip of tail [m]
Xeng = 0; %distance of the eng CG wrt fuselage tip
%Trial Variables to Save data
%Can't seem to find an easy/worthwhile way to preallocate for structs
n = 100000; %number of trials to run

%The Loop to Run Trials
for i = 1:n
    %% Calculation Control
    %if you add to here, also add to save
    j = ceil(rand()*length(We)); %Chooses a random engine
    PE = Pe(j); %engine power
    Winv = ((PE-0.75)/(450-0.75))*297 + 3; %inverter weight estimate based on alibaba specs
    Wbat = Wbat + Winv; %adjusted batt weight

    Df = rand()*2 + 1; %Df range control
    bw = rand()*5 + 10; %wingspan randomizer
    t = bw; %approximate max horizontal thickness along the vertical (m)
    bt = rand()*(bw-7) + 5; %tail wingspan
    c = rand()*1 + 1; %wing chord randomizer
    ct = rand()*(c-0.7) + 0.5; %tail chord randomizer

    Xeng = rand()*(fuselageL - 0.6*fuselageL) + 0.6*fuselageL; %randomized the location of the motor cg between 0.6 and fuselageL
    YZmotor = fuselageL-Lmot(j)*0.5; %distance between motor CG and YZ plane (total length minus half motor length, assuming motor CG is 1/2way)

    %% Dealing with Planform area, AR, & B; Also adjusts Sref if input is 0
    if(Sw == 0 && Aw ~= 0 && bw ~= 0)
       Sw = bw^2/Aw;
    end
    if(Sw ~= 0 && Aw == 0 && bw ~= 0)
       Aw = bw^2/Sw;
    end
    if(Sw ~= 0 && Aw ~= 0 && bw == 0)
       bw = sqrt(Sw*Aw);
    end

    if(St == 0 && At ~= 0 && bt ~= 0)
       St = bt^2/At;
    end
    if(St ~= 0 && At == 0 && bt ~= 0)
       At = bt^2/St;
    end
    if(St ~= 0 && At ~= 0 && bt == 0)
       bt = sqrt(St*At);
    end
    if(Sref == 0)
        Sref = Sw;
    end
    %Recalculating Swet for changing fuselage diameter
    Swet = pi*(Df/2)^2 + T*bw + T*bt + T*vtailh;
    %% Niccolai Estimate
    [W, Ww, Wf, Wht, Wvt, Weng] = weight_viability(W_guess*2.205,Wfilters*2.205,W_avionics*2.205,W_landgear*2.205,We(j)*2.205,Aw,Sw*3.281^2,St*3.281^2,St*3.281^2,bt*3.281,bt*3.281,tapw,T/c,V_max*3.281,c*3.281, acw*3.281, act*3.281, vtailac*3.281, ct*3.281, vtailc*3.281,fuselageL*3.281,Df*3.281,depthf*3.281);
    %conversion from lb to kg
    W = W*.4536; %total weight
    Ww = Ww*.4536; %wing weight
    Wf = Wf*.4536; %fuselage weight
    Wht = Wht*.4536; %horizontal tail weight
    Wvt = Wvt*.4536; %vertical tail weight
    Weng = Weng*.4536; %weight of propulsion system (eng and air intake etc)

    %nose cone: length,  x position (CG),  z position (CG),  weight
    nosearr =  [0.5,     .75*0.5,             Df/2,             Wnc];
    %tail cone: length,  x position (CG),               z position (CG),  weight
    tailarr =  [0.5,     fuselageL-nosearr(1)+.25,      Df/2,             Wtc];
    %avioncis: length, x position (CG), z position (CG), weight
    avioarr = [0,      nosearr(2),      Df/2,            W_avionics];
    %filters:  length,                          x position (CG),  z position (CG), weight
    filtarr = [fuselageL-nosearr(1)-tailarr(1), nosearr(1)+(fuselageL-nosearr(1)-tailarr(1)/2),      Df/2,            Wfilters];
    %fuselage: length,                          x position (CG),  z position (CG), weight
    fusearr = [fuselageL-nosearr(1)-tailarr(1), fuselageL/2,      Df/2,            Wf];
    %h. tail:   length, x position (CG), z position (CG), weight
    htailarr = [htailc, tailarr(2),      Df/2,            Wht];
    %engine:  length,  x position (CG),  z position (CG), weight
    engarr = [Lmot(j), Xeng,      Df/2,            Weng];
    %v. tail:   length, x position (CG),  z position (CG), weight
    vtailarr = [vtailc, htailarr(2),      Df/2,            Wvt];
    %wing:     length, x position (CG),                    z position (CG), weight
    wingarr = [c,      .4*fuselageL+(.4*fuselageL)/2,      Df/2,            Ww];
    %landing gear: length, x position (CG), z position (CG), weight
    geararr =     [0,      wingarr(2),      Df/2,            W_landgear];
    %battery: length, x position (CG), z position (CG), weight
    battarr =      [0.5,      nosearr(1),   Df/2,            Wbat];
    % doubling up filters if leftover diameter with one filter is more than
    % 5/8ths the diameter (somewhat arbitrary, can be changed to another
    % preference)
    if Df-.33 >= Df*5/8
        filtarr(2) = nosearr(1)+(fuselageL-nosearr(1)-tailarr(1)/4); % halves effect of filters on CG location
        doubledup = doubledup+1; % counts one doubling up
    end
    %               landing gear  nose cone    avionics       filters       fuselage       h. tail         engine          v. tail          wing            tail cone     battery
    Xarmarray =   [ geararr(2)    nosearr(2)   avioarr(2)     filtarr(2)    fusearr(2)     htailarr(2)     engarr(2)       vtailarr(2)      wingarr(2)      tailarr(2)    battarr(2)];     % x position from nose of masses's listed in weight array USED FOR INERTIAL AND CG CALC
    Zarmarray =   [ geararr(3)    nosearr(3)   avioarr(3)     filtarr(3)    fusearr(3)     htailarr(3)     engarr(3)       vtailarr(3)      wingarr(3)      tailarr(3)    battarr(3)];     % z position from aircraft centerline along bottom of fuselage of masses listed in weight array USED FOR INERTIAL AND CG CALC
    weightarray = [ geararr(4)    nosearr(4)   avioarr(4)     filtarr(4)    fusearr(4)     htailarr(4)     engarr(4)       vtailarr(4)      wingarr(4)      tailarr(4)    battarr(4)];     % masses of subsystems in aircraft USED FOR INERTIAL AND CG CALC
    downwash = 0; %downwash effect on tail
    tailac = wingarr(2)-htailarr(2)-.25*c; % Position of tail AC relative to wing LE USED IN INERTIAL AND NEUTRAL POINT CALC

    %% Lift
    [Sw,St,CLwa,CLta,CLw,CLt,CL,CLmax,Lw,Lt,bw,bt,Aw,At,ew,et] = lift(rho,Clwa,Clta,Clwo,Clto,Sw,St,Sref,bw,bt,Aw,At,Df,tapw,tapt,phiw,phit,V,aoarange);
    L = Lw + Lt;
    %planform surface area of plane
    S = Sw+St; %assuming only wings provide lift
    %k calc for performance input
    k = 1/(pi*Aw*S);

    %% Drag, second Swet is input for Sref in DPT
    [CDi,CDo,CD,D,Ds,Di,Dis,Do,Tr,np,Pav,Tav,Pr] = DPT(length(aoarange),PE,CL,W,Swet,Swet,Aw,Sw,At,St,ew,et,t,c,phiw,CLw,CLt,xdivc,V);
    %% Performance
    [Sto,Sl,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,Rmin,Vstall] = Perf(length(aoarange),C,rho,Tav,V,D,S,L,k,np,CL,CD,CDo,Pav,Pr,W,Vhead);
    %% CG
    [XCG,ZCG,Wtotal] = CG_calc(Xarmarray,Zarmarray,weightarray);

    %% Neutral Point
    [hn] = neutral_point(acw, act + (Xtail-Xwing), St, Sw, CLta, CLwa, downwash);

    %% Static Margin
    staticmargin = (XCG-Xwing)/c - hn;

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

    % Wing Inertial Contributions
    % left wing
    [Ixw, Iyw, Izw, Ixzw] =     wing_inertia(Ww/2,false,phiw,phiw,wingroot_T,wingtip_T,c,Xwing,(bw-Df)/2,Df/2,bw/4 + Df/4);
    % right wing
    [Ixw2, Iyw2, Izw2, Ixzw2] = wing_inertia(Ww/2,false,phiw,phiw,wingroot_T,wingtip_T,c,Xwing,(bw-Df)/2,Df/2,-bw/4 - Df/4);
    % left tail
    [Ixt, Iyt, Izt, Ixzt] =     wing_inertia(Wht/2,true,phit,phit,tailroot_T,tailtip_T,ct,Xtail,(bt-Df)/2,Df/2,-bt/4 - Df/4);
    % right tail
    [Ixt2, Iyt2, Izt2, Ixzt2] = wing_inertia(Wht/2,true,phit,phit,tailroot_T,tailtip_T,ct,Xtail,(bt-Df)/2,Df/2,-bt/4 - Df/4);
    % v. tail
    [Ixvt, Iyvt, Izvt, Ixzvt] = wing_inertia(Wht/2,true,phit,phit,tailroot_T,tailtip_T,ct,Xtail,(bt-Df)/2,Df/2,-bt/4 - Df/4);
    % engine
    [Ixe, Iye, Ize, Ixze] = engine_inertia(Weng,Rnac(j),XZmotor,XYmotor,YZmotor,Lmot(j));
    % fuselage (includes nose and tail cones)
    [Ixf, Iyf, Izf, Ixzf] = fuselage_inertia(Df/2,Wnc,Wtc,Wf,XYfuselage,Lnc,Ltc,fuselageL);
    % filters
    [Ixfilt, Iyfilt, Izfilt, Ixzfilt] = box_inertia(filtarr(4),filtarr(1),Df/2,filtarr(3),filtarr(2));
    % avionics
    [Ixavi, Iyavi, Izavi, Ixzavi] = point_inertia(avioarr(4),avioarr(2),0,avioarr(3));
    % landing gear
    [Ixlgear, Iylgear, Izlgear, Ixzlgear] = point_inertia(avioarr(4),avioarr(2),0,avioarr(3));
    % batteries
    [Ixbatt, Iybatt, Izbatt, Ixzbatt] = box_inertia(battarr(4),battarr(1),Df,battarr(3),battarr(2));

    Ixtot = Ixw+Ixw2+Ixt+Ixt2+Ixvt+Ixe+Ixf+Ixfilt+Ixavi+Ixlgear+Ixbatt;
    Iytot = Iyw+Iyw2+Iyt+Iyt2+Iyvt+Iye+Iyf+Iyfilt+Iyavi+Iylgear+Iybatt;
    Iztot = Izw+Izw2+Izt+Izt2+Izvt+Ize+Izf+Izfilt+Izavi+Izlgear+Izbatt;
    Ixztot = Ixzw+Ixzw2+Ixzt+Ixzt2+Ixzvt+Ixze+Ixzf+Ixzfilt+Ixzavi+Ixzlgear+Ixzbatt;
    [Ixcg,Iycg,Izcg,Ixzcg] = inertial_calc(Ixtot,Iytot,Iztot,Ixztot,weightarray,Xarmarray,Zarmarray);

    %% Saving the Data, considering
    %Change the static stab. var when ryan's functions function
    Data(i) = Save(Df,Motors(j),Rmot(j),Lmot(j),We(j),Pe(j),Rnac(j),Sw,St,CLwa,CLta,CLw,CLt,CL,CLmax,Lw,Lt,bw,bt,Aw,At,ew,et,CDi,CDo,CD,D,Ds,Di,Dis,Do,Tr,np,Pav,Tav,Pr,Sto,Sl,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,Rmin,Vstall,XCG,ZCG,Wtotal,hn,Xeng,staticmargin);
end

%% Spec Verification
for i = 1:n
    Data(i).result = test(Data(i));
end
%% Result Plotting
%Preallocate arrays to hold histogram data here
%Successful Data arrays
HDf = []; %empty arrays bc. i dont want it to saturate the 0 mark
HWe = [];
HPe = [];
Hbw = [];
HXe = [];
%All Data Arrays
PDf = zeros(1,n);
PWe = zeros(1,n);
PPe = zeros(1,n);
Pbw = zeros(1,n);
PXe = zeros(1,n);

i2 = 1; %for successful trials
for i = 1:n
    if Data(i).result == true
        HDf(i2) = Data(i).Df;
        HWe(i2) = Data(i).We;
        HPe(i2) = Data(i).Pe;
        Hbw(i2) = Data(i).bw;
        HXe(i2) = Data(i).Xeng;
        i2 = i2+1;
    end
    PDf(i) = Data(i).Df;
    PWe(i) = Data(i).We;
    PPe(i) = Data(i).Pe;
    Pbw(i) = Data(i).bw;
    PXe(i) = Data(i).Xeng;
end
% Add more figures following the format to plot other data
dataAnalysis(HDf,PDf,'Df');
dataAnalysis(HWe,PWe,'We');
dataAnalysis(HPe,PPe,'Pe');
dataAnalysis(Hbw,Pbw,'bw');
dataAnalysis(HXe,PXe,'Xeng');
