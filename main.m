clc; clear; close all;

%CONSTANTS
rho = 1.18; % Air density kg/m^3
V = [0:47];       %aircraft velocity (m/s)
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

%AIRCRAFT DATA
Wb = 400; % weight of the body(minus engine) [kg] 
%engine weight seems to inc by ~53 every 10 hp increase. 365kg for 72 hp
h = 1.8;       %height of aircraft (m)
t = bw;       %approximate max horizontal thickness along the vertical (m)
Pe = [20*10^3:5*10^3:75*10^3]; %engine power [W]
We = [192.5:27.5:495];
xdivc = .7;   %chord wise position of max thickness (m)
Swet = 0;    %wetted area (m^2) CALC AGAIN BELOW
Df = 0;      %diameter of fuselage
Vstall = 15.5;%m/s
Vhead = 0; %headwind
sspan = 10;

%Actual Loop Portion
%Will iterate through several different possibilities for certain variables
for Df = [.1:.2:1.5]
    %Recalculating Swet for changing fuselage diameter
    Swet = pi*(Df/2)^2;
  %  for phiw = [0:1:15] %wing sweep iteration
  %      phit = phiw;
        for j = [1:1:length(Pe)] %engine power iteration [W]
            W = We(j) + Wb;
            PE = Pe(j);
            for tapw = [0.2:0.05:0.4] %taper ratio iteration
                tapt = tapw;
                %Lift Call        
                [W,Sw,St,CLw,CLt,CLa,CLwa,CLta,L,Vstall,bw,bt,Aw,At,ew,et] = lift(rho,Clw,Clt,Clwa,Clta,aoaw,aoat,ew,et,W,Vstall,Sw,St,bw,bt,Aw,At,Df,tapw,tapt,phiw,phit,V);
                %Getting a CL of the whole plane at stall speed for drag calc
                dex = find(V==round(Vstall));
                CL = CLa(dex);
                %Drag Call
                [Wfilters,CDiw,CDit,CDi,tdivc,Q,K,Cf,CDmisc,CDleak,CDprot,CDo,CD,q,D,Di,Do,Tr,np,Pav,Tav,Pr] = DPT(PE,CLa,W,Swet,Sw,Aw,Sw,At,St,ew,t,c,phiw,CLw,CLt,xdivc,h,V);
                %Performance Call
                %planform surface area of plane
                S = Sw+St; %assuming only wings provide lift
                %robust alt. but possibly error prone
                %S = L(1)/(.5 * arCL(1) * V(1)^2 * rho);
                %k calc
                kw = 1/(pi*Aw*ew);
                [Sto,Sl,C,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,Rmin] = Perf(rho,Tav,V,D,S,L,kw,np,CL,CD(dex),CDo(dex),Pav,Pr,W,Vhead,Vstall);
                jj = find(min(D));
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
                    pause;
                    close all;
                else
                    fprintf('end of test \n') %for debug
                    close all;
                end
            end
        end
    %end
end