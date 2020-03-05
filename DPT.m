function [CDi,CDo,CD,D,Di,Do,Tr,np,Pav,Tav,Pr] = DPT(n,PE,CL,W,Swet,Sref,Aw,Sw,At,St,ew,et,t,c,Y,CLw,CLt,xdivc,V)
%Use length of aircraft for c

%%%154A Drag Block
%%%Air Constants
    gamma = 1.4;
    Rair = 287;     %(J/KgK)
    rho = 1.18883546;     %avg. air density (kg/m^3), range of 400-1700ft (121.92-518.16 m) above sea level
    Tair = 339.77598;        %(K), range of 400-1700ft (121.92-518.16 m) above sea level
%%%Aircraft Constants
%     W = ;
%    g = 9.81;       %m/s^2
%     Swet =;        %wetted area (m^2)
%     Sref =;        %reference area. ie wing area (m^2)
%     Aw = ;      % wing aspect ratio
%     Sw = Sref;      % wing area (m^2)
%     At = ;      % tail aspect ratio
%     St = ;      % tail area (m^2)
%     e = ;       % oswald efficiency factor
%     h = ;       %height of aircraft (m)
%     t = ;       %approximate max horizontal thickness along the vertical (m)
%     c = ;       %chord of aircraft (m)
%     Y = ;       %max thickness line sweep angl (degrees)
%    V = [0:47];       %aircraft velocity (m/s)
    a = sqrt(gamma*Rair*Tair);      %speed of sound (m/s)
    M = V./a;        %mach number
    mew = 1.77907876*10^(-6);     %dynamic viscosity of fluid (N s/m^2)
    v = mew/rho;     %kinematic viscosity (m^2/s)
    Re = V.*c/v;       %reynold's number
    q = 0.5*rho*Swet*(V.^2);
    m = length(V);

%%
CDo = zeros(n,m);
CDi = zeros(n,m);
CD = zeros(n,m);
Tr = zeros(n,m);
D = zeros(n,m);
Di = zeros(n,m);
Do = zeros(n,m);
for nn = 1:n
%%%Calculating Induced Drag Coefficient
    CDiw = (CLw(nn,:).^2)/pi/Aw/ew;      % CLw = coef. of lift of wing
    CDit = (CLt(nn,:).^2)/pi/At/et;      % CLt = coef. of lift of tail
    %%%Total Induced Drag Coefficient, CDi
    CDi(nn,:) = CDiw + (St/Sw)*CDit;
%%%Calculating Parasitic Drag Coefficient
    %Form Factor, K
    tdivc = t/c;        %thickness to chord ratio
    K = (1 + (0.6/xdivc)*tdivc + 100*(tdivc^4))*1.34*(M.^0.18)*((cos(Y))^0.28);
    %Interference Factor, Q
    Q = 1;      %fuselages and wings assumed negligible interference
    %Flat Plate Skin Friction, Cf
    Cf = 0.455./((log(Re).^2.58).*(1 + 0.144*M.^2).^0.65);
    %Miscellaneous Drag, CDmisc (due to flaps or unretracted landing gear)
    CDmisc = 0;
    %Leakage Drag, CDleak (due to 'inhalation' and 'exhalation' of filtered air)
    intakearea = 0.3302^2;  %intake area in m^2
    A = intakearea;
    air = 330;  %air delivered cfm
    vel = air/60*(.3048^3)/A;  %inlet velocity in m/s
    mdot = rho*A*vel;
    Rdrag = mdot*vel;
%     CDleak = Rdrag/q;
%%old CDleak calculation (resulted in 10^6 drag)
%     fricfrac1 = 0.28;   %prefilter friction factor
%     fricfrac2 = 0.6;    %carbon filter friction factor
%     fricfrac3 = fricfrac1;  %HyperHepa filter friction factor
%     Pair = g*10^3;    %pressure of still air in Force/m^2
%     intakearea = 0.3302^2;  %intake area in m^2
%     N = Pair*intakearea;    %normal force of still air
%     CDleakpercubicmeter = (fricfrac1+fricfrac2+fricfrac3)*N;
%     FilterSA1 = 2.787;  %PreFilter surface area in m^2
%     FilterSA2 = 3.894;  %Carbon filter surface area in m^2
%     FilterSA3 = 4.924;  %HyperHepa filter surface area in m^2
%     SingleModuleSA = FilterSA1+FilterSA2+FilterSA3; %total surface area of one module (3 air filters) in m^2
    %Protuberance Drag, CDprot (due to things that stick out of the body,
    %like the booms for example)
    CDprot = 0;
    %%%Total Parasitic Drag Coefficient, CDo
    CDo(nn,:) =  K.*Q.*Cf.*(Swet/Sref)+CDmisc+CDprot+Rdrag./q;
%%%Finding Max L/D
%   CL([n],:) = CLw([n],:) + CLt([n],:);
    Tr(nn,:) = D(nn,:);         %thrust required
%     figure
%     hold on
%     plot(CD([nn],:), CL([nn],:))
%     plot(CD([nn],:), max(diff(CL)), '--')      %max cl/cd is at the tangent that crosses the zero (steepest tangent?
%     legend('CL versus CD', 'Tangent of L/D')
%     hold off
%%%Propellor
    npmin = 0;
    npmax = 0.725;
    % a = [npmin:npmax];
    np = npmin:(npmax/(length(V)-1)):npmax;      %propellor efficiency
%     PE = 55*10^3;       %Power that the engine provides to the propellor(W)
    Pav = np.*PE;        %Power available
    Tav = Pav./V;       %thrust available
%%%Power Required
    top = 2*(Sref^2).*CD([nn],:).^2;
    bottom = 1./(rho*(CL([nn],:).^3));
    right = (W/Sref)^3;
    Pr = sqrt(top.*bottom.*right);     %power required
%%% Plotting moved to the main
%     figure
%     hold on
%     plot(V,Pr)
%     xlabel('Velocity, V (m/s)')
%     ylabel('Power Required, Pr (W)')
%     hold off
%%% Plotting moved to the main
%%% Plotting Parasitic, Induced, and Total Drag
%     figure
%     hold on
%     plot(V, Di, 'r')
%     plot(V, Do, 'g')
%     plot(V, D, 'b')
%     plot(V, Tav)
%     %%add vertical dashed lines for vmin and vmax (point of intersection)
%     legend('Induced Drag', 'Parasitic Drag', 'Total Drag', 'Thrust Available')
%     xlabel('Velocity, V (m/s)')
%     ylabel('Force, D or T(N)')
%     hold off
end
%%%Total Drag
    CD = CDi + CDo;
    D = q.*CD;
    Di = q.*CDi;
    Do = q.*CDo;
end
