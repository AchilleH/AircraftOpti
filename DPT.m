function[Wfilters,CDiw,CDit,CDi,tdivc,Q,K,Cf,CDmisc,CDleak,CDprot,CDo,CD,q,D,Di,Do,Tr,np,Pav,Tav,Pr] = DPT(PE,CL,W,Swet,Sref,Aw,Sw,At,St,e,t,c,Y,CLw,CLt,xdivc,h,V)
%Use length of aircraft for c

%%%154A Drag Block
%%%Air Constants
    gamma = 1.4;
    Rair = 287;     %(J/KgK)
    rho = 1.18883546;     %avg. air density (kg/m^3), range of 400-1700ft (121.92-518.16 m) above sea level
    Tair = 339.77598;        %(K), range of 400-1700ft (121.92-518.16 m) above sea level
%%%Aircraft Constants
%     W = ;
    g = 9.81;       %m/s^2
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
    M = V/a;        %mach number
    mew = 1.77907876*10^(-6);     %dynamic viscosity of fluid (N s/m^2)
    v = mew/rho;     %kinematic viscosity (m^2/s)
    Re = V*c/v;       %reynold's number

%%    
%%%Calculating Induced Drag Coefficient
    CDiw = (CLw.^2)/pi/Aw/e;      % CLw = coef. of lift of wing
    CDit = (CLt.^2)/pi/At/e;      % CLt = coef. of lift of tail 
    %%%Total Induced Drag Coefficient, CDi
    CDi = CDiw + (St/Sw)*CDit;
    
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
    fricfrac1 = 0.28;   %prefilter friction factor
    fricfrac2 = 0.6;    %carbon filter friction factor
    fricfrac3 = fricfrac1;  %HyperHepa filter friction factor
    Pair = g*10^3;    %pressure of still air in Force/m^2
    intakearea = 0.3302^2;  %intake area in m^2
    N = Pair*intakearea;    %normal force of still air
    CDleakpercubicmeter = fricfrac1*N+fricfrac2*N+fricfrac3*N;
    FilterSA1 = 2.787;  %PreFilter surface area in m^2
    FilterSA2 = 3.894;  %Carbon filter surface area in m^2
    FilterSA3 = 4.924;  %HyperHepa filter surface area in m^2
    SingleModuleSA = FilterSA1+FilterSA2+FilterSA3; %total surface area of one module (3 air filters) in m^2
    ModNum = 26;    %Num of modules needed to filter 10^6 cf2h based on modules capable of 330cfm
    ModuleSA = SingleModuleSA*ModNum;   %surface area of all filters in the body (m^2)
    CDleak = CDleakpercubicmeter*ModuleSA;  %leakage drag for all modules
    Wprefilter = 1.2224;
    Wcarbonfilter = 2.5;
    Whyperhepafilter = 2.1044;
    Wfilters = ModNum*(Wprefilter+Wprefilter+Whyperhepafilter); %weight of all filters in kg
    %Protuberance Drag, CDprot (due to things that stick out of the body,
    %like the booms for example)
    CDprot = 0;
    %%%Total Parasitic Drag Coefficient, CDo
    CDo =  K.*Q.*Cf.*(Swet/Sref)+CDmisc+CDleak+CDprot; 

    %TEMPORARYY
    %CDo = 0;
%%%Total Drag
    CD = CDi + CDo;
    q = 0.5*rho*Sref*(V.^2);
    D = q.*CD;
    Di = q.*CDi;
    Do = q.*CDo;
    
%%
%%%Finding Max L/D
%     CL = CLw + CLt;
    Tr = D;         %thrust required
    figure
    hold on
    plot(CD, CL)
%     plot(CD, max(diff(CL)), '--')      %max cl/cd is at the tangent that crosses the zero (steepest tangent?
    legend('CL versus CD', 'Tangent of L/D')
    hold off
    
%%%Propellor 
    npmin = 0;
    npmax = 0.725;
    a = [npmin:npmax];
    np = [npmin:(npmax/(length(V)-1)):npmax];      %propellor efficiency
%     PE = 55*10^3;       %Power that the engine provides to the propellor(W)
    Pav = np.*PE;        %Power available
    Tav = Pav./V;       %thrust available
%%%Power Required
    top = 2*(Sref^2).*CD.^2;
    bottom = 1./(rho*(CL.^3));
    right = (W/Sref)^3;
    Pr = sqrt(top.*bottom.*right);     %power required
    figure
    hold on
    plot(V,Pr)
    xlabel('Velocity, V (m/s)')
    ylabel('Power Required, Pr (W)')
    hold off
    
%%% Plotting Parasitic, Induced, and Total Drag
    figure
    hold on
    plot(V, Di, 'r')
    plot(V, Do, 'g')
    plot(V, D, 'b')
    plot(V, Tav)
    %%add vertical dashed lines for vmin and vmax (point of intersection)
    legend('Induced Drag', 'Parasitic Drag', 'Total Drag', 'Thrust Available')
    xlabel('Velocity, V (m/s)')
    ylabel('Force, D or T(N)')
    hold off
end