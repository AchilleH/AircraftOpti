function [Sw,St,CLwa,CLta,CLw,CLt,CL,CLe,CLmax,CLemax,Lw,Lt,Le,bw,bt,ARw,ARt,ew,et] = lift(rho,Clwa,Clta,Clwo,Clto,Sw,St,Sref,bw,bt,ARw,ARt,Df,tapw,tapt,phiw,phit,Vrange,aoa,tau,ealph)
%Lift Calculation Block,[W,Sw,St,CLwa,CLta,Lw,Lt,bw,bt,ARw,ARt,ew,et] = lift(rho,Clwa,Clta,ew,et,W,Sw,St,bw,bt,ARw,ARt,Df,tapw,tapt,phiw,phit,Vrange)  Leave a 0 in 1 argument to have it return the design
%value. For CL leave a 0 in Cl to have it generate CL. Df is fuselage
%diameter and will be used to approximate e if e == 0, along with
%tap(taper). phi is sweep angle measured from the horizontal. Assiming V of < 333 ft/s and general aviation
%If e cannot be calculated, will assume a rectangular wing where e == 0.7.
%Leave a 0 in Cla to have it not generate a lift polar.

%e = 0.7 for rect, 1 for ellipse or taper of ~0.3-0.4

%% Assigning e
if(bw ~= 0 && ARw ~= 0 && Df ~= 0)
    dtap = -0.357 + 0.45*exp(0.0375*phiw);
    tap = tapw - dtap;
    f = 0.0524*tap^4 - 0.15*tap^3 + 0.1659*tap^2 - 0.0706*tap + 0.0119;
    eth = 1/(1+f*ARw);
    kf = 1 - 2*(Df/bw)^2;
    ew = abs(eth * kf * 0.804);
else
    ew = 0.7;
end

if(bt ~= 0 && ARt ~= 0 && Df ~= 0)
    dtap = -0.357 + 0.45*exp(0.0375*phit);
    tap = tapt - dtap;
    f = 0.0524*tap^4 - 0.15*tap^3 + 0.1659*tap^2 - 0.0706*tap + 0.0119;
    eth = 1/(1+f*ARt);
    kf = 1 - 2*(Df/bt)^2;
    et = abs(eth * kf * 0.804);
else
    et = 0.7;
end

%% Getting 3d CLa and Plotting CLa over specified range
%using equation from design textbook
%beta = 1-M^2; % M is
%nu = Cla/(2*pi/beta)
M = Vrange./343; %Mach number were c is in m/s
beta = 1-M.^2;
nuw = Clwa./(2.*pi./beta);
nut = Clta./(2.*pi./beta);
CLwa = getCL(beta,nuw,phiw,Sw,Sref,Df,bw,ARw);
CLta = getCL(beta,nut,phit,St,Sref,Df,bt,ARt);

%% running lift Equation calculations
qw = .5 .* Vrange.^2 .* rho .* Sw; %L = .5 * rho * v^2 * S * CL
qt = .5 .* Vrange.^2 .* rho .* St;
Clew = tau * ealph * Clwa; %Elevator effect on CL at max deflection for wing
Clet = tau * ealph * Clta; %uses lowercase l but for 3d case, assumes the deflection is essentially = to deflection of wing
CLw = meshgrid(CLwa,aoa);
CLt = meshgrid(CLta,aoa);
CLew = meshgrid(CLwa,aoa);
CLet = meshgrid(CLta,aoa);
Lw = meshgrid(qw,aoa);
Lew = meshgrid(qw,aoa);
Lt = meshgrid(qt,aoa);
Let = meshgrid(qt,aoa);

    for I = 1:length(aoa)
     CLw(I,:) = CLw(I,:) * aoa(I) + Clwo;
     CLt(I,:) = CLt(I,:) * aoa(I) + Clto;
     CLew(I,:) = CLw(I,:) + Clew;
     CLet(I,:) = CLt(I,:) + Clet;
     Lew(I,:) = Lw(I,:) .* CLew(I,:);
     Lw(I,:) = Lw(I,:) .* CLw(I,:);
     Let(I,:) = Lt(I,:) .* CLet(I,:);
     Lt(I,:) = Lt(I,:) .* CLw(I,:);
    end
    Le = Lew + Let;
    S = Sw+St;
    CLwmax = CLw(length(aoa),:);
    CLtmax = CLt(length(aoa),:);
    CLewmax = CLew(length(aoa),:);
    CLetmax = CLet(length(aoa),:);
    CLmax = Sw/S*CLwmax + St/S*CLtmax;
    CLemax = Sw/S*CLewmax + St/S*CLetmax;
    CL = Sw/S*CLw + St/S*CLt;
    CLe = Sw/S*CLew + St/S*CLet;
end


function [CLa] = getCL(beta,nu,phi,S,Sref,df,b,A)
F = 1.07*S/Sref*(1+df/b)^2;
CLa = (2.*pi.*A).*F./(2+sqrt(4 + (A^2.*beta.^2./nu.^2).*((1+tan(phi)^2)./beta.^2)));
end

%% Retired Code
% OLD CLa CALCULATION METHOD
% if(Clwa ~= 0 && ew ~= 0)
%     CLwa = Clwa/(1+(Clwa/(pi*ARw*ew)));
%     PCLa = CLwa * aoarange;
%     figure(1)
%     plot(aoarange,PCLa)
%     title('Lift Polar: Wing');
%     xlabel('Angle of Attack [deg]');
%     ylabel('CL');
% end
% if(Clw ~= 0 && Clwa ~= 0)
%     a = Clw/Clwa;
%     CLw = Clwa * a;
% end
% %
% if(Clta ~= 0 && et ~= 0)
%     CLta = Clta/(1+(Clta/(pi*ARt*et)));
%     PCLa = CLta * aoarange;
%     figure(2)
%     plot(aoarange,PCLa)
%     title('Lift Polar: Tail');
%     xlabel('Angle of Attack [deg]');
%     ylabel('CL');
% end
% if(Clt ~= 0 && Clta ~= 0)
%     a = Clt/Clta;
%     CLt = Clta * a;
% end
%OLD LIFT EQUATION CALCULATIONS
% elseif(Sw == 0)
%     Sw = W./(.5.*rho.*CLw.*Vrange.^2) - St*CLt/CLw/Sw;
% elseif(St == 0)
%     St = W/(.5*rho*CLt*Vstall^2) - Sw*CLw/CLt/St;
% % elseif(Clw == 0)
% %     Clw = W/(.5*rho*Clw*Vstall^2) - St*CLt/Sw;
% % elseif(Clt == 0)
% %     Clt = W/(.5*rho*CLt*Vstall^2) - Sw*CLw/St;
% elseif(W==0)
%     W = .5*rho*Vstall^2*(Sw*CLw + St*CLt/Sw);
%
%  Getting CL values as they change with respect to the velocity range given,
% %-the effect Re has on CL for wing, tail and body
% %Since we're flying below 0.3 mach we should be abe to assume a constant CL
% %just based on the aoa of the wing(s)
% CLt = CLta * aoat;
% CLw = CLwa * aoaw;
%
% CLa = zeros(size(Vrange));
% CLa(1:end) = CLw + St*CLt/Sw;
% %CLa = [0:0.05:0.05*47]; uncomment to use the old CL array that produced
% %somewhat ok drag polars
%
% L = .5 .* CLa .* Vrange.^2 .* rho .* (Sw+St);
% figure(3)
% plot(Vrange,L)
% title('Lift Vs Velocity')
% xlabel('Velocity [m/s]')
% ylabel('Lift [N]')
% end
