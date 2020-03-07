function [i_t, delta_e, CL_tot] = it_and_deltae(at, aw, Vstall, Sw, St, W, staticmargin, rho, downwash, aoa, tau)

% function calculates i_t and delta_e for TRIMMED FLIGHT (L = W)
% still needs air density, range for angle of attack, and elevator airfoil
% CM and CL to be user-defined

%at lift-curve slope for tail
%aw lift-curve slope for wing
epsal = downwash; % downwash effect
% static margin
% wing surface area
% tail surface area
v = Vstall; % speed (for L = W)
Wtotal = W; % total aircraft weight

% preliminary calculations
q = .5*rho*v.^2; % dynamic pressure
VH = -staticmargin*St/Sw; % tail volume ratio
CLdelta_e = tau*at*St/Sw; % lift-curve slope for elevator
CMdelta_e = -tau.*at.*VH; % Moment-curve slope for elevator

% deterimining i_t from total lift equation on Lecture 10 Slide 251
i_t = 1./at.*Sw./St.*(Wtotal./(q.*Sw)-(aw+at.*St./Sw.*(1-epsal)).*aoa);

% determining elevator deflection from Lecture 10 Slide 253

% total lift-curve slope
CLa_tot = aw+at.*St./Sw.*(1-epsal); 
% total moment-curve slope
CMa = CLa_tot.*staticmargin;
% lift-curve slope incidence contribution
CLi = -at.*St./Sw;
% moment-curve slope incidence contribution
CMi = at.*VH;
% total coefficent of lift for plane with incidence
CL_tot = CLa_tot.*aoa+CLi.*i_t;
% moment coefficient around plane AC
CMac = 1./CLa_tot.*(CMa.*CLi-CLa_tot.*CMi-CMa.*CL_tot);
% moment coefficient (not aoa dependent)
CM0 = VH.*i_t.*at+CMac;

% elevator deflection
delta_e = -(CM0.*CLa_tot+CMa.*CL_tot)./(CLa_tot.*CMdelta_e-CMa.*CLdelta_e);
