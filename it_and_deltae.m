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

for I = 1:length(aoa)
    for J = 1:length(v)

% preliminary calculations
q = .5*rho*v(J).^2; % dynamic pressure
VH = -staticmargin*St/Sw; % tail volume ratio
CLdelta_e(I,:) = tau*at*St/Sw; % lift-curve slope for elevator
CMdelta_e(I,:) = -tau.*at.*VH; % Moment-curve slope for elevator

% deterimining i_t from total lift equation on Lecture 10 Slide 251
i_t(I,:) = 1./at.*Sw./St.*(Wtotal./(q.*Sw)-(aw+at.*St./Sw.*(1-epsal)).*aoa(I));

% determining elevator deflection from Lecture 10 Slide 253

% total lift-curve slope
CLa_tot(I) = aw+at.*St./Sw.*(1-epsal); 
% total moment-curve slope
CMa(I,:) = CLa_tot(I).*staticmargin;
% lift-curve slope incidence contribution
CLi(I,:) = -at.*St./Sw;
% moment-curve slope incidence contribution
CMi(I,:) = at.*VH;
% total coefficent of lift for plane with incidence
CL_tot(I) = CLa_tot(I,:).*aoa(I)+CLi(I,:).*i_t(I,:);
% moment coefficient around plane AC

CMac(I,:) = 1./CLa_tot(I,:).*(CMa(I,:).*CLi(I,:)-CLa_tot(I,:).*CMi(I,:)-CMa(I,:).*CL_tot(I));
% moment coefficient (not aoa dependent)
CM0(I,:) = VH.*i_t(I,:).*at+CMac(I,:);

% elevator deflection
delta_e(I,:) = -(CM0(I,:).*CLa_tot(I,:)+CMa(I,:).*CL_tot(I))./(CLa_tot(I,:).*CMdelta_e(I,:)-CMa(I,:).*CLdelta_e(I,:));
    end
end
end