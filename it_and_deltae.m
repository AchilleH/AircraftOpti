function [i_t, delta_e] = it_and_deltae(CLwa, CLta, Sw, St, rho, v, Wtotal, downwash, aoa, staticmargin, CMdelta_e, CLdelta_e)

% function calculates i_t and delta_e for TRIMMED FLIGHT (L = W)

at = CLta; % lift-curve slope for tail
q = .5*rho*v^2; % dynamic pressure
aw = CLwa; % lift-curve slope for wing
epsal = downwash; % downwash effect
VH = -staticmargin*St/Sw; % tail volume ratio

% deterimining i_t from total lift equation on Lecture 10 Slide 251

i_t = 1/at*Sw/St*(Wtotal/(q*Sw)-(aw+at*St/Sw*(1-epsal))*aoa);

% determining elevator deflection from Lecture 10 Slide 253 (needs a lot of
% intermediate calculations, NEEDS CM AND CL for elevator to be user-defined)

% total lift-curve slope
CLa_tot = aw+at*St/Sw*(1-epsal); 
% total moment-curve slope
CMa = CLa_tot*staticmargin;
% lift-curve slope incidence contribution
CLi = -at*St/Sw;
% moment-curve slope incidence contribution
CMi = at*VH;
% total coefficent of lift for plane with incidence
CL_tot = CLa_tot*aoa+CLi*i_t;
% moment coefficient around plane AC
CMac = 1/CLa_tot*(CMa*CLi-CLa_tot*CMi-CMa*CL_tot);
% moment coefficient (not aoa dependent)
CM0 = VH*i_t*at+CMac;

% elevator deflection
delta_e = -(CM0*CLa_tot+CMa*CL_tot)/(CLa_tot*CMdelta_e-CMa*CLdelta_e);

end
