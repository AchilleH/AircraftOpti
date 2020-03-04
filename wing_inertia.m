function [Ixw, Iyw, Izw, Ixzw] = wing_inertia(wingweight,type, leadingsweep,trailingsweep,roottaper,tiptaper,chord,wingstart,span,fuselageradius,planetoCG)

%why not just rename the variables in the argument section of the function?
%you can always use comments to descibe what they are instead of writing
%your thesis in the arguments -achille
W_w = wingweight;
Lam_l = leadingsweep;
Lam_t = trailingsweep;
t_r = roottaper;
t_t = tiptaper;
c = chord;
b = span;
wing = 0; % control for calculating wing_inertia for an actual wing
tail = 1; % control for calculating wing_inertia for a tail wing
YS1dot = span-planetoCG; % wing span (one wing) subtracted by the distance from the aircraft fuselage centerline to the wing CG (half span plus fuselage D for rectangular wings)
YS4 = fuselageradius; % y position of wing start ie fuselage radius
XS4 = wingstart; % x position of start of wing

V = b*(t_r(c+b/2*(tan(Lam_t)-tan(Lam_l)))-(t_r-t_t)*(c/2+b/3*(tan(Lam_t)-tan(Lam_l))));

I_1x = W_s*b^3/V*(((t_r-t_t)*(c/4+b*tan(Lam_t)/5-b*tan(Lam_l)/5))+(t_r*(c/3+b*tan(Lam_t)/4-b*tan(Lam_l)/4)));
I_1y = W_s*b/V*((t_r*(c^3/3+b*c*tan(Lam_t)*(c/2+b*tan(Lam_t)/3)+b^3/12*(tan(Lam_t)^3-tan(Lam_l)^3)))-((t_r-t_t)*(c^3/6+b*c*tan(Lam_t)*(c/3+b*tan(Lam_t)/4)+b^3/15*(tan(Lam_t)^3-tan(Lam_l)^3))));
I_1z = I_1x+I_1y;
I_1xz = 0; % != 0 for anhedral or dihedral angles; This isn't even used ?


%This appears to be some control structure that uses the fuselage parallel
%sweep distance compared to chord to set some future varibles up
%also seems as though 1 possibility was not considered, and the else only
%sets up a value for k leaving the C_values to be random garbage in memory.
%-achille
if c < b*tan(Lam_l) && c < b*tan(Lam_l)+c
    C_a = c;
    C_b = b*tan(Lam_l);
    C_c = b*tan(Lam_l)+c;
elseif c > b*tan(Lam_l) && c < b*tan(Lam_l)+c
    C_a = b*tan(Lam_l);
    C_b = c;
    C_c = b*tan(Lam_l)+c;
elseif c > b*tan(Lam_l) && c > b*tan(Lam_l)+c
    C_a = b*tan(Lam_l);
    C_b = b*tan(Lam_l)+c;
    C_c = c;
else
    K_0 = -1;
end

if type == wing
    K_0 = 0.703;
elseif type == tail
    K_0 = 0.771;
else
    K_0 = -1;
end
    
%seems as though, if the else clause in 1st control struct is executed that
%this would evaluate to random garbage -achille
XS1 = (-C_a^2+C_b^2+C_c*C_b+C_c^2)/(3*(C_b+C_c-C_a))*sqrt(K_0);

YS1 = b^2/V*((t_r*(c/2+b/3*(tan(Lam_t)-tan(Lam_l)))-(t_r-t_t)*(c/3+b/4*(tan(Lam_t)-tan(Lam_l)))));

Ixw = I_1x-W_w*(YS1dot)^2+W_w*(YS1dot+YS4)^2;
Iyw = I_1y-W_w*(XS1)^2+W_w*(XS1+XS4)^2;
Izw = I_1z-W_w*(XS1^2+YS1dot^2)+W_w*(XS1+XS4)^2+W_w*(YS1dot+YS4)^2;
Ixzw = 0;

end