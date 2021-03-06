function [Ixw, Iyw, Izw, Ixzw] = wing_inertia(wingweight,type, leadingsweep,trailingsweep,roottaper,tiptaper,chord,wingstart,span,fuselageradius,planetoCG)


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

V = b*(t_r*(c+b/2*(tan(Lam_t)-tan(Lam_l)))-(t_r-t_t)*(c/2+b/3*(tan(Lam_t)-tan(Lam_l))));

I_1x = W_w*b^3/V*(((t_r-t_t)*(c/4+b*tan(Lam_t)/5-b*tan(Lam_l)/5))+(t_r*(c/3+b*tan(Lam_t)/4-b*tan(Lam_l)/4)));
I_1y = W_w*b/V*((t_r*(c^3/3+b*c*tan(Lam_t)*(c/2+b*tan(Lam_t)/3)+b^3/12*(tan(Lam_t)^3-tan(Lam_l)^3)))-((t_r-t_t)*(c^3/6+b*c*tan(Lam_t)*(c/3+b*tan(Lam_t)/4)+b^3/15*(tan(Lam_t)^3-tan(Lam_l)^3))));
I_1z = I_1x+I_1y;
I_1xz = 0; % != 0 for anhedral or dihedral angles; This isn't even used ?




temp1 = b*tan(Lam_l);
temp2 = b*tan(Lam_l)+c;

if c <= temp1 && c <= temp2
    C_a = c;
    C_b = b*tan(Lam_l);
    C_c = b*tan(Lam_l)+c;
elseif c >= temp1 && c < temp2
    C_a = b*tan(Lam_l);
    C_b = c;
    C_c = b*tan(Lam_l)+c;
elseif c >= temp1 && c >= temp2
    C_a = b*tan(Lam_l);
    C_b = b*tan(Lam_l)+c;
    C_c = c;
else
    K_0 = -1;
end

% different constants if it's a wing or tail
if type == wing
    K_0 = 0.703;
elseif type == tail
    K_0 = 0.771;
else
    K_0 = -1;
end
    

XS1 = (-C_a^2+C_b^2+C_c*C_b+C_c^2)/(3*(C_b+C_c-C_a))*sqrt(K_0);

YS1 = b^2/V*((t_r*(c/2+b/3*(tan(Lam_t)-tan(Lam_l)))-(t_r-t_t)*(c/3+b/4*(tan(Lam_t)-tan(Lam_l)))));

Ixw = I_1x-W_w*(YS1dot)^2+W_w*(YS1dot+YS4)^2;
Iyw = I_1y-W_w*(XS1)^2+W_w*(XS1+XS4)^2;
Izw = I_1z-W_w*(XS1^2+YS1dot^2)+W_w*(XS1+XS4)^2+W_w*(YS1dot+YS4)^2;
Ixzw = 0;

end