function [Ixe, Iye, Ize, Ixze] = engine_inertia(engineweight,nacelleradius,XZtoengineCG,XYtoengineCG,enginelength)

W_e = engineweight;
R_e = nacelleradius;
YP = XZtoengineCG;
ZP = XYtoengineCG;
l_e = enginelength;

Ixe = W_e*R_e^2/2+W_e*(YP^2+ZP^2);
Iye = W_e/12*(3*R_e^2+l_e^2)+W_e*(XP^2+ZP^2);
Ize = W_e/12*(3*R_e^2+l_e^2)+W_e*(XP^2+YP^2);

Ixze = W_e*XP*ZP;

end