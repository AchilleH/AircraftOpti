function [Ix, Iy, Iz, Ixz] = box_inertia(volumeweight,volumelength,averageradius,Zposition,Xposition)

W_vo = volumeweight;
l_v = volumelength;
R_v = averageradius;
z = Zposition;
x = Xposition;

Ix = W_vo/12*(2*l_v^2+2*R_v^2)+W_vo*z^2;
Iy = W_vo/12*(l_v^2+2*R_v^2)+W_vo*(x^2+z^2);
Iz = W_vo/12*(l_v^2+2*R_v^2)+W_vo*(x^2);
Ixz = W_vo*x*z;

end