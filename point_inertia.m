function [Ix, Iy, Iz, Ixz] = point_inertia(pointweight,Xposition,Yposition,Zposition)

% Point Mass Solid
% W_p = weight of point mass
% x = x-location of point mass
% y = y-location of point mass
% z = z-location of point mass

W_p = pointweight;
x = Xposition;
y = Yposition;
z = Zposition;

Ix = W_p*(y^2+z^2);
Iy = W_p*(x^2+z^2);
Iz = W_p*(x^2+y^2);
Ixz = W_p*x*z;

end
