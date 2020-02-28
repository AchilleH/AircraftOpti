function [Ixcg,Iycg,Izcg,Ixzcg] = inertial_calc(Ixarray,Iyarray,Izarray,Ixzarray,weightarray,Xarmarray,Zarmarray)

Ix = 0;
Iy = 0;
Iz = 0;
Ixz = 0;

for i=0:length(Ixarray)
    Ix = Ix+Ixarray(i);
    Iy = Iy+Iyarray(i);
    Iz = Iz+Izarray(i);
    Ixz = Ixz+Ixzarray(i);
end

[xcg, zcg, Wtot] = CG_calc(Xarmarray,Zarmarray,weightarray);

Ixcg = Ix-Wtot*zcg^2;
Iycg = Iy-Wtot*(zcg^2+xcg^2);
Izcg = Iz-Wtot*xcg^2;
Ixzcg = Ixz-Wtot*xcg^2;

end