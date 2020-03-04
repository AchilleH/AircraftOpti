function [Ixcg,Iycg,Izcg,Ixzcg] = inertial_calc(Ixarray,Iyarray,Izarray,Ixzarray,weightarray,Xarmarray,Zarmarray)

Ix = 0;
Iy = 0;
Iz = 0;
Ixz = 0;

for i=1:length(Ixarray)
    Ix = Ix+Ixarray(i); % adding up Ix arrays
    Iy = Iy+Iyarray(i); % adding up Iy arrays
    Iz = Iz+Izarray(i); % adding up Iz arrays
    Ixz = Ixz+Ixzarray(i); % adding up Ixz arrays
end

[xcg, zcg, Wtot] = CG_calc(Xarmarray,Zarmarray,weightarray);

% recentering inertias around the CG
Ixcg = Ix-Wtot*zcg^2; 
Iycg = Iy-Wtot*(zcg^2+xcg^2);
Izcg = Iz-Wtot*xcg^2;
Ixzcg = Ixz-Wtot*xcg^2;

end