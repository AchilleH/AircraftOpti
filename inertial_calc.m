function [Ixcg,Iycg,Izcg,Ixzcg] = inertial_calc(Ixtotal,Iytotal,Iztotal,Ixztotal,weightarray,Xarmarray,Zarmarray)


[xcg, zcg, Wtot] = CG_calc(Xarmarray,Zarmarray,weightarray);

% recentering inertias around the CG
Ixcg = Ixtotal-Wtot*zcg^2; 
Iycg = Iytotal-Wtot*(zcg^2+xcg^2);
Izcg = Iztotal-Wtot*xcg^2;
Ixzcg = Ixztotal-Wtot*xcg^2;

end