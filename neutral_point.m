function [hn] = neutral_point(wingac, tailac, tailarea, wingarea, tailliftcurve, wingliftcurve, downwasheffect)

hacw = wingac;
hact = tailac;
St = tailarea;
Sw = wingarea; 
at = tailliftcurve; % tail lift curve slope
aw = wingliftcurve; % wing lift curve slope
epsal = downwasheffect;

% neutral point calculation from slides (dimensionless)
hn = (hacw+hact*(St/Sw)*(at/aw)*(1-epsal))/(1+(St/Sw)*(at/aw)*(1-epsal));

end
