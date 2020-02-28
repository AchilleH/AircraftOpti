function [hn] = neutral_point(wingac, tailac, tailarea, wingarea, tailliftcurve, wingliftcurve, downwasheffect)

hacw = wingac;
hact = tailac;
St = tailarea;
Sw = wingarea;
at = tailliftcurve;
aw = wingliftcurve;
epsal = downwasheffect;

hn = (hacw+hact*(St/Sw)*(at/aw)*(1-epsal))/(1+(St/Sw)*(at/aw)*(1-epsal));

end
