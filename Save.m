function [Data] = Save(Sw,St,CLwa,CLta,CLw,CLt,CL,CLmax,Lw,Lt,bw,bt,Aw,At,ew,et,CDi,CDo,CD,D,Di,Do,Tr,np,Pav,Tav,Pr,Sto,Sl,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,Rmin,Vstall,XCG,ZCG,Wtotal,hn,staticmargin)
%% Lift Variables and Output
Data.Sw = Sw;
Data.St = St;
Data.CLwa = CLwa;
Data.CLta = CLta;
Data.CLw = CLw;
Data.CLt = CLt;
Data.CL = CL;
Data.CLmax = CLmax;
Data.Lw = Lw;
Data.Lt = Lt;
Data.bw = bw;
Data.bt = bt;
Data.Aw = Aw;
Data.At = At;
Data.ew = ew;
Data.et = et;
%% Drag
Data.CDi = CDi;
Data.CDo = CDo;
Data.CD = CD;
Data.D = D;
Data.Di = Di;
Data.Do = Do;
Data.Tr = Tr;
Data.np = np;
Data.Pav = Pav;
Data.Tav = Tav;
Data.Pr = Pr;
%% Performance
Data.Sto = Sto;
Data.Sl = Sl;
Data.Emax = Emax;
Data.Rmax = Rmax;
Data.RCmin = RCmin;
Data.RCmax = RCmax;
Data.gamMin = gamMin;
Data.gamMax = gamMax;
Data.Rmin = Rmin;
Data.Vstall = Vstall;
%% CG
Data.Xcg = XCG;
Data.Zcg = ZCG;
Data.Wtotal = Wtotal;
%% Stability
Data.hn = hn;
Data.S_Margin = staticmargin;
end
