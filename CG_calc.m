function [XCG,ZCG,Wtotal] = CG_calc(Xarmarray,Zarmarray,weightarray)

XCG = 0;
ZCG = 0;
Wtotal = 0;
MXtotal = 0;
MZtotal = 0;

for i=0:length(cgarray)
    Wtotal = Wtotal+weightarray(i);
    MXtotal = MXtotal+weightarray(i)*Xarmarray(i);
    MZtotal = MZtotal+weightarray(i)*Zarmarray(i);
end

XCG = MXtotal/Wtotal;
ZCG = MZtotal/Wtotal;

end