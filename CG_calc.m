function [XCG,ZCG,Wtotal] = CG_calc(Xarmarray,Zarmarray,weightarray)

Wtotal = 0;
MXtotal = 0;
MZtotal = 0;

for i=0:length(weightarray)
    Wtotal = Wtotal+weightarray(i); % sum all the weights
    MXtotal = MXtotal+weightarray(i)*Xarmarray(i); % sum all the x direction moment arms
    MZtotal = MZtotal+weightarray(i)*Zarmarray(i); % sum all the z direction moment arms
end

XCG = MXtotal/Wtotal; % divide x moment arms by total weight to get CG x position
ZCG = MZtotal/Wtotal; % divide z moment arms by total weight to get CG z position

end