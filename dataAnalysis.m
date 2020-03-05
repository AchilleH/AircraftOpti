function [] = dataAnalysis(Success,All,name)
figure()
hold
histogram(Success);
histogram(All);
legend('Successful','All');
title(name);
hold
end