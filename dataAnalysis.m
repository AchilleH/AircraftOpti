function [] = dataAnalysis(Success,All,name)
figure()
hold
histogram(Success,80);
histogram(All,80);
legend('Successful','All');
title(name);
hold
end