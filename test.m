function [Return] = test(UAV)
    %Tests the input UAV Strucutre if it meets mission specs and can reasonably fly
    %TODO: maybe have this function output the success histograms
    % Test for FAIL if it FAILS return false and return from function, This makes it easier to add/remove conditions. Also more flexible to edit if we want it to be auditable.
    if UAV.Sto > 122 || UAV.Sl > 122
        Return = false;
        return;
    end
    if UAV.Emax < 2
        Return = false;
        return;
    end

end
