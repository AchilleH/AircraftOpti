function [Return] = test(UAV)
    %Tests the input UAV Strucutre if it meets mission specs and can reasonably fly
    %TODO: maybe have this function output the success histograms
    % Test for FAIL if it FAILS return false and return from function, This makes it easier to add/remove conditions. Also more flexible to edit if we want it to be auditable.
    if min(UAV.Sto) > 122 || min(UAV.Sl) > 122 || min(UAV.Sto)==0 || min(UAV.Sl)==0
        Return = false;
        return;
    end
    if max(UAV.Emax) < 2 || max(UAV.Emax)==0
        Return = false;
        return;
    end
    if UAV.S_Margin > 0
        Return = false;
        return;
    end
    Return = true;
end
