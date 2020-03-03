function [Sto,Sl,Emax,Rmax,RCmin,RCmax,gamMin,gamMax,TRmin,Vstall] = Perf(n,C,rho,T,V,D,S,L,k,nu,Clg,Cdg,Cdo,Pav,Preq,W,Vhead)
%Every output is a 1xn vector where n corresponds to the length of the aoa array used to create the other mutidimensional arrays
%will output 0 as Vstall if there is no stall condition at that aoa.
%general Performance notes
%Tr = W/(L/D)
%Pr = Tr*V = D * V
%Tr=D for level flight
%Battery Related Notes
%Rt[hr] is the time that the battery of choice was evaluated over(i.e 1Ah = 1
%amp for 1 hour so Rt is 1 hour.
%n is a discharge parameter specific to each battery and changes throughout
%it's lifespan
%C is battery capacity in ampere hours

%Constants and Preallocation
g = 9.81; %metric g
Vstall = zeros(1,n);
Sto = zeros(1,n);
Sl = zeros(1,n);
Emax = zeros(1,n);
Rmax = zeros(1,n);
RCmin = zeros(1,n);
RCmax = zeros(1,n);
gamMin = zeros(1,n);
gamMax = zeros(1,n);
TRmin = zeros(1,n);

%for loop to handle row by row calculations using matrix inputs
    for i = 1:n %Iterates through the length of aoarange in main(rows)
        i_stall = 1;
        while 1==1
            CompL = L(i,i_stall);
            if (CompL>=W)
                break;
            elseif i_stall==length(V)
                i_stall=1;
                break;
            end
            i_stall = i_stall+1;
        end
        Vstall(i) = V(i_stall);
        Vlof = 1.1*Vstall(i);
        Vtd = 0.7*Vstall(i);
        Va = 1.3*Vstall(i);
        V2 = 1.2*Vstall(i);
        if 1.1*Vstall(i) <= max(V)
            in1 = find(V==round(1.1*Vstall(i))); %gets the index of the closest V for avg takeoff V
        else
            in1 = 1;
        end
        if 1.25*Vstall(i) <= max(V)
            in2 = find(V==round(1.25*Vstall(i))); %same but landing
        else
            in2 = 1;
        end
        if 0.7*Vtd <= max(V)
            in3 = find(V==round(0.7*Vtd)); %for touchdown acceleration assumption
        else
            in3 = 1;
        end
        T1 = T(in1); D1 = D(i,in1);
        T2 = T(in2); D2 = D(i,in2);
        T3 = T(in3); D3 = D(i,in3); L3 = L(i,in3);

        mu = 0.35; %hard braking friction
        q = .5.*rho.*V.^2.*S;
        Rt = 1; %For Ah batteries rated over 1 hour.
        n = 1.3; %Typical Lithium Polymer battery value


        %Climbing Performance
        %maybe take array inputs to get a spectrum and with it our min and max
        %Cl(maxRC) = sqrt(3*Cdo*pi*AR*e); %For Propeller Aircraft
        RC = (Pav(i_stall:end) - Preq(i_stall:end))./W;
        RCmin(i) = min(RC);
        RCmax(i) = max(RC);
        gam = (Pav(i_stall:end) - Preq(i_stall:end))./(W.*V(i_stall:end));
        gamMin(i) = min(gam);
        gamMax(i) = max(gam);

        %To/L
        %Sg is Ground roll for lift off, denominator is evaluated at V = 0.7Vlof
        %Sa is air distance: use average T-D, V2 = 1.2 Vstall
        %Vlof is 1.1*Vmu, Vmu is the min speed aircraft could TO w/o stiking tail
        %on ground
        %TO and landing calc not working, TO goest to inf, landing hovers around 0
        Sg = ((Vlof - Vhead)^2 ./ (2*g*(T1./W - mu - (Cdg(i,i_stall) - mu*Clg(i,i_stall)).*q(in1)./(W./S))));
        Sa = W ./ (T1 - 0.01.*D1) .* floor((V2.^2 - Vlof.^2)./(2.*g) + 15.24); %Relies on faa screen height converted to meters
        Sto(i) = abs(Sg) + abs(Sa);

        %landing, a is avg acceleration at 0.7Vtd
        a = T3-D3-mu.*(W-L3); %use
        Slg = (Vtd - Vhead).^2./(2.*a);
        Sa = W ./ (T2-D2) .* floor((Va.^2 - Vtd.^2)/(2.*g) + 15.24); %15.24 is faa screen height in meters
        Sl(i) = abs(Slg) + abs(Sa);

        %Endurance and Range
        %metric units, km for distance
        %t = Rt/i^n *(c/Rt)^n; %battery endurance
        Emax(i) = Rt^(1-n)*((nu(i_stall)*Vstall(i)*C)/((2/sqrt(rho*S))*Cdo(i,i_stall)^.25*(2*W*sqrt(k/3)^1.5)))^n;
        Rmax(i) = Rt^(1-n)*(nu(i_stall)*Vstall(i)*C/((2/sqrt(rho*S))*Cdo(i,i_stall)^.25*(2*W*sqrt(k)^1.5)))^n * sqrt(2*W/(rho*S) * sqrt(k/Cdo(i,i_stall))) * 3.6;

        %Turning
        %nlc = .5*rho*V^2 * (Clmax/(W/S)); %lift constrained load factor
        n = L(i,i_stall:end)./W; %fyi
        R = 2.*q(i_stall:end)./(g.*rho.*sqrt(n.^2 - 1));
        TRmin(i) = min(R);
    end
end
