function tc = getToverCratio(Mdiv,gamma,airfoiltype)

%% Equations from Figure 1A and 1B, chosen depending on airfoiltype and gamma
if(airfoiltype == 1)
    if(gamma > 25 && gamma < 30)
        tc1 = -0.5339*Mdiv + 0.5202;
        tc2 = -0.5098*Mdiv + 0.5095;
        tc = tc1 + ((tc2 - tc1)/(30-25))*(gamma - 25);
    elseif(gamma > 20 && gamma < 25)
        tc1 = -0.5679*Mdiv + 0.5385;
        tc2 = -0.5339*Mdiv + 0.5202;
        tc = tc1 + ((tc2 - tc1)/(25-20))*(gamma - 20);
    elseif(gamma > 30 && gamma < 35)
        tc1 = -0.5098*Mdiv + 0.5095;
        tc2 = -0.4714*Mdiv + 0.4884;
        tc = tc1 + ((tc2 - tc1)/(35-30))*(gamma - 30);
    elseif(gamma > 35 && gamma < 40)
        tc1 = -0.4714*Mdiv + 0.4884;
        tc2 = -0.4268*Mdiv + 0.4632;
        tc = tc1 + ((tc2 - tc1)/(40-35))*(gamma - 35);
    else
        switch gamma
            case 0
                tc = -0.6348*Mdiv + 0.5732;
            case 10
                tc = -0.6161*Mdiv + 0.5635;
            case 15
                tc = -0.5902*Mdiv + 0.5495;
            case 20
                tc = -0.5679*Mdiv + 0.5385;
            case 25
                tc = -0.5339*Mdiv + 0.5202;
            case 30
                tc = -0.5098*Mdiv + 0.5095;
            case 35
                tc = -0.4714*Mdiv + 0.4884;
            case 40
                tc = -0.4268*Mdiv + 0.4632;
        end
    end
    %% Figure 1B equations
else
    if(gamma > 25 && gamma < 30)
        tc1 = -19.565*Mdiv.^3 + 50.237*Mdiv.^2 - 43.52*Mdiv + 12.796;
        tc2 = -16.608*Mdiv.^3 + 44.743*Mdiv.^2 - 40.674*Mdiv + 12.544;
        tc = tc1 + ((tc2 - tc1)/(30-25))*(gamma - 25);
    elseif(gamma > 20 && gamma < 25)
        tc1 = -14.046*Mdiv.^3 + 36.115*Mdiv.^2 - 31.419*Mdiv + 9.314;
        tc2 = -19.565*Mdiv.^3 + 50.237*Mdiv.^2 - 43.52*Mdiv + 12.796;
        tc = tc1 + ((tc2 - tc1)/(25-20))*(gamma - 20);
    elseif(gamma > 30 && gamma < 35)
        tc1 = -16.608*Mdiv.^3 + 44.743*Mdiv.^2 - 40.674*Mdiv + 12.544;
        tc2 = -44.192*Mdiv.^3 + 117.01*Mdiv.^2 - 104.05*Mdiv + 31.167;
        tc = tc1 + ((tc2 - tc1)/(35-30))*(gamma - 30);
    elseif(gamma > 35 && gamma < 40)
        tc1 = -44.192*Mdiv.^3 + 117.01*Mdiv.^2 - 104.05*Mdiv + 31.167;
        tc2 = 13.75*Mdiv.^2 - 25.693*Mdiv + 12.088;
        tc = tc1 + ((tc2 - tc1)/(40-35))*(gamma - 35);
    else
        switch gamma
            case 0
                tc = -36.001*Mdiv.^3 + 84.138*Mdiv.^2 - 66.194*Mdiv + 17.612;
            case 5
                tc = -34.059*Mdiv.^3 + 79.922*Mdiv.^2 - 63.172*Mdiv + 16.9;
            case 10
                tc = -22.468*Mdiv.^3 + 54.439*Mdiv.^2 - 44.541*Mdiv + 12.377;
            case 15
                tc = -19.085*Mdiv.^3 + 47.25*Mdiv.^2 - 39.52*Mdiv + 11.237;
            case 20
                tc = -14.046*Mdiv.^3 + 36.115*Mdiv.^2 - 31.419*Mdiv + 9.314;
            case 25
                tc = -19.565*Mdiv.^3 + 50.237*Mdiv.^2 - 43.52*Mdiv + 12.796;
            case 30
                tc = -16.608*Mdiv.^3 + 44.743*Mdiv.^2 - 40.674*Mdiv + 12.544;
            case 35
                tc = -44.192*Mdiv.^3 + 117.01*Mdiv.^2 - 104.05*Mdiv + 31.167;
            case 40
                tc = 13.75*Mdiv.^2 - 25.693*Mdiv + 12.088;
                %% tc = -583.33*Mdiv.^3 + 1527.5*Mdiv.^2 - 1335*Mdiv + 389.53;
        end
    end
end



end

function dMdiv = getdMdiv(Cl,airfoiltype)
%% Figure 2 Equations
if(airfoiltype == 1)
    dMdiv = -0.0539*Cl.^3 - 0.1824*Cl.^2 - 0.0837*Cl + 0.1104;
else
    dMdiv = 0.8531*Cl.^3 - 1.8037*Cl.^2 + 1.0531*Cl - 0.1755;
end
end

function [ClmaxTO,ClmaxL] = getClmaxAtTOandL(k)
%% Figure 3 Equations
ClmaxTO = 88.229*k.^3 - 65.956*k.^2 + 17.062*k + 1.0456;
ClmaxL = 95.361*k.^3 - 61.901*k.^2 + 15.861*k + 2.0182;

end

function WfWto = getFuelRatio(R)
%% Figure 4 using JT8D estimate
WfWto = (-1.14114471699129)*(10^-16)*R.^4 + 2.43384668306607*(10^-12)*R.^3 - 2.21881698023386*(10^-8)*R.^2 + 0.000143909330264091*R - 0.000902992777565004;
end

function K = getK(tofl,numofengines)
%% Figure 5 Equations
switch numofengines
    case 2
        K = 0.0281*tofl - 8.5091;
    case 3
        K = 0.031*tofl - 4.8909;
    case 4
        K = 0.0319*tofl + 5.0909;
end
end

function netthrust = getNetThrust(M,enginetype)
%% Sea Level Thrust Values at Low Mach Number
%% For use when determining TOFL using JT9D Engines
if(enginetype == 1)
    netthrust = 7500*M^2 - 9150*M + 14430;
else
    if(M == 0)
        netthrust = 45500;
    else
        netthrust = 36778*M.^2 - 46917*M + 45457;
    end
end




end

function solutions = solveForWto(A,B,C,D)
    %%syms x
    %%solutions = vpasolve(A*x^1.195 + B*x^0.235 + C*x + D == 0, x);
    %%k = find(values > 100000 & values < 400000,1);
    %%solutions = values(k);
    x = 50000:50:500000;
    %%checkmatrix = zeros(size(x));
    test = abs(A*x.^1.195 + B*x.^0.235 + C*x + D);
    %%solutions = 500;
    [lookfor,index] = min(test,[],'omitnan');
    solutions = x(index);
    %%for i = 50000:500:900000
       %% test2 = abs(A*i^1.195 + B*i^0.235 + C*i + D);
        %%if(abs(test2 - lookfor) < 10)
          %%  solutions = i;
       %% end
   %% end
    
    
   %% for x = 50000:500:900000
     %%   test = A*x^1.195 + B*x^0.235 + C*x + D;
       %% checkmatrix(x,2) = abs(test);
        %%if(abs(test) < 1000)
           %% solutions = x;
       %% end
        
    %%end
    
    
end

function Cf = getCf(RNk)

Cf = 0.0756*(RNk)^-0.19;

end

function K = getFormFactor(M0,gamma,tc)
    Z = ((2-M0^2)*cosd(gamma))/(((1-M0^2)*(cosd(gamma))^2)^0.5);
    K = (1 + Z*tc + 100*tc^4);
end

function K = getFuselageFormFactor(LoD)

    K = -0.0019*LoD^3 + 0.0477*LoD^2 - 0.4157*LoD + 2.4156;
end

function [Ta,sfc] = getMaxClimbThrustandSFC(M,enginetype)

if enginetype == 1
    Ta = -928.57*M + 6917.9;
    sfc = 0.4357*M + 0.4782;
else
    T1 = 8571.4*M^2 - 19914*M + 27263;
    sfc1 = 0.44*M + 0.308;
    T2 = -7142.9*M^2 + 10729*M + 10184;
    sfc2 = 0.37*M + 0.341;
    Ta = (T1+T2)/2;
    sfc = 1.1*((sfc1+sfc2)/2);
end


end

function c = getC(T,enginetype)
    
if enginetype == 1
    c = -0.0000000000123894176*T^3 + 0.00000012928744101277*T^2 - 0.000432156078048275*T + 1.2436537595825;
else
    c = -0.00000000000457949485*T^3 + 0.00000007376752381776*T^2 - 0.00040729645082369*T + 1.4091620718046;
    c = 1.1*c;
end

end

function dCdp = figure6(CloverClmax,segment)

    if segment == 1
        dCdp = 0.1107*CloverClmax^3 - 0.0437*CloverClmax^2 - 0.0483*CloverClmax + 0.0321;
    else
        dCdp = 0.0821*CloverClmax^3 + 0.003*CloverClmax^2 - 0.066*CloverClmax + 0.0408;
    end

end

function Clmaxclean = cleanairfoil(tc,gamma)

    if(gamma == 35)
        Clmaxclean = -338*tc^3 + 88.403*tc^2 - 3.0204*tc + 0.8582;
    else
        y1 = -308.86*tc^3 + 83.566*tc^2 - 2.9097*tc + 0.9141;
        y2 = -338*tc^3 + 88.403*tc^2 - 3.0204*tc + 0.8582;
        Clmaxclean = y1 + ((y2 - y1)/(35-15))*(gamma - 15);
    end

end

function takeoffthrust = getMaxTakeoff(M,enginetype)
    if(enginetype == 1) 
     takeoffthrust = 7500*M^2 - 9150*M + 14430;
    else
        takeoffthrust = 39286*M^2 - 48214*M + 45486;
    end
end

%% use for segment 1

function T = getMaxClimbThrust(M,enginetype)
    if(enginetype == 1)
        T = 6547.6*M^2 - 8464.3*M + 12602;
    else
        T = 22500*M^2 - 38607*M + 37843;
    end


end

%% use for segment 3

