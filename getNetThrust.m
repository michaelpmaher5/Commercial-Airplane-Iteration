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