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