function T = getMaxClimbThrust(M,enginetype)
    if(enginetype == 1)
        T = 6547.6*M^2 - 8464.3*M + 12602;
    else
        T = 22500*M^2 - 38607*M + 37843;
    end


end

%% use for segment 3