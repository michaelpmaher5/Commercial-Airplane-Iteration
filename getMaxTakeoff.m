function takeoffthrust = getMaxTakeoff(M,enginetype)
    if(enginetype == 1) 
     takeoffthrust = 7500*M^2 - 9150*M + 14430;
    else
        takeoffthrust = 39286*M^2 - 48214*M + 45486;
    end
end

%% use for segment 1