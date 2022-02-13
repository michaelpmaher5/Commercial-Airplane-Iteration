function [ClmaxTO,ClmaxL] = getClmaxAtTOandL(k)
%% Figure 3 Equations
ClmaxTO = 88.229*k.^3 - 65.956*k.^2 + 17.062*k + 1.0456;
ClmaxL = 95.361*k.^3 - 61.901*k.^2 + 15.861*k + 2.0182;

end