function WfWto = getFuelRatio(R)
%% Figure 4 using JT8D estimate
WfWto = (-1.14114471699129)*(10^-16)*R.^4 + 2.43384668306607*(10^-12)*R.^3 - 2.21881698023386*(10^-8)*R.^2 + 0.000143909330264091*R - 0.000902992777565004;
end