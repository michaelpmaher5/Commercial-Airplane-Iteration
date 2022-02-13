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