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