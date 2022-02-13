function Clmaxclean = cleanairfoil(tc,gamma)

    if(gamma == 35)
        Clmaxclean = -338*tc^3 + 88.403*tc^2 - 3.0204*tc + 0.8582;
    else
        y1 = -308.86*tc^3 + 83.566*tc^2 - 2.9097*tc + 0.9141;
        y2 = -338*tc^3 + 88.403*tc^2 - 3.0204*tc + 0.8582;
        Clmaxclean = y1 + ((y2 - y1)/(35-15))*(gamma - 15);
    end

end