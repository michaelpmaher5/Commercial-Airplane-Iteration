function K = getFormFactor(M0,gamma,tc)
    Z = ((2-M0^2)*cosd(gamma))/(((1-M0^2)*(cosd(gamma))^2)^0.5);
    K = (1 + Z*tc + 100*tc^4);
end