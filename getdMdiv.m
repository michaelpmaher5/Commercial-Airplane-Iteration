function dMdiv = getdMdiv(Cl,airfoiltype)
%% Figure 2 Equations
if(airfoiltype == 1)
    dMdiv = -0.0539*Cl.^3 - 0.1824*Cl.^2 - 0.0837*Cl + 0.1104;
else
    dMdiv = 0.8531*Cl.^3 - 1.8037*Cl.^2 + 1.0531*Cl - 0.1755;
end
end