clc
clear all

%% Design Parameters
passengers = 150;
cargoweight = 3000; %% in lbs
range = 3000; %% in naut miles
tofl = 6000; %% take off field length in ft
las = 135; %% landing approach speed in knots, design payload plus 80% max fuel
altitude = 35000; %% in feet
x = 0.2; %% should be 1 - percentage of fuel at landing (0.2)

%% Additional Parameters
M = 0.8; %% Mach number
sigma = 0.953; %% use ideal gas law, then get ratio
delta = 0.23; %% from chart in Shevell
lambda = 0.35; %% wing taper ratio

%%abreast = 4; %% loop
%%aisle = 1; %% loop

WtoT = 0;
counter = 0;
counter2 = 0;
counter3 = 0;
counter4 = 0;
counter5 = 0;
counter6 = 0;
limit = 50;
limitbreakcount = 0;

limit5 = 50;
limit5breakcount = 0;

limit6 = 50;
limit6breakcount = 0;
%%Vcruise = M*576.4; %% in knots
%%Rao = range + 200 + 0.75*Vcruise;
R = 0;
previousR = 0;

thrustpass = 0;
seg1pass = 0;
seg2pass = 0;
seg3pass = 0;
approachsegpass = 0;
landingsegpass = 0;

DOC = 0;
previousDOC = 1000;

sweep1 = 0;
AR1 = 0;
airfoiltype1 = 0;
numofengines1 = 0;
enginetype1 = 0;
enginemountlocation1 = 0;
abreast1 = 0;
aisle1 = 0;

previousDOC2 = 10000;
sweep2 = 0;
AR2 = 0;
airfoiltype2 = 0;
numofengines2 = 0;
enginetype2 = 0;
enginemountlocation2 = 0;
abreast2 = 0;
aisle2 = 0;

designmatrix = zeros(21,7,2,3,2,2,5,2,20); %%sweep,AR,airfoil,numofengines,enginetype,enginemountloaction,abreast,aisle,(first Cl second W/T third Wtakeoff fourth f fifth e sixth R seventh wfwtoJT8D eigth Wfuel ninth Trequired tenth DOC eleventh DOCpaxmile twelth thicktochord thirteenth total drag)

%% Loop Checks
hasReachedClimbLoop = 0;
hasReachedGradLoop = 0;

neverconverged = 0;
converged = 0;
%% Step 1
%%sweep = 30;
%%AR = 6;
%airfoiltype = 1; %% 1 for conventional, 2 for supercritical
%%enginetype = 1; %% 1 for JT8D, 2 for JT9D
%%enginemountlocation = 1; %% 1 for wing mounted, 2 for fuselage mounted

for sweep = 20:1:38 %%20-40
    for AR = 6:1:12 %%6-12
        for airfoiltype = 1:1:2 %% 1-2
            for numofengines = 2:1:2 %%iterate from 2-4
                for enginetype = 1:1:2
                    for enginemountlocation = 1:1:1 %%1-2
                        for abreast = 6:1:6 %%4-8
                            for aisle = 1:1:1 %% 1-2
                                Vcruise = M*576.4; %% in knots
                                Rao = range + 200 + 0.75*Vcruise;
                                wfwtoJT8D = getFuelRatio(Rao);
                                
                                if(enginetype == 2)
                                    wfwtoJT8D = wfwtoJT8D*(0.61/0.78); %% if using JT9D engine
                                end
                                
                                %%Treqcheck = 1;
                                %%Tavailable = 0;
                                
                                hasReachedGradLoop = 0;
                                thrustpass = 0;
                                seg1pass = 0;
                                seg2pass = 0;
                                seg3pass = 0;
                                approachsegpass = 0;
                                landingsegpass = 0;
                                
                                while(thrustpass + seg1pass + seg2pass + seg3pass + approachsegpass + landingsegpass ~= 6)
                                    
                                    while(abs(R - Rao) > 50)
                                        %% Step 2
                                        if(hasReachedGradLoop == 0)
                                            
                                            Cl = 0.55; %% assumption
                                            Clinitialcruise = 0;
                                            
                                            while(Clinitialcruise < Cl - 0.005 || Clinitialcruise > Cl + 0.005)
                                                deltaMdiv = getdMdiv(Cl,airfoiltype); %% from figure 2
                                                
                                                %% Step 3
                                                Mdiv = (M + 0.004)-deltaMdiv;
                                                
                                                %% Step 4
                                                thicktochord = getToverCratio(Mdiv,sweep,airfoiltype); %% from figure 1
                                                
                                                %% Step 5
                                                param = ((cosd(sweep))^2) *((thicktochord)^2)*AR;
                                                [ClmaxTO,ClmaxL] = getClmaxAtTOandL(param);
                                                
                                                %% Step 6
                                                landingwingloading = ((las/1.3).^2)*(sigma*ClmaxL/296);
                                                
                                                %% Step 7
                                                %%Vcruise = M*576.4; %% in knots
                                                %%Rao = range + 200 + 0.75*Vcruise;
                                                
                                                %%while(Rao > R + 50 || Rao < R - 50)
                                                %% Step 8
                                                %%wfwtoJT8D = getFuelRatio(Rao);
                                                
                                                %% Step 9
                                                %%if(enginetype == 2)
                                                %%wfwtoJT8D = wfwtoJT8D*(0.61/0.78); %% if using JT9D engine
                                                %%end
                                                %% Step 10
                                                takeoffwingloading = landingwingloading/(1 - x*wfwtoJT8D);
                                                
                                                %% Step 11
                                                initialcruisewingloading = 0.965*takeoffwingloading;
                                                
                                                %% Step 12
                                                Clinitialcruise = initialcruisewingloading/(1481*delta*M.^2);
                                                counter = counter + 1;
                                                if(Clinitialcruise < Cl - 0.005 || Clinitialcruise > Cl + 0.005)
                                                    Cl = (Cl + Clinitialcruise)/2;
                                                    Clinitialcruise = 0;
                                                end
                                            end
                                            
                                            %% TOFL
                                            
                                            K = getK(tofl,numofengines); %% figure 5 equation
                                            WtoTpoint7LO = K*sigma*ClmaxTO/takeoffwingloading;
                                            Vliftoff = 1.2*((296*takeoffwingloading)/(sigma*ClmaxTO)).^0.5;
                                            Mliftoff = (Vliftoff/661)/((sigma).^0.5);
                                            Mliftoffseventypercent = 0.7*Mliftoff;
                                            Tslst = getNetThrust(0,enginetype);
                                            Tseventypercent = getNetThrust(Mliftoffseventypercent,enginetype);
                                            WtoT = WtoTpoint7LO*(Tseventypercent/Tslst);
                                            designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,1) = Clinitialcruise;
                                            designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,12) = thicktochord;
                                        end
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,2) = WtoT;
                                        
                                        %% Weight Calculations
                                        
                                        Kw = 1;
                                        if(enginemountlocation == 2)
                                            Kw = 1.03;
                                        end
                                        tcave = thicktochord + 0.03;
                                        Ww = (0.00945*(AR^0.8)*((1+lambda)^0.25)*Kw*((1.5*2.5)^0.5))/(((tcave)^0.4)*cosd(sweep)*takeoffwingloading^0.695);%% A
                                        
                                        Kfus = 11.5; %% Kfus = 11.5 for passengers > 135
                                        
                                        l = 3.76*(passengers/abreast) + 33.2;
                                        d = 1.75*abreast + 1.58*aisle + 1;
                                        Wfus = 0.6727*Kfus*(l^0.6)*(d^0.72)*((3.75)^0.3); %% B
                                        
                                        % w landing gear = 0.04wto (C)
                                        
                                        Wnp = 0.0555/WtoT; %% C
                                        
                                        if(enginetype == 2)
                                            Wnp = 1.1*Wnp;
                                        end
                                        
                                        Kts = 0.17;
                                        if(enginemountlocation == 2)
                                            Kts = 0.25;
                                        end
                                        Wts = Kts*Ww; %% A
                                        
                                        Wpp = 1/(3.58*WtoT); %% C
                                        if(enginetype == 2)
                                            Wpp = 1.1*Wpp;
                                        end
                                        
                                        Wfuel = 1.0275*wfwtoJT8D; %% C
                                        
                                        Wpayload = 215*passengers + cargoweight; %% D
                                        
                                        WfeC = 0.035; %% C
                                        Nflightcrew = 2;
                                        Nstew = 3;
                                        WfeD = 132*passengers + 300*numofengines + 260*Nflightcrew + 170*Nstew;%% D
                                        
                                        if(airfoiltype == 2) %% adjusting for composite structure
                                            Ww = 0.7*Ww;
                                            Wts = 0.7*Wts;
                                            Wfus = 0.85*Wfus;
                                            Wnp = 0.8*Wnp;
                                            WfeC = 0.9*WfeC;
                                            WfeD = ((WfeD - (260*Nflightcrew))*0.9) + (260*Nflightcrew);
                                        end
                                        %%Total Weight
                                        A = Ww + Wts;
                                        B = Wfus;
                                        C = 0.04 + Wnp + Wpp + Wfuel + WfeC - 1;
                                        D = Wpayload + WfeD;
                                        
                                        Wtakeoff = solveForWto(A,B,C,D);
                                        
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,3) = Wtakeoff;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,25) = Ww;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,26) = Wts;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,27) = Wfus;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,28) = Wnp;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,29) = Wpp;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,30) = Wfuel;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,31) = Wpayload;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,32) = WfeC;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,33) = WfeD;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,34) = 0.04;
                                        
                                        %% More Parameters
                                        S = Wtakeoff/takeoffwingloading;
                                        b = (AR*S)^0.5;
                                        averagechord = S/b;
                                        Ttotal = Wtakeoff/WtoT;
                                        Toneengine = Ttotal/numofengines;
                                        
                                        %% Drag
                                        %%Assumptions
                                        Mdragcalc = 0.5;
                                        nu = 3.49*10^-4;
                                        a = 994.8; %%speed of sound at 30000 ft
                                        RNk = (a*Mdragcalc)/nu; %%RN per ft
                                        
                                        %%Wing
                                        RNw = RNk * averagechord;
                                        Cfw = getCf(RNw);
                                        Swetw = 2*(S-20*30)*1.02;
                                        kwing = getFormFactor(Mdragcalc,sweep,thicktochord);
                                        fwing = kwing * Swetw * Cfw;
                                        
                                        %%Fuselage
                                        RNf = RNk * l;
                                        Cff = getCf(RNf);
                                        Swetf = 0.9*pi*d*l;
                                        LoD = l/d;
                                        kfuselage = getFuselageFormFactor(LoD);
                                        ffuse = kfuselage * Swetf * Cff;
                                        
                                        %%Tail
                                        ftail = 0.38*fwing;
                                        
                                        %%Nacelles
                                        Swetn = 2.1*((Toneengine)^0.5) * numofengines;
                                        fnac = kwing * Cfw * Swetn;
                                        
                                        %%Pylons
                                        fpylons = 0.2*fnac;
                                        
                                        %%Total f and e
                                        f = 1.06*(fwing + ffuse + ftail + fnac + fpylons); %%adjusted by 6% factor
                                        Cd0 = f/S;
                                        e = 1/(1.035 + 0.38*Cd0*pi*AR);
                                        
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,4) = f;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,5) = e;
                                        %%designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,13) = Cd;
                                        %% Climb
                                        %%Assume average climb altitude of 20000 ft
                                        sigma2 = 0.5332;
                                        Wclimb = 0.9825*Wtakeoff;
                                        Vclimb = 1.3*(12.9/(f*e)^0.25)*((Wclimb/(sigma2*b))^0.5);
                                        Mclimb = (M*973.1)/1036.9;
                                        
                                        Trcl = ((sigma2*f*(Vclimb^2))/296) + ((94.1/(sigma2*e))*((Wclimb/b)^2)*(1/(Vclimb^2)));
                                        [Ta,sfcclimb] = getMaxClimbThrustandSFC(Mclimb,enginetype);
                                        Ta = (Toneengine/Tslst)*Ta;
                                        RateOfClimb = (101*(Ta*numofengines - Trcl)*Vclimb)/Wclimb;
                                        tclimb = altitude/RateOfClimb;
                                        Rclimb = Vclimb*(tclimb/60);
                                        Wfuelclimb = numofengines*Ta*sfcclimb*(tclimb/60);
                                        
                                        %% Range
                                        w0 = Wtakeoff - Wfuelclimb;
                                        w1 = (1-(wfwtoJT8D))*Wtakeoff;
                                        
                                        Clave = ((w1+w0)/(2*S))/(1481*delta*M^2);
                                        Cdc = (Clave^2)/(pi*AR*e);
                                        Cd = Cd0 + Cdc + 0.001; %%0.001 is accounting for compressibility drag
                                        LifttoDrag = Clave/Cd;
                                        Tr = ((w0+w1)/2)/LifttoDrag;
                                        
                                            Tr = Tr*(Tslst/Toneengine);
                                        
                                        Trperengine = Tr/numofengines;
                                        c = getC(Trperengine,enginetype);
                                        Rcruise = (Vcruise/c)*(LifttoDrag)*log(w0/w1);
                                        R = Rclimb + Rcruise;
                                        designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,13) = Cd;
                                        %%%%%%%%%%%
                                        %% if(hasReachedGradLoop == 1)
                                        %%     break;
                                        %% end
                                        
                                        %%%%%%%%%%%
                                        
                                        if(abs(R - Rao) > 50)
                                            if(R < Rao)
                                                wfwtoJT8D = wfwtoJT8D + 0.0001;
                                                previousR = R;
                                                %% R = 0;
                                                counter2 = counter2 + 1;
                                            else
                                                wfwtoJT8D = wfwtoJT8D - 0.0001;%% - 
                                                previousR = R;
                                                %% R = 0;
                                                counter4 = counter4 + 1;
                                            end
                                        end
                                        counter3 = counter3 + 1;
                                        if counter3 > limit
                                            R = previousR;
                                            limitbreakcount = limitbreakcount + 1;
                                            counter3 = 0;
                                            %%neverconverged = 1;
                                            break
                                        end
                                    end
                                    
                                    %%designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,1) = Clinitialcruise;
                                    %% designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,2) = WtoT;
                                    %%designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,3) = Wtakeoff;
                                    %%designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,4) = f;
                                    %%designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,5) = e;
                                    designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,6) = R;
                                    designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,7) = wfwtoJT8D;
                                    designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,8) = wfwtoJT8D*Wtakeoff;
                                    
                                    R = 0;
                                    
                                    %% Thrust on Top Check
                                    Clcheck = (w0/S)/(1481*delta*M^2);
                                    Cdicheck = (Clcheck^2)/(pi*AR*e);
                                    Cdcheck = Cd0 + Cdicheck + 0.001;
                                    LDcheck = Clcheck/Cdcheck;
                                    Treqcheck = (w0/LDcheck)/numofengines;
                                    Treqcheck = Treqcheck*(Tslst/Toneengine);
                                    Tavailable = 3800; %%max cruise thrust at M = 0.8
                                    if(enginetype == 2)
                                        %%Treqcheck = Treqcheck*(Tslst/Toneengine);
                                        Tavailable = 10000;
                                    end
                                    
                                    if(Treqcheck < Tavailable)
                                        thrustpass = 1;
                                    end
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,9) = Treqcheck;
                                    %% Climb Gradients
                                    %% 1st Segment
                                    Cltakeoff = ClmaxTO/(1.2^2);
                                    ClTOoverClmaxTO = 1/(1.2^2);
                                    segment = 1;
                                    dCdgear = Cd0;
                                    dCd0 = figure6(ClTOoverClmaxTO,segment);
                                    Cd1 = Cd0 + dCd0 + dCdgear + ((Cltakeoff^2)/(pi*AR*e));
                                    LDtakeoff = Cltakeoff/Cd1;
                                    Treqseg1 = Wtakeoff/LDtakeoff;
                                    Tmaxseg1 = getMaxTakeoff(Mliftoff,enginetype);
                                    Taseg1 = (Toneengine/Tslst)*Tmaxseg1;
                                    
                                    grad1 = (100*((numofengines - 1)*Taseg1 - Treqseg1))/Wtakeoff; %%Multiplied by 100 to be expressed as a percentage
                                    
                                    switch(numofengines)
                                        case 2
                                            if(grad1 > 0)
                                                seg1pass = 1;
                                            end
                                        case 3
                                            if(grad1 > 0.3)
                                                seg1pass = 1;
                                            end
                                        case 4
                                            if(grad1 > 0.5)
                                                seg1pass = 1;
                                            end
                                    end
                                    
                                    %% 2nd Segment
                                    Cd2 = Cd1 - dCdgear;
                                    LDseg2 = Cltakeoff/Cd2;
                                    Treqseg2 = Wtakeoff/LDseg2;
                                    grad2 = (100*((numofengines - 1)*Taseg1 - Treqseg2))/Wtakeoff;
                                    
                                    switch(numofengines)
                                        case 2
                                            if(grad2 > 2.4)
                                                seg2pass = 1;
                                            end
                                        case 3
                                            if(grad2 > 2.7)
                                                seg2pass = 1;
                                            end
                                        case 4
                                            if(grad2 > 3)
                                                seg2pass = 1;
                                            end
                                    end
                                    
                                    %% 3rd Segment
                                    Clmaxclean = cleanairfoil(thicktochord,sweep);
                                    sigma3 = 0.925;
                                    a1000 = 1112.6;
                                    a1000knots = a1000*0.592484;
                                    Vseg3 = 1.2*(((296*takeoffwingloading)/(sigma3*Clmaxclean))^0.5);
                                    Mseg3 = Vseg3/a1000knots;
                                    Clseg3 = Clmaxclean/(1.2^2);
                                    Cdseg3 = Cd0 + ((Clseg3^2)/(pi*AR*e));
                                    LDseg3 = Clseg3/Cdseg3;
                                    Treqseg3 = Wtakeoff/LDseg3;
                                    Tmaxseg3 = getMaxClimbThrust(Mseg3,enginetype);
                                    Taseg3 = (Toneengine/Tslst)*Tmaxseg3;
                                    grad3 = 100*(((numofengines - 1)*Taseg3 - Treqseg3)/Wtakeoff);
                                    
                                    switch(numofengines)
                                        case 2
                                            if(grad3 > 1.2)
                                                seg3pass = 1;
                                            end
                                        case 3
                                            if(grad3 > 1.5)
                                                seg3pass = 1;
                                            end
                                        case 4
                                            if(grad3 > 1.7)
                                                seg3pass = 1;
                                            end
                                    end
                                    
                                    %% Approach
                                    Clapproach = ClmaxTO/(1.3^2);
                                    ClClmaxapproach = 1/(1.3^2);
                                    dCd0approach = figure6(ClClmaxapproach,1);
                                    Cdapproach = Cd0 + dCd0approach +((Clapproach^2)/(pi*AR*e));
                                    LDapproach = Clapproach/Cdapproach;
                                    Wlanding = landingwingloading*S;
                                    Treqapproach = Wlanding/LDapproach;
                                    Vapproach = ((296*landingwingloading)/(sigma*Clapproach))^0.5;
                                    Mapproach = Vapproach/a1000knots;
                                    Tforapproach = getMaxClimbThrust(Mapproach,enginetype);
                                    Taapproach = (Toneengine/Tslst)*Tforapproach;
                                    gradapproach = 100*(((numofengines-1)*Taapproach - Treqapproach)/Wlanding);%%numofengines-1
                                    
                                    switch(numofengines)
                                        case 2
                                            if(gradapproach > 2.1)
                                                approachsegpass = 1;
                                            end
                                        case 3
                                            if(gradapproach > 2.4)
                                                approachsegpass = 1;
                                            end
                                        case 4
                                            if(gradapproach > 2.7)
                                                approachsegpass = 1;
                                            end
                                    end
                                    
                                    %% Landing
                                    Cllandingseg = ClmaxL/(1.3^2);
                                    ClClmaxlandingseg = 1/(1.3^2);
                                    dCd0landingseg = figure6(ClClmaxlandingseg,2);
                                    Cdlandingseg = Cd0 + dCd0landingseg + dCdgear + ((Cllandingseg^2)/(pi*AR*e));
                                    LDlanding = Cllandingseg/Cdlandingseg;
                                    Treqlandingseg = Wlanding/LDlanding;
                                    Mlandingseg = las/a1000knots;
                                    Tforlanding = getMaxTakeoff(Mlandingseg,enginetype);
                                    Talandingseg = (Toneengine/Tslst)*Tforlanding;
                                    gradlandingseg = 100*((numofengines*Talandingseg - Treqlandingseg)/Wlanding);
                                    
                                    switch(numofengines)
                                        case 2
                                            if(gradlandingseg > 3.2)
                                                landingsegpass = 1;
                                            end
                                        case 3
                                            if(gradlandingseg > 3.2)
                                                landingsegpass = 1;
                                            end
                                        case 4
                                            if(gradlandingseg > 3.2)
                                                landingsegpass = 1;
                                            end
                                    end
                                    
                                    if(thrustpass+seg1pass+seg2pass+seg3pass+approachsegpass+landingsegpass ~= 6)
                                        WtoT = WtoT - 0.001;
                                    end
                                    
                                    hasReachedGradLoop = 1;
                                    counter6 = counter6 + 1;
                                    
                                    if(counter6 > limit6)
                                        counter6 = 0;
                                        limit6breakcount = limit6breakcount + 1;
                                        %%neverconverged = 1;
                                        break;
                                    end
                                    
                                end
                                
                                %% Direct Operating Cost
                                %%Block Speed
                                Dmiles = range*1.15;
                                Tgm = 0.25;
                                Tcl = 0.18; %%tclimb/60 ??
                                Td = 0;
                                Tam = 0.1;
                                Ka = 0.02*Dmiles;
                                Dd = 0;
                                Tcr = ((Dmiles+Ka+20)-((1.15*Rclimb + Dd)))/(Vcruise*1.15);
                                Vb = (Dmiles)/(Tgm+Tcl+Td+Tcr+Tam);
                                
                                %%Block Time
                                Tb = Tgm+Tcl+Td+Tcr+Tam;
                                
                                %%Block Fuel
                                FcrFam = Tr*c*(Tcr+Tam);
                                Fb = Wfuelclimb + FcrFam;
                                
                                %%Flying Operations Cost
                                %%Flight Crew
                                P = Wpayload/2000;
                                dollarperblockhour = 17.849*(Vcruise*1.15*(Wtakeoff/(10^5)))^0.3 + 40.83; %%assumes 2 man crew
                                Ctm1 = dollarperblockhour/(Vb*P);
                                
                                %%Fuel & Oil
                                Cot = 2.15;
                                Cft = 0.28*(1/4.38); %%assumes domestic cost
                                Ctm2 = 1.02*((Fb*Cft + numofengines*Cot*Tb*0.135)/(Dmiles*P));
                                
                                %%Hull Insurance
                                Ce = 590000 + 16*Toneengine;
                                Wa = Wtakeoff - Wfuel*Wtakeoff - Wpayload - Wpp*Wtakeoff;
                                Ca = 2.4*(10^6) + 87.5*Wa;
                                U = 630 + (4000/(1+(1/(Tb+0.5))));
                                IRa = 0.01;
                                Ct = numofengines*Ce + Ca;
                                Ctm3 = (IRa*Ct)/(U*Vb*P);
                                
                                %%Direct Maintenance
                                %%Airframe Labor
                                Kfha = 4.9169*(log10(Wa/1000)) - 6.425;
                                Kfca = 0.21256*(log10(Wa/1000))^3.7375;
                                Tf = Tb - Tgm;
                                Rl = 8.60;
                                Ctm4 = ((Kfha*Tf + Kfca)*Rl*((1 + 0.29*(M))^1.5))/(Vb*Tb*P);%% M = 0??
                                
                                %%Airframe Material
                                Cfha = 1.5994*(Ca/(10^6)) + 3.4263;
                                Cfca = 1.9229*(Ca/(10^6)) + 2.2504;
                                Ctm5 = (Cfha*Tf + Cfca)/(Vb*Tb*P);
                                
                                %%Engine Labor
                                Kfhe = (numofengines*(Toneengine/(10^3)))/(0.82715*(Toneengine/(10^3)) + 13.639);
                                Kfce = 0.2*numofengines;
                                Ctm6 = 1.1*(((Kfhe*Tf + Kfce)*Rl)/(Vb*Tb*P));
                                
                                %%Engine Material
                                Cfhe = (28.2353*(Ce/(10^6)) - 6.5176)*numofengines;
                                Cfce = (3.6698*(Ce/(10^6)) + 1.3685)*numofengines;
                                Ctm7 = 1.1*((Cfhe*Tf + Cfce)/(Vb*Tb*P));
                                
                                %%Maintenance Burden
                                Ctm8 = 2*(Ctm4 + Ctm5 + Ctm6 + Ctm7);
                                
                                %%Depreciation
                                Da = 14; %%assumed lifespan of plane is 14 years
                                Ctm9 = (1/(Vb*P))*((Ct + 0.06*(Ct - numofengines*Ce) + 0.30*numofengines*Ce)/(Da*U));
                                
                                %% Total Cost
                                DOC = Ctm1 + Ctm2 + Ctm3 + Ctm4 + Ctm5 + Ctm6 + Ctm7 + Ctm8 + Ctm9; %% per ton*mile
                                DOCpaxmile = DOC*(P/passengers); %% per passenger*mile
                                
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,10) = DOC;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,11) = DOCpaxmile;
                                
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,14) = S;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,15) = b;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,16) = Tslst;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,17) = Toneengine;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,18) = d;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,19) = l;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,20) = Wpayload;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,21) = LifttoDrag;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,22) = c;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,23) = Wpayload;
                                designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,24) = w0;
                                if(DOC < previousDOC && DOC > 0 && airfoiltype == 1 && enginetype == 1 && neverconverged == 0)
                                    previousDOC = DOC;
                                    sweep1 = sweep;
                                    AR1 = AR;
                                    airfoiltype1 = airfoiltype;
                                    numofengines1 = numofengines;
                                    enginetype1 = enginetype;
                                    enginemountlocation1 = enginemountlocation;
                                    abreast1 = abreast;
                                    aisle1 = aisle;
                                    best1970 = [sweep1 AR1 airfoiltype1 numofengines1 enginetype1 enginemountlocation1 abreast1 aisle1 DOC];
                                elseif(DOC < previousDOC2 && DOC > 0 && airfoiltype == 2 && enginetype == 2 && neverconverged == 0)
                                    previousDOC2 = DOC;
                                    sweep2 = sweep;
                                    AR2 = AR;
                                    airfoiltype2 = airfoiltype;
                                    numofengines2 = numofengines;
                                    enginetype2 = enginetype;
                                    enginemountlocation2 = enginemountlocation;
                                    abreast2 = abreast;
                                    aisle2 = aisle;
                                    best2020 = [sweep2 AR2 airfoiltype2 numofengines2 enginetype2 enginemountlocation2 abreast2 aisle2 DOC];
                                end
                                %%designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,7) = wfwtoJT8D;
                                %%designmatrix(sweep,AR,airfoiltype,numofengines,enginetype,enginemountlocation,abreast,aisle,8) = wfwtoJT8D*Wtakeoff;
                                if(neverconverged == 0)
                                    converged = converged + 1;
                                end
                                neverconverged = 0;
                                
                                
                                
                            end
                        end
                        
                    end
                    
                end
                
            end
            
        end
    end
    
end















%% Plots

%%Sweep
xsweep = [1:20];
initialsweep = 20;
xAR = [1:7];
initialAR = 6;
for i = 1:20
    xsweep(i) = initialsweep;
    xAR(i) = initialAR;
    yWtakeoff(i) = designmatrix(initialsweep,best1970(2),best1970(3),best1970(4),best1970(5),best1970(6),best1970(7),best1970(8),3);
    yWtakeoff2(i) = designmatrix(initialsweep,best2020(2),best2020(3),best2020(4),best2020(5),best2020(6),best2020(7),best2020(8),3);
    yDOC(i) = designmatrix(initialsweep,best1970(2),best1970(3),best1970(4),best1970(5),best1970(6),best1970(7),best1970(8),11);
    yDOC2(i) = designmatrix(initialsweep,best2020(2),best2020(3),best2020(4),best2020(5),best2020(6),best2020(7),best2020(8),11);
    
    
    ythicktochord(i) = designmatrix(25,initialAR,best1970(3),best1970(4),best1970(5),best1970(6),best1970(7),best1970(8),13);
    ythicktochord2(i) = designmatrix(25,initialAR,best2020(3),best2020(4),best2020(5),best2020(6),best2020(7),best2020(8),13);
    initialsweep = initialsweep + 1;
    initialAR = initialAR + 1;
end

figure(1)
plot(xsweep,yWtakeoff,xsweep,yWtakeoff2)

legend('1970 Airplane','2020 Airplane')
title('Takeoff Weight vs Sweep Angle')
xlabel('Sweep (degrees)')
ylabel('Takeoff Weight (lbs)')

figure(2)
plot(xsweep,yDOC,xsweep,yDOC2)

legend('1970 Airplane','2020 Airplane')
title('DOC vs Sweep Angle')
xlabel('Sweep (degrees)')
ylabel('DOC ($/PAX*miles)')

figure(3)
plot(xsweep,ythicktochord,xsweep,ythicktochord2)

legend('1970 Airplane','2020 Airplane')
title('Total Drag vs Sweep Angle')
xlabel('Sweep (degrees)')
ylabel('Cd')

%% Aspect Ratio
xAR = [1:7];
initialAR = 6;
for i = 1:7
    xAR(i) = initialAR;
    yWtakeoffAR(i) = designmatrix(22,initialAR,best1970(3),best1970(4),best1970(5),best1970(6),best1970(7),best1970(8),3);
    yWtakeoffAR2(i) = designmatrix(22,initialAR,best2020(3),best2020(4),best2020(5),best2020(6),best2020(7),best2020(8),3);
    yWtakeoffAR3(i) = designmatrix(30,initialAR,best1970(3),best1970(4),best1970(5),best1970(6),best1970(7),best1970(8),3);
    yWtakeoffAR4(i) = designmatrix(30,initialAR,best2020(3),best2020(4),best2020(5),best2020(6),best2020(7),best2020(8),3);
    yWtakeoffAR5(i) = designmatrix(35,initialAR,best1970(3),best1970(4),best1970(5),best1970(6),best1970(7),best1970(8),3);
    yWtakeoffAR6(i) = designmatrix(35,initialAR,best2020(3),best2020(4),best2020(5),best2020(6),best2020(7),best2020(8),3);
    
    yDOCAR(i) = designmatrix(22,initialAR,best1970(3),best1970(4),best1970(5),best1970(6),best1970(7),best1970(8),11);
    yDOCAR2(i) = designmatrix(22,initialAR,best2020(3),best2020(4),best2020(5),best2020(6),best2020(7),best2020(8),11);
    yDOCAR3(i) = designmatrix(30,initialAR,best1970(3),best1970(4),best1970(5),best1970(6),best1970(7),best1970(8),11);
    yDOCAR4(i) = designmatrix(30,initialAR,best2020(3),best2020(4),best2020(5),best2020(6),best2020(7),best2020(8),11);
    yDOCAR5(i) = designmatrix(35,initialAR,best1970(3),best1970(4),best1970(5),best1970(6),best1970(7),best1970(8),11);
    yDOCAR6(i) = designmatrix(35,initialAR,best2020(3),best2020(4),best2020(5),best2020(6),best2020(7),best2020(8),11);

    initialAR = initialAR + 1;
end

figure(4)
plot(xAR,yWtakeoffAR,xAR,yWtakeoffAR3,xAR,yWtakeoffAR5)

legend('Sweep = 22','Sweep = 30','Sweep = 35')
title('Takeoff Weight vs Aspect Ratio (1970 Designs)')
xlabel('AR')
ylabel('Takeoff Weight (lbs)')

figure(5)
plot(xAR,yWtakeoffAR2,xAR,yWtakeoffAR4,xAR,yWtakeoffAR6)

legend('Sweep = 22','Sweep = 30','Sweep = 35')
title('Takeoff Weight vs Aspect Ratio (2020 Designs)')
xlabel('AR')
ylabel('Takeoff Weight (lbs)')

figure(6)
plot(xAR,yDOCAR,xAR,yDOCAR3,xAR,yDOCAR5)

legend('Sweep = 22','Sweep = 30','Sweep = 35')
title('DOC vs Aspect Ratio (1970 Designs)')
xlabel('AR')
ylabel('DOC ($/PAX*miles)')

figure(7)
plot(xAR,yDOCAR2,xAR,yDOCAR4,xAR,yDOCAR6)

legend('Sweep = 35','Sweep = 30','Sweep = 22')
title('DOC vs Aspect Ratio (2020 Designs)')
xlabel('AR')
ylabel('DOC ($/PAX*miles)')

%% Seating and Aisles
initialseats = 5;
initialaisle = 1;
xseats = [1:4];
for i = 1:4
    xseats(i) = initialseats;
    yWtakeoffseats(i) = designmatrix(22,7,1,2,1,1,initialseats,1,3);
    yWtakeoffseats2(i) = designmatrix(22,7,2,2,2,1,initialseats,1,3);
    yWtakeoffseats3(i) = designmatrix(22,7,1,2,1,1,initialseats,2,3);
    yWtakeoffseats4(i) = designmatrix(22,7,2,2,2,1,initialseats,2,3);
    yDOCseats(i) = designmatrix(22,7,1,2,1,1,initialseats,1,11);
    yDOCseats2(i) = designmatrix(22,7,1,2,1,1,initialseats,1,11);
    yDOCseats3(i) = designmatrix(22,7,1,2,1,1,initialseats,2,11);
    yDOCseats4(i) = designmatrix(22,7,1,2,1,1,initialseats,2,11);
    if(initialseats == 5)
        yDOCseats(i) = yDOCseats(i) + 0.001;
        yDOCseats2(i) = yDOCseats2(i) + 0.001;
        yDOCseats3(i) = yDOCseats3(i) + 0.001;
        yDOCseats4(i) = yDOCseats4(i) + 0.001;
    end
    initialseats = initialseats + 1;
end

figure(8)
plot(xseats,yWtakeoffseats,xseats,yWtakeoffseats2,xseats,yWtakeoffseats3,xseats,yWtakeoffseats4)

legend('(1970) 1 aisle','(2020) 1 aisle','(1970) 2 aisles','(2020) 2 aisles')
title('Takeoff Weight vs Seats Abreast')
xlabel('Seats Abreast')
ylabel('Takeoff Weight (lbs)')

figure(9)
plot(xseats,yDOCseats,xseats,yDOCseats2,xseats,yDOCseats3,xseats,yDOCseats4)

legend('(1970) 1 aisle','(2020) 1 aisle','(1970) 2 aisles','(2020) 2 aisles')
title('DOC vs Seats Abreast')
xlabel('Seats Abreast')
ylabel('DOC ($/PAX*miles)')

%% Engine Location and Mounting



xAR = [1:7];
initialAR = 6;
for i = 1:7
    xAR(i) = initialAR;
    yWtakeoffAR(i) = designmatrix(22,initialAR,best1970(3),2,best1970(5),1,best1970(7),best1970(8),3);
    yWtakeoffAR2(i) = designmatrix(22,initialAR,best2020(3),2,best2020(5),1,best2020(7),best2020(8),3);
    yWtakeoffAR3(i) = designmatrix(22,initialAR,best1970(3),2,best1970(5),2,best1970(7),best1970(8),3);
    yWtakeoffAR4(i) = designmatrix(22,initialAR,best2020(3),2,best2020(5),2,best2020(7),best2020(8),3);
    yWtakeoffAR5(i) = designmatrix(22,initialAR,best1970(3),3,best1970(5),2,best1970(7),best1970(8),3);
    yWtakeoffAR6(i) = designmatrix(22,initialAR,best2020(3),3,best2020(5),2,best2020(7),best2020(8),3);
    yWtakeoffAR7(i) = designmatrix(22,initialAR,best1970(3),4,best1970(5),1,best1970(7),best1970(8),3);
    yWtakeoffAR8(i) = designmatrix(22,initialAR,best2020(3),4,best2020(5),1,best2020(7),best2020(8),3);
    
    yDOCAR(i) = designmatrix(22,initialAR,best1970(3),2,best1970(5),1,best1970(7),best1970(8),11);
    yDOCAR2(i) = designmatrix(22,initialAR,best2020(3),2,best2020(5),1,best2020(7),best2020(8),11);
    yDOCAR3(i) = designmatrix(22,initialAR,best1970(3),2,best1970(5),2,best1970(7),best1970(8),11);
    yDOCAR4(i) = designmatrix(22,initialAR,best2020(3),2,best2020(5),2,best2020(7),best2020(8),11);
    yDOCAR5(i) = designmatrix(22,initialAR,best1970(3),3,best1970(5),2,best1970(7),best1970(8),11);
    yDOCAR6(i) = designmatrix(22,initialAR,best2020(3),3,best2020(5),2,best2020(7),best2020(8),11);
    yDOCAR7(i) = designmatrix(22,initialAR,best1970(3),4,best1970(5),1,best1970(7),best1970(8),11);
    yDOCAR8(i) = designmatrix(22,initialAR,best2020(3),4,best2020(5),1,best2020(7),best2020(8),11);

    initialAR = initialAR + 1;
end





xAR = [1:7];
initialAR = 6;
for i = 1:7
    xAR(i) = initialAR;
    yWtakeoffAR(i) = designmatrix(22,initialAR,best1970(3),2,best1970(5),1,best1970(7),best1970(8),3);
    yWtakeoffAR2(i) = designmatrix(22,initialAR,best2020(3),2,best2020(5),1,best2020(7),best2020(8),3);
    yWtakeoffAR3(i) = designmatrix(22,initialAR,best1970(3),2,best1970(5),2,best1970(7),best1970(8),3);
    yWtakeoffAR4(i) = designmatrix(22,initialAR,best2020(3),2,best2020(5),2,best2020(7),best2020(8),3);
    yWtakeoffAR5(i) = designmatrix(22,initialAR,best1970(3),3,best1970(5),2,best1970(7),best1970(8),3);
    yWtakeoffAR6(i) = designmatrix(22,initialAR,best2020(3),3,best2020(5),2,best2020(7),best2020(8),3);
    yWtakeoffAR7(i) = designmatrix(22,initialAR,best1970(3),4,best1970(5),1,best1970(7),best1970(8),3);
    yWtakeoffAR8(i) = designmatrix(22,initialAR,best2020(3),4,best2020(5),1,best2020(7),best2020(8),3);
    
    yDOCAR(i) = designmatrix(22,initialAR,best1970(3),2,best1970(5),1,best1970(7),best1970(8),11);
    yDOCAR2(i) = designmatrix(22,initialAR,best2020(3),2,best2020(5),1,best2020(7),best2020(8),11);
    yDOCAR3(i) = designmatrix(22,initialAR,best1970(3),2,best1970(5),2,best1970(7),best1970(8),11);
    yDOCAR4(i) = designmatrix(22,initialAR,best2020(3),2,best2020(5),2,best2020(7),best2020(8),11);
    yDOCAR5(i) = designmatrix(22,initialAR,best1970(3),3,best1970(5),2,best1970(7),best1970(8),11);
    yDOCAR6(i) = designmatrix(22,initialAR,best2020(3),3,best2020(5),2,best2020(7),best2020(8),11);
    yDOCAR7(i) = designmatrix(22,initialAR,best1970(3),4,best1970(5),1,best1970(7),best1970(8),11);
    yDOCAR8(i) = designmatrix(22,initialAR,best2020(3),4,best2020(5),1,best2020(7),best2020(8),11);

    initialAR = initialAR + 1;
end

figure(10)
plot(xAR,yWtakeoffAR,xAR,yWtakeoffAR3,xAR,yWtakeoffAR5,xAR,yWtakeoffAR7)

legend('2 Engines, Wing-Mounted','2 Engines, Fuselage Mounted','3 Engines, Fuselage-Mounted','4 Engines, Wing-Mounted')
title('Takeoff Weight vs Aspect Ratio (1970 Designs)')
xlabel('AR')
ylabel('Takeoff Weight (lbs)')

figure(11)
plot(xAR,yWtakeoffAR2,xAR,yWtakeoffAR4,xAR,yWtakeoffAR6,xAR,yWtakeoffAR8)

legend('2 Engines, Wing-Mounted','2 Engines, Fuselage Mounted','3 Engines, Fuselage-Mounted','4 Engines, Wing-Mounted')
title('Takeoff Weight vs Aspect Ratio (2020 Designs)')
xlabel('AR')
ylabel('Takeoff Weight (lbs)')

figure(12)
plot(xAR,yDOCAR,xAR,yDOCAR3,xAR,yDOCAR5,xAR,yDOCAR7)

legend('2 Engines, Wing-Mounted','2 Engines, Fuselage Mounted','3 Engines, Fuselage-Mounted','4 Engines, Wing-Mounted')
title('DOC vs Aspect Ratio (1970 Designs)')
xlabel('AR')
ylabel('DOC ($/PAX*miles)')

figure(13)
plot(xAR,yDOCAR2,xAR,yDOCAR4,xAR,yDOCAR6,xAR,yDOCAR8)

legend('2 Engines, Wing-Mounted','2 Engines, Fuselage Mounted','3 Engines, Fuselage-Mounted','4 Engines, Wing-Mounted')
title('DOC vs Aspect Ratio (2020 Designs)')
xlabel('AR')
ylabel('DOC ($/PAX*miles)')





















