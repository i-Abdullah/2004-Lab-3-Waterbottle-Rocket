function [ derivatives ] = ThrustODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz,Thrust_test,time_test,Mwater_i)
% This's the ODE45 function for analyzing the water bottle rocket.
%    Done by:
%            1- Brendan Palmer, id : 108102169
%            2- Abdulla AlAmeri id : 109364560
%
% For more info, read Project2.m .
% 
% --------- (Inputs)---------------------------
%   1- Time
%   2- States (explained below in order): the states of the rocket
%   3:end : the initial conditions, not assigned as globals thus has to be inputs
%
%
% --------- (Outputs)---------------------------
%
%   derivatives : The derivatives of the states.
%
% --------- (States In order (OLD!! CHECK NEW BELOW!!)------------------
% 1- Mass of rocket;
% 2- Mass of Air
% 3- Volume of Air;
% 4- Velocity x;
% 5- Velocity z;
% 6- downrange (X location);
% 7- Height (Z location);
% 8- Velocity y: 
% 9- Crossrange (Y location);
% 
%
%
%
% Purpose: Discover how to use derivatives and numerical integration to
% determine the rocket thrust and trajectory. Learn how rocket specifications
% and change trajectory.
% 
% Givens: Needed x displacement of the rocket. Derivatives modeling change in
% mass, pressure, volume, and velocity
% 
% Required: Determine the launch angle, , and bottle rocket dimensions to
% reach a particular distance.
% 
% Assumptions: Expansion of air is isentropic, no change of rocket deflection
% (in y- direction), rocket adjusts flight path to match the heading
% instantaneously, and static air around the rocket.
%
%
% --------- ( UPDATED ORDER OF STATES )------------------
% 1- Mass of rocket;
% 2- Mass of Air
% 3- Volume of Air;
% 4- Velocity x;
% 5- Velocity z;
% 6- Velocity y;
% 7- downrange (X location);
% 8- Height (Z location);
% 9- Crossrange (Y location);


%% define states:



RocketMass = States(1);
AirMass = States(2);
AirVolume = States(3);

dxdt = States(4);
dzdt = States(5);
dydt = States(6);

x = States(7);
z = States(8);
y = States(9);


%% Phase 1: 

if AirVolume < Volbottle
    

% Check if we still on stand or left:
if sqrt((x^2)+(z-z0)^2) <= TestStandLength
    %TotalVeloc = sqrt( (dzdt.^2) + (dxdt.^2) + (dydt.^2) );
    TotalVeloc = [ dxdt - Vwx ; dzdt - Vwz ; dydt - Vwy ] ;
    HeadingX = cosd(Theta);
    HeadingZ = sind(Theta);
    HeadingY = 0;
else
    TotalVeloc = [ dxdt - Vwx ; dzdt - Vwz ; dydt - Vwy ] ;
    HeadingX = TotalVeloc(1)/norm(TotalVeloc);
    HeadingZ = TotalVeloc(2)/norm(TotalVeloc);
    HeadingY = TotalVeloc(3)/norm(TotalVeloc);
end

Pressure = ( ( VAirInit ./ AirVolume ) .^ GammaGas ) .* (Pgage+Pamb) ; 
%Thrust = 2.* Cd .* ThroatArea .* ( Pressure - Pamb) ;

%find first min:
diff = abs(time_test-Time);
ind1 = find(diff==min(diff));

% zero and find the second min:
diff(ind1) = [];
ind2 = find(diff==min(diff));


Thrust = abs((Thrust_test(ind1) + Thrust_test(ind2) ))/2;


Drag = ( Rhoairamb / 2) .* (norm(TotalVeloc)).^2 * CD*BottleArea; 

% Define Derivatives

% Accelration in x and Z
dadt_X = ( (Thrust - Drag) * HeadingX) ./ RocketMass ;
dadt_Z =  ( ((Thrust - Drag) * HeadingZ) - RocketMass*g ) ./ RocketMass ;
dadt_Y = ( (Thrust - Drag) * HeadingY) ./ RocketMass ; 


%Volume (how volume changes with time)
DVolume_Dt = Cd * ThroatArea * sqrt ( (2/RhoWater) * ( ( (Pgage+Pamb) * (( VAirInit/AirVolume ) ^ (GammaGas))) - Pamb ));

%Mass (how mass of Rocket changes with time, it's )
DMass_Dt = - Cd .* ThroatArea .* sqrt ( 2.*RhoWater.* ( Pressure - Pamb ) );

%DMass_Dt = -(VWaterInit*RhoWater)


derivatives = [ DMass_Dt; 0; DVolume_Dt; dadt_X; dadt_Z; dadt_Y ; dxdt ; dzdt ; dydt ] ;
%% Phase 2:

elseif AirVolume>= Volbottle
%T and P of end states
Tend = TAirInit * (( VAirInit/Volbottle) ^ (GammaGas-1) );
Pend = (Pgage+Pamb) * (( VAirInit/Volbottle) ^ (GammaGas) );

    TotalVeloc = [ dxdt - Vwx ; dzdt - Vwz ; dydt - Vwy ] ;
    HeadingX = TotalVeloc(1)/norm(TotalVeloc);
    HeadingZ = TotalVeloc(2)/norm(TotalVeloc);
    HeadingY = TotalVeloc(3)/norm(TotalVeloc);
    
PressureCond = Pend * (AirMass/MassAirInit)^(GammaGas) ;

if PressureCond>Pamb
Density = AirMass / Volbottle;
Temp = PressureCond/(Density*R);
CriticalP = (PressureCond) * (2./(GammaGas+1)).^(GammaGas/(GammaGas-1));
  
if CriticalP > Pamb
    
    Mach  = 1;
    Texit = (2/(GammaGas+1))*Temp ;
    Vexit = sqrt(GammaGas*Texit*R);
    Pexit = CriticalP;
    Densityexit = CriticalP/(R*Texit) ;
    
elseif CriticalP <= Pamb
    
   Mach = sqrt(( (PressureCond/Pamb)^( ( (GammaGas-1)/GammaGas)) - 1 ) * (2/(GammaGas-1)));
   Texit = Temp/(1+((GammaGas-1)/2)*Mach^2);
   Pexit = Pamb;
   Densityexit = Pamb/(R*Texit) ;
   Vexit = Mach * sqrt(GammaGas*Texit*R);
end

% how mass of air and rocket change with time

MassAirFlowRate = Cd*Densityexit*ThroatArea*Vexit;
MassRocketFlowRate = -MassAirFlowRate;

Thrust = MassAirFlowRate *Vexit + (Pexit-Pamb)*ThroatArea ;
Drag = ( Rhoairamb / 2) .* (norm(TotalVeloc)).^2 * CD*BottleArea; 

dadt_X = ( (Thrust - Drag) * HeadingX) ./ RocketMass ;
dadt_Z =  ( ((Thrust - Drag) * HeadingZ) - RocketMass*g ) ./ RocketMass ;
dadt_Y = ( (Thrust - Drag) * HeadingY) ./ RocketMass ;



derivatives = [ MassRocketFlowRate; -MassAirFlowRate; 0; dadt_X; dadt_Z; dadt_Y; dxdt ; dzdt; dydt ] ;


%% Phase 3: 

else

    TotalVeloc = [ dxdt - Vwx ; dzdt - Vwz ; dydt - Vwy ] ;
    HeadingX = TotalVeloc(1)/norm(TotalVeloc);
    HeadingZ = TotalVeloc(2)/norm(TotalVeloc);
    HeadingY = TotalVeloc(3)/norm(TotalVeloc);

Thrust = 0 ;
Drag = ( Rhoairamb / 2) .* (norm(TotalVeloc)).^2 * CD*BottleArea; 

dadt_X = ( (Thrust - Drag) * HeadingX) ./ RocketMass ;
dadt_Z =  ( ((Thrust - Drag) * HeadingZ) - RocketMass*g ) ./ RocketMass ;
dadt_Y = ( (Thrust - Drag) * HeadingY) ./ RocketMass ;


derivatives = [ 0; 0; 0; dadt_X; dadt_Z; dadt_Y; dxdt ; dzdt ; dydt ] ;

end



end
