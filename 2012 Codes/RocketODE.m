function [ derivatives ] = RocketODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,y0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R)
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
% --------- (States In order)------------------
% 1- Mass of rocket;
% 2- Mass of Air
% 3- Volume of Air;
% 4- Velocity x;
% 5- Velocity z;
% 6- Range (X location);
% 7- Height (Z location);
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


% PREDEFINE Constants:

global t2 t1 t3

%% Phase 1: 

if States(3) < Volbottle
    

% Check if we still on stand or left:
if sqrt((States(6)^2)+(States(7)-y0)^2) <= TestStandLength
    TotalVeloc = sqrt( (States(5).^2) + (States(4).^2) );
    HeadingX = cosd(Theta);
    HeadingZ = sind(Theta);
else
    TotalVeloc = sqrt( (States(5).^2) + (States(4).^2) );
    HeadingX = States(4)/TotalVeloc;
    HeadingZ = States(5)/TotalVeloc;
end

Pressure = ( ( VAirInit ./ States(3) ) .^ GammaGas ) .* (Pgage+Pamb) ; 
Thrust = 2.* Cd .* ThroatArea .* ( Pressure - Pamb) ;
Drag = ( Rhoairamb / 2) .* (TotalVeloc).^2 * CD*BottleArea; 

% Define Derivatives

% Accelration in x and Z
dadt_X = ( (Thrust - Drag) * HeadingX) ./ States(1) ;
dadt_Z =  ( ((Thrust - Drag) * HeadingZ) - States(1)*g ) ./ States(1) ;

%Volume (how volume changes with time)
DVolume_Dt = Cd * ThroatArea * sqrt ( (2/RhoWater) * ( ( (Pgage+Pamb) * (( VAirInit/States(3) ) ^ (GammaGas))) - Pamb ));

%Mass (how mass of Rocket changes with time, it's )
DMass_Dt = - Cd .* ThroatArea .* sqrt ( 2.*RhoWater.* ( Pressure - Pamb ) );


derivatives = [ DMass_Dt; 0; DVolume_Dt; dadt_X; dadt_Z; States(4) ; States(5) ] ;
t1 = [ t1 ; [Time] ];
%% Phase 2:

elseif States(3)>= Volbottle
%T and P of end states
Tend = TAirInit * (( VAirInit/Volbottle) ^ (GammaGas-1) );
Pend = (Pgage+Pamb) * (( VAirInit/Volbottle) ^ (GammaGas) );

    TotalVeloc = sqrt( (States(5).^2) + (States(4).^2) );
    HeadingX = States(4)/TotalVeloc;
    HeadingZ = States(5)/TotalVeloc;
    
PressureCond = Pend * (States(2)/MassAirInit)^(GammaGas) ;

if PressureCond>Pamb
Density = States(2) / Volbottle;
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
Drag = ( Rhoairamb / 2) .* (TotalVeloc).^2 * CD*BottleArea; 

dadt_X = ( (Thrust - Drag) * HeadingX) ./ States(1) ;
dadt_Z =  ( ((Thrust - Drag) * HeadingZ) - States(1)*g ) ./ States(1) ;


derivatives = [ MassRocketFlowRate; -MassAirFlowRate; 0; dadt_X; dadt_Z; States(4) ; States(5) ] ;
t2 = [ t2 ; Time ];


%% Phase 3: 

else

TotalVeloc = sqrt( (States(5).^2) + (States(4).^2) );
HeadingX = States(4)/TotalVeloc;
HeadingZ = States(5)/TotalVeloc;

Thrust = 0 ;
Drag = ( Rhoairamb / 2) .* (TotalVeloc).^2 * CD*BottleArea; 

dadt_X = ( (Thrust - Drag) * HeadingX) ./ States(1) ;
dadt_Z =  ( ((Thrust - Drag) * HeadingZ) - States(1)*g ) ./ States(1) ;


derivatives = [ 0; 0; 0; dadt_X; dadt_Z; States(4) ; States(5) ] ;
t3 = [ t3 ; Time ];

end



end



