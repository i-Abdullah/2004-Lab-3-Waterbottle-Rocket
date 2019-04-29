function [ derivatives ] = IspODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz)
% This's the ODE45 function for analyzing the water bottle rocket.
%
% For more info, read Assignment information.
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



%% Phase 3: 


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
