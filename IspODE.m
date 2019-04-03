function [ derivatives ] = RocketODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,y0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz)
% This's the ODE45 function for analyzing the water bottle rocket.
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
% 8- Velocity y: 
% 9- Cross Range: (Y locations);
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

if sqrt((States(6)^2)+(States(7)-y0)^2) <= TestStandLength
    %TotalVeloc = sqrt( (States(5).^2) + (States(4).^2) + (States(8).^2) );
    TotalVeloc = [ States(4) - Vwx ; States(5) - Vwz ; States(8) - Vwy ] ;
    HeadingX = cosd(Theta);
    HeadingZ = sind(Theta);
    HeadingY = 0;
else
    TotalVeloc = [ States(4) - Vwx ; States(5) - Vwz ; States(8) - Vwy ] ;
    HeadingX = TotalVeloc(1)/norm(TotalVeloc);
    HeadingZ = TotalVeloc(2)/norm(TotalVeloc);
    HeadingY = TotalVeloc(3)/norm(TotalVeloc);
end

Thrust = 0 ;
Drag = ( Rhoairamb / 2) .* (norm(TotalVeloc)).^2 * CD*BottleArea; 

dadt_X = ( (Thrust - Drag) * HeadingX) ./ States(1) ;
dadt_Z =  ( ((Thrust - Drag) * HeadingZ) - States(1)*g ) ./ States(1) ;
dadt_Y = ( (Thrust - Drag) * HeadingY) ./ States(1) ;


derivatives = [ 0; 0; 0; dadt_X; dadt_Z; States(4) ; States(5) ; dadt_Y ; States(8) ] ;

end
