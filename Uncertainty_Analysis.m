%% info:


%{

this code will run uncertainty analysis using monter carlo simulations.

%}

%% housekeeping

clear
clc
close all;

%% ODE conditions


g = 9.81; % m/s2, acceleration due to gravity,
Cd= 0.8; % discharge coefficient
Rhoairamb = 0.961; % kg/m^3 ambient air density
Volbottle= 0.002; % m^3 volume of empty bottle
Pamb= 12.1*6894.76; % converted to Pa, atmospheric pressure
GammaGas = 1.4; % ratio of specific heats for air
RhoWater = 1000; % kg/m^3, density of water
DThroat= 2.1; % cm, diameter of throat
DBottle= 10.5; % in cm, diameter of bottle
R = 287; %J/kgK, gas constant of air
MBottle= 0.15; % kg mass of empty 2-liter bottle with cone and fins
MBottle= 0.110; % kg mass of empty 2-liter bottle with cone and fins
CD= 0.35; % drag coefficient
Pgage= 40*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
VWaterInit= 0.0006; % m^3, initial volume of water inside bottle
TAirInit = 300; % K, initial temperature of
TAirInit = 290.15; % K, initial temperature of
Airv0 = 0.0 ;% m/s, initial velocity of rocket
Theta= 45 ; % initial angle of rocket in degress
X0 = 0.0; % in meters, initial horizontal distance
z0 = 0.25; % in m, initial vertical height
TestStandLength= 0.5; % in m, length of test stand
VAirInit = Volbottle - VWaterInit ; %initial volume of Air.
ThroatArea = pi * ((DThroat*10^-2)/2)^2; %Area of throat
BottleArea  = pi * ((DBottle*10^-2)/2)^2; %Bottle Area
PayLoad = 0 ;
Fins = 0 ;
TotalMass0 = PayLoad + Fins + MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air

TotalMass0 = 694e-3 ;

% redefine the terms that has uncertintiy

N = 100;

CD = randn(N,1)*0.02 + CD ;
Pgage = randn(N,1)*2*6894.76 + Pgage ;
