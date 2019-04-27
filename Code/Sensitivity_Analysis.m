%% info:

% this code runs Sensitivity Analysis for waterbottle rocket.
% It varies the parameters to achieve the optimal range.

%% initial parameters:

clear
clc
close all;


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
CD= 0.324  ; % drag coefficient
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
%MassAirInit = TotalMass0 - 600e-3 - MBottle  ;

%% velocity of the wind as initial conditions

% x y z respectively:

%x = downrange.
%z = height.
%y = cross range.

% the following from TA's Launch

[ xwind ywind ] = WindLaunch(204, 'NNE',4);
Vwx = xwind*0.44704 ; %multiplication to convert from mph to m/s
Vwy = ywind*0.44704 ;
Vwz = 0 ;



%% 

%initial conditions for ode:

VelX0 = 0;
VelZ0 = 0;
VelY0 = 0;

x0 = 0;
z0 = z0;
y0 = 0; %intial condition for location into the page

Opts = odeset('Events',@HitGround);


%% for-loops :( !

% i = change in angle
% k = change in water volume
% p = change in pressure (gage)

% s = change in temp % REMOVED!

maxi = 70;
maxk = 0.0011; % max
maxp = 40*6894.76 ; 

% counter
store = 1;

for i = 24:1:maxi
    
    for k = 1e-4:1e-6:maxk
        
        for p = 30*6894.76:2*6894.76:maxp
            
                            
                % re-define stuff that'll change again and again.
                
VAirInit = Volbottle - k ; %initial volume of Air.
TotalMass0 =  MBottle + (k*RhoWater) + (((p+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((p+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air

                
% Call ODE
[ Time Results ] = ode45(@(Time,States) ThermoODE(Time,States,TestStandLength,i,p,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 VelY0 x0 z0 y0 ],Opts);

        
Parameters(store,:) = [ i k p ];
RangeObtained (store,:) = max(Results(:,7));
                
store = store +1;
            end
        
        end
        
end


%% printout results;

index = find(max(RangeObtained)==RangeObtained);

fprintf('Your maximum range is %0.3f \n', max(RangeObtained));
fprintf('angle (deg): %0.3f \n', Parameters(index,1));
fprintf('Water Volume (m^3): %0.3f \n', Parameters(index,2));
fprintf('Gage pressure (Pa): %0.3f \n', Parameters(index,3));


