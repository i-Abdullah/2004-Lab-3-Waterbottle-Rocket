%% info:

% this code runs Sensitivity Analysis for waterbottle rocket.
% It varies the parameters to achieve the optimal range.

%% initial parameters:

clear
clc
close all;


g = 9.81; % m/s2, acceleration due to gravity,
Cd= 0.8; % discharge coefficient
RhoAmb = 0.961; % kg/m^3 ambient air density
Volbottle= 0.002; % m^3 volume of empty bottle
Pamb= 12.1*6894.76; % converted to Pa, atmospheric pressure
GammaGas = 1.4; % ratio of specific heats for air
RhoWater = 1000; % kg/m^3, density of water
DThroat= 2.1; % cm, diameter of throat
DBottle= 10.5; % in cm, diameter of bottle
R = 287; %J/kgK, gas constant of air
MBottle= 0.15; % kg mass of empty 2-liter bottle with cone and fins
CD= 0.5; % drag coefficient
Pgage= 40*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
VWaterInit= 0.001; % m^3, initial volume of water inside bottle
TAirInit = 300; % K, initial temperature of
Airv0 = 0.0 ;% m/s, initial velocity of rocket
Theta= 45 ; % initial angle of rocket in degress
X0 = 0.0; % in meters, initial horizontal distance
y0 = 0.25; % in m, initial vertical height
TestStandLength= 0.5; % in m, length of test stand
VAirInit = Volbottle - VWaterInit ; %initial volume of Air.
ThroatArea = pi * ((DThroat*10^-2)/2)^2; %Area of throat
BottleArea  = pi * ((DBottle*10^-2)/2)^2; %Bottle Area
PayLoad = 25*10-3; % in kg
Fins = 10*10^-3 ; % in kg
TotalMass0 = Fins + PayLoad + MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air


VelX0 = 0;
VelZ0 = 0;
Range0 = 0;
Height0 = y0;


%% for-loops :( !

% i = change in angle
% k = change in water volume
% p = change in pressure (gage)
% s = change in coefficient of drag

maxi = 90;
maxk = 0.0015; % max
maxp = 50*6894.76 ; 
maxs = 1;
j = 1;

for i = 1:1:maxi
    
    for k = 0.0001:0.0001:maxk
        
        for p = 1*6894.76:1*6894.76:maxp
            
            
            for s = 0.01:0.03:maxs
                
                % re-define stuff that'll change again and again.
                
VAirInit = Volbottle - k ; %initial volume of Air.
TotalMass0 = Fins + PayLoad + MBottle + (k*RhoWater) + (((p+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((p+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air

                
                % call function
    [ Time Results ] = ode45(@(Time,States) RocketODE(Time,States,TestStandLength,i,p,Pamb,Cd,ThroatArea,s,BottleArea,RhoAmb,RhoWater,Volbottle,y0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R), [ 0 6],[TotalMass0 MassAirInit...
        VAirInit VelX0 VelZ0 Range0 y0 ]);

    
    % find maximum range by figuring out where the indexing changes.
    %find where y value is first negative
    h = 1;
            while Results(h,7)>0
                   h=h+1;
            end

%interpolate max distance
m=(Results(h,7)-Results(h-1,7))/(Results(h,6)-Results(h-1,6));
%concatinate the max distance with max distance for specific theta
maxR=Results(h-1,6)+(0-Results(h-1,7))/m;



                Range_results(1,j) = maxR;
                Conditions(j,:) = [ i k p s ] ;
                    
                j = j +1;
                
                
            end
        
        end
        
    end
    
end